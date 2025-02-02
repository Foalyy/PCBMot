from typing import Self
from enum import Enum
from config import Config
import svgwrite as svg
from geometry import sin, cos, tan, asin, acos, atan, atan2
from geometry import DrawableObject, Vector, Point, Line, Segment, Circle, Arc, PathSegment, PathArc, Path

class CoilSide(Enum):
    """Enum that describes the four sides of a coil: OUTER, RIGHT, INNER, and LEFT"""
    OUTER = 1
    RIGHT = 2
    INNER = 3
    LEFT = 4

class CoilConnection(Enum):
    """Enum that describes all the possible connection points for a coil, either via or terminal"""
    TERMINAL = 1
    OUTSIDE_OUTER_LEFT_VIA = 2
    OUTSIDE_OUTER_RIGHT_VIA = 3
    OUTSIDE_INNER_LEFT_VIA = 4
    OUTSIDE_INNER_RIGHT_VIA = 5
    INSIDE_OUTER_VIA = 6
    INSIDE_INNER_VIA = 7
    INSIDE_LEFT_VIA = 8
    INSIDE_RIGHT_VIA = 9

    def match_side(connection: Self, side: CoilSide) -> bool:
        """Return True if the given connection point could be connected to the given side, assuming a clockwise coil"""
        return \
            connection == CoilConnection.TERMINAL and side == CoilSide.OUTER or \
            connection == CoilConnection.OUTSIDE_OUTER_LEFT_VIA and side == CoilSide.OUTER or \
            connection == CoilConnection.OUTSIDE_OUTER_RIGHT_VIA and side == CoilSide.RIGHT or \
            connection == CoilConnection.OUTSIDE_INNER_LEFT_VIA and side == CoilSide.LEFT or \
            connection == CoilConnection.OUTSIDE_INNER_RIGHT_VIA and side == CoilSide.INNER or \
            connection == CoilConnection.INSIDE_OUTER_VIA and side == CoilSide.OUTER or \
            connection == CoilConnection.INSIDE_INNER_VIA and side == CoilSide.INNER or \
            connection == CoilConnection.INSIDE_LEFT_VIA and side == CoilSide.LEFT or \
            connection == CoilConnection.INSIDE_RIGHT_VIA and side == CoilSide.RIGHT
    
    def mirrored_y(connection: Self) -> Self:
        if connection == CoilConnection.OUTSIDE_OUTER_LEFT_VIA:
            return CoilConnection.OUTSIDE_OUTER_RIGHT_VIA
        elif connection == CoilConnection.OUTSIDE_OUTER_RIGHT_VIA:
            return CoilConnection.OUTSIDE_OUTER_LEFT_VIA
        elif connection == CoilConnection.OUTSIDE_INNER_LEFT_VIA:
            return CoilConnection.OUTSIDE_INNER_RIGHT_VIA
        elif connection == CoilConnection.OUTSIDE_INNER_RIGHT_VIA:
            return CoilConnection.OUTSIDE_INNER_LEFT_VIA
        elif connection == CoilConnection.INSIDE_OUTER_VIA:
            return CoilConnection.INSIDE_OUTER_VIA
        elif connection == CoilConnection.INSIDE_INNER_VIA:
            return CoilConnection.INSIDE_INNER_VIA
        elif connection == CoilConnection.INSIDE_LEFT_VIA:
            return CoilConnection.INSIDE_RIGHT_VIA
        elif connection == CoilConnection.INSIDE_RIGHT_VIA:
            return CoilConnection.INSIDE_LEFT_VIA

class Via:
    """A via allowing connections between layers on the board"""

    def __init__(self, center: Point, diameter: float, hole_diameter: float, tag=None):
        self.center: Point = center
        self.diameter: float = diameter
        self.hole_diameter: float = hole_diameter
        self.tag = tag
    
    def distance(self, object) -> float:
        """Calculate the distance between this Via and the given object

        Supported objects : see Point.distance().
        Throws a TypeError if any other object type is given.
        """
        return self.center.distance(object) - self.diameter / 2.0
    
    def rotated(self, rotation_center: Self, angle: float) -> Self:
        """Create a copy of this Via rotated around the given center point by the given angle"""
        return Via(self.center.rotated(rotation_center, angle), self.diameter, self.hole_diameter)
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Via mirrored about the X axis"""
        return Via(self.center.mirrored_x(), self.diameter, self.hole_diameter)
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Via mirrored about the Y axis"""
        return Via(self.center.mirrored_y(), self.diameter, self.hole_diameter)
    
    def draw_svg(self, drawing: svg.Drawing, via_color: str, via_hole_color: str):
        """Draw this Via on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        # Draw the pad
        drawing.add(drawing.circle(
            self.center.to_viewport().as_tuple(),
            self.diameter / 2.0,
            stroke = "none",
            fill = via_color,
            opacity = 1.0,
        ))

        # Draw the hole
        drawing.add(drawing.circle(
            self.center.to_viewport().as_tuple(),
            self.hole_diameter / 2.0,
            stroke = "none",
            fill = via_hole_color,
            opacity = 1.0,
        ))

        return self
    
    def closest_in_list(vias: list[Self], object) -> tuple[Self, float]:
        """Class method that returns the Via in the given list that is closest to the given object
        
        See distance() for the list of supported objects."""
        min_distance = 0.0
        closest = None
        for via in vias:
            distance = via.distance(object)
            if closest is None or distance < min_distance:
                min_distance = distance
                closest = via
        return (closest, min_distance)

class Coil:
    """A single coil on the board"""

    def __init__(self, path: Path, n_turns: int):
        self.path: Path = path
        self.n_turns: int = n_turns

    def generate(
            angle: float,
            outer_radius: float,
            inner_radius: float,
            anticlockwise: bool,
            trace_width: float,
            trace_spacing: float,
            outside_vias: dict[CoilConnection, Via],
            inside_vias: dict[CoilConnection, Via],
            outside_connection: CoilConnection,
            inside_connection: CoilConnection,
            max_turns: int,
            construction_geometry: list,
        ) -> Self:
        """Generate a coil centered around the vertical axis based on the given parameters"""

        # If the coil should be anticlockwise, generate a clockwise coil and mirror it along the Y axis
        # at the end. The requested connections should therefore be mirrored first in order to match.
        if anticlockwise:
            outside_connection = CoilConnection.mirrored_y(outside_connection)
            inside_connection = CoilConnection.mirrored_y(inside_connection)

        # Distance between the centerline of each adjacent coil turn
        loop_offset = trace_spacing + trace_width
        
        # Calculate the best outer and inner fillet radii of the outermost coil turn to fit the vias
        # Twice the trace width is a sane value to start for any track width
        outer_fillet_radius = trace_width * 2
        inner_fillet_radius = trace_width * 2
        outer_via = outside_vias[CoilConnection.OUTSIDE_OUTER_RIGHT_VIA]
        inner_via = outside_vias[CoilConnection.OUTSIDE_INNER_RIGHT_VIA]
        construction_line_right = Line.from_two_points(Point.origin(), Point.polar(angle/2.0, outer_radius)).offset(-loop_offset * 0.5)
        construction_arc_outer = Arc(Point.polar(-angle/2.0, outer_radius), Point.polar(angle/2.0, outer_radius), outer_radius)
        construction_arc_inner = Arc(Point.polar(-angle/2.0, inner_radius), Point.polar(angle/2.0, inner_radius), inner_radius)
        construction_point_outer = construction_line_right.intersect(construction_arc_outer)
        construction_point_inner = construction_line_right.intersect(construction_arc_inner)
        start_point = construction_arc_outer.midpoint()
        end_point = construction_arc_inner.midpoint()
        while outer_fillet_radius < outer_radius and inner_fillet_radius < inner_radius:
            path = Path(start_point)
            path.append_arc(construction_point_outer, outer_radius, anticlockwise=False)
            path.append_segment(construction_point_inner, fillet_radius=outer_fillet_radius)
            outer_fillet = Arc(path.elements[-3].p2, path.elements[-2].p2, radius=outer_fillet_radius)
            path.append_arc(end_point, inner_radius, anticlockwise=True, fillet_radius=inner_fillet_radius)
            inner_fillet = Arc(path.elements[-3].p2, path.elements[-2].p2, radius=inner_fillet_radius)
            modified = False
            if outer_via.distance(outer_fillet) <= trace_width / 2.0 + trace_spacing:
                outer_fillet_radius += 0.1
                modified = True
            if inner_via.distance(inner_fillet) <= trace_width / 2.0 + trace_width:
                inner_fillet_radius += 0.1
                modified = True
            if not modified:
                break

        outer_radius_initial = outer_radius
        inner_radius_initial = inner_radius
        line_left = Line.from_two_points(Point.origin(), Point.polar(-angle/2.0, outer_radius)).offset(loop_offset * 0.5 + outer_fillet_radius)
        arc_outer = Arc(Point.polar(-angle/2.0, outer_radius), Point.polar(angle/2.0, outer_radius), outer_radius)
        point_start = line_left.intersect(arc_outer)
        path = Path(point_start)
        n_turns = 0
        collision_via = None
        for i in range(max_turns or 1000):
            # Construction geometry
            line_left = Line.from_two_points(Point.origin(), Point.polar(-angle/2.0, outer_radius)).offset(loop_offset * (i + 0.5))
            line_right = Line.from_two_points(Point.origin(), Point.polar(angle/2.0, outer_radius)).offset(-loop_offset * (i + 0.5))
            arc_outer = Arc(Point.polar(-angle/2.0, outer_radius), Point.polar(angle/2.0, outer_radius), outer_radius)
            inner_radius = inner_radius_initial + i * loop_offset
            arc_inner = Arc(Point.polar(-angle/2.0, inner_radius), Point.polar(angle/2.0, inner_radius), inner_radius)

            # Arc to outer right
            point_outer_right = line_right.intersect(arc_outer)
            arc = Arc(path.end_point, point_outer_right, outer_radius)
            closest_via, distance = Via.closest_in_list(inside_vias.values(), arc)
            if distance < trace_width / 2.0 + trace_spacing: # Collision
                collision_via = closest_via
            path.append_arc(point_outer_right, outer_radius, anticlockwise=False, fillet_radius=outer_fillet_radius, tag=CoilSide.OUTER)
            if collision_via is not None:
                break

            # Segment to inner right
            point_inner_right = line_right.intersect(arc_inner)
            segment = Segment(path.end_point, point_inner_right)
            closest_via, distance = Via.closest_in_list(inside_vias.values(), segment)
            if distance < trace_width / 2.0 + trace_spacing: # Collision
                collision_via = closest_via
            path.append_segment(point_inner_right, fillet_radius=outer_fillet_radius, tag=CoilSide.RIGHT)
            if collision_via is not None:
                break

            # Arc to inner left
            point_inner_left = line_left.intersect(arc_inner)
            arc = Arc(path.end_point, point_inner_left, outer_radius, reverse=True)
            closest_via, distance = Via.closest_in_list(inside_vias.values(), arc)
            if distance < trace_width / 2.0 + trace_spacing: # Collision
                collision_via = closest_via
            path.append_arc(point_inner_left, inner_radius, anticlockwise=True, fillet_radius=inner_fillet_radius, tag=CoilSide.INNER)
            if collision_via is not None:
                break

            # Construction outer arc with offset
            outer_radius = outer_radius_initial - (i + 1) * loop_offset
            p1 = Point.polar(-angle/2.0, outer_radius)
            p2 = Point.polar(angle/2.0, outer_radius)
            arc_outer = Arc(p1, p2, outer_radius)

            # Segment to outer left with offset
            point_outer_left = line_left.intersect(arc_outer)
            segment = Segment(path.end_point, point_outer_left)
            closest_via, distance = Via.closest_in_list(inside_vias.values(), segment)
            if distance < trace_width / 2.0 + trace_spacing: # Collision
                collision_via = closest_via
            path.append_segment(point_outer_left, fillet_radius=inner_fillet_radius, tag=CoilSide.LEFT)
            if collision_via is not None:
                break

            # Reduce the fillet radius for the next loop
            outer_fillet_radius -= loop_offset
            if outer_fillet_radius < 0.1:
                outer_fillet_radius = 0.1
            inner_fillet_radius -= loop_offset
            if inner_fillet_radius < 0.1:
                inner_fillet_radius = 0.1
            
            # Count the number of turns in the coil
            n_turns += 1

        # Connect the start of the coil to the requested via or terminal
        if outside_connection is not None:
            if outside_connection == CoilConnection.TERMINAL:
                # TODO
                pass
            else:
                # Via to connect to
                try:
                    target_via = outside_vias[outside_connection]
                except KeyError:
                    raise ValueError("Invalid outside_connection")

                # Pop the first elements from the path until we find the one that should be connected to the target via
                while not outside_connection.match_side(path.first().tag):
                    path.pop_first()
                    path.pop_first() # Fillet
                    if len(path.elements) == 1:
                        raise ValueError("Cannot find path element to connect to requested outside connection")
                
                # Connect the start of the coil to the target via
                if outside_connection in [CoilConnection.OUTSIDE_OUTER_RIGHT_VIA, CoilConnection.OUTSIDE_INNER_LEFT_VIA]:
                    segment = Segment(path.first().p2, path.start_point)
                    tangent_arc = segment.tangent_arc_through_point(target_via.center)
                    path.prepend_arc(target_via.center, tangent_arc.radius, anticlockwise = segment.p2 == tangent_arc.p1)
                elif outside_connection in [CoilConnection.OUTSIDE_OUTER_LEFT_VIA, CoilConnection.OUTSIDE_INNER_RIGHT_VIA]:
                    first = path.first()
                    arc = Arc(path.first().p2, path.start_point, first.radius, reverse = not first.anticlockwise)
                    tangent_arc = arc.tangent_arc_through_point(target_via.center, at_start = not first.anticlockwise)
                    path.prepend_arc(target_via.center, tangent_arc.radius, anticlockwise = first.anticlockwise)

        # Connect the end of the coil to the requested via
        if inside_connection is not None:
            # Via to connect to
            try:
                target_via = inside_vias[inside_connection]
            except KeyError:
                raise ValueError("Invalid inside_connection")

            # Pop the last elements from the path until we find the one that should be connected to the target via
            while not inside_connection.match_side(path.last().tag):
                path.pop()
                if len(path.elements) == 1:
                    raise ValueError("Cannot find path element to connect to requested inside connection")
            
            # Replace the last element in the path with a corner connected to the target via
            projected = target_via.center.projected(path.last_geometry())
            replaced_element = path.pop()
            match replaced_element:
                case PathSegment():
                    path.append_segment(projected)
                    path.append_segment(target_via.center, fillet_radius=0.2)

                case PathArc():
                    path.append_arc(projected, replaced_element.radius, replaced_element.anticlockwise)
                    path.append_segment(target_via.center, fillet_radius=0.15)

        # Return the coil based on this path, mirrored relative to the Y axis if anticlockwise
        if anticlockwise:
            path = path.mirrored_y()
        return Coil(path, n_turns)

    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Coil rotated around the given center point by the given angle"""
        path = self.path.rotated(center, angle)
        return Coil(path, self.n_turns)
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Coil mirrored about the X axis"""
        path = self.path.mirrored_x()
        return Coil(path, self.n_turns)
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Coil mirrored about the Y axis"""
        path = self.path.mirrored_y()
        return Coil(path, self.n_turns)
    
    def draw_svg(self, drawing: svg.Drawing, color=None, opacity=None, thickness=None, dashes=None):
        """Draw this Path on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        self.path.draw_svg(drawing, color, opacity, thickness, dashes)
        return self

class PCB:
    """A PCB containing layers"""

    def __init__(self, config: Config, layers: dict, vias: list[Via], construction_geometry: list):
        self.config = config
        self.layers = layers
        self.vias = vias
        self.construction_geometry = construction_geometry
    
    def generate(config: Config):
        """Generate a new PCB based on the given config"""

        # Center the board on the origin
        board_center = Point.origin()

        # Board outline
        outline = [
            Circle(board_center, config.board_radius),
            Circle(board_center, config.hole_radius),
        ]

        # Construction geometry
        construction_geometry = [
            # Base coil center line
            svg.shapes.Line(
                board_center.to_viewport().as_tuple(), Point(0, config.board_radius).to_viewport().as_tuple(),
                stroke = config.construction_geometry_color,
                stroke_width = config.construction_geometry_thickness,
                stroke_dasharray=config.construction_geometry_dashes,
            ),
            # Base coil left line
            svg.shapes.Line(
                board_center.to_viewport().as_tuple(), Point.polar(-config.coil_angle/2.0, config.board_radius).to_viewport().as_tuple(),
                stroke = config.construction_geometry_color,
                stroke_width = config.construction_geometry_thickness,
                stroke_dasharray = config.construction_geometry_dashes,
            ),
            # Base coil right line
            svg.shapes.Line(
                board_center.to_viewport().as_tuple(), Point.polar(config.coil_angle/2.0, config.board_radius).to_viewport().as_tuple(),
                stroke = config.construction_geometry_color,
                stroke_width = config.construction_geometry_thickness,
                stroke_dasharray = config.construction_geometry_dashes,
            ),
            # Outer circle
            svg.shapes.Circle(
                board_center.to_viewport().as_tuple(), config.board_radius - config.board_outer_margin,
                stroke = config.construction_geometry_color,
                stroke_width = config.construction_geometry_thickness,
                stroke_dasharray = config.construction_geometry_dashes,
            ),
            # Inner circle
            svg.shapes.Circle(
                board_center.to_viewport().as_tuple(), config.hole_radius + config.board_inner_margin,
                stroke = config.construction_geometry_color,
                stroke_width = config.construction_geometry_thickness,
                stroke_dasharray = config.construction_geometry_dashes,
            ),
        ]

        # Inside vias
        # TODO : replace the "sep" parameter with a computation to make the line between via_inner_2 and via_inner_3 parallel to the side
        coil_B_center = Point(0, ((config.board_radius - config.board_outer_margin) + (config.hole_radius + config.board_inner_margin)) / 2.0)
        sep = 0.4
        inside_outer_via = Via(coil_B_center + Vector(0, (config.via_diameter_w_spacing + sep) / 2.0), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
        inside_inner_via = Via(coil_B_center - Vector(0, (config.via_diameter_w_spacing + sep) / 2.0), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
        c1 = Circle(inside_outer_via.center, config.via_diameter_w_spacing)
        c2 = Circle(inside_inner_via.center, config.via_diameter_w_spacing)
        points = c1.intersect(c2)
        if points[0].x < points[1].x:
            tag3, tag4 = CoilConnection.INSIDE_LEFT_VIA, CoilConnection.INSIDE_RIGHT_VIA
        else:
            tag3, tag4 = CoilConnection.INSIDE_RIGHT_VIA, CoilConnection.INSIDE_LEFT_VIA
        inside_via_3 = Via(points[0], config.via_diameter, config.via_hole_diameter, tag=tag3)
        inside_via_4 = Via(points[1], config.via_diameter, config.via_hole_diameter, tag=tag4)
        inside_vias_list = [inside_outer_via, inside_inner_via, inside_via_3, inside_via_4]
        inside_vias = {}
        for via in inside_vias_list:
            inside_vias[via.tag] = via

        # Outside vias
        left_line = Line.from_two_points(
            board_center,
            Point.polar(-config.coil_angle/2.0, config.board_radius)
        ).offset(config.via_diameter_w_spacing / 2.0)
        right_line = Line.from_two_points(
            board_center,
            Point.polar(config.coil_angle/2.0, config.board_radius)
        ).offset(-config.via_diameter_w_spacing / 2.0)
        outer_arc = Arc(
            Point.polar(-config.coil_angle/2.0, config.board_radius - config.board_outer_margin),
            Point.polar(config.coil_angle/2.0, config.board_radius - config.board_outer_margin),
            config.board_radius - config.board_outer_margin
        ).offset(-config.via_diameter / 2.0 + config.outer_vias_offset)
        inner_arc = Arc(
            Point.polar(-config.coil_angle/2.0, config.hole_radius + config.board_inner_margin),
            Point.polar(config.coil_angle/2.0, config.hole_radius + config.board_inner_margin),
            config.hole_radius + config.board_inner_margin
        ).offset(config.via_diameter / 2.0 - config.inner_vias_offset)
        points = [
            (left_line.intersect(outer_arc), CoilConnection.OUTSIDE_OUTER_LEFT_VIA),
            (right_line.intersect(outer_arc), CoilConnection.OUTSIDE_OUTER_RIGHT_VIA),
            (left_line.intersect(inner_arc), CoilConnection.OUTSIDE_INNER_LEFT_VIA),
            (right_line.intersect(inner_arc), CoilConnection.OUTSIDE_INNER_RIGHT_VIA),
        ]
        outside_vias = {}
        for point, tag in points:
            outside_vias[tag] = Via(point, config.via_diameter, config.via_hole_diameter, tag=tag)

        # Copy the vias on all coils
        vias = []
        for i in range(config.n_coils):
            for via in inside_vias.values():
                vias.append(via.rotated(board_center, 360.0 * i / config.n_coils))
            for via in outside_vias.values():
                vias.append(via.rotated(board_center, 360.0 * i / config.n_coils))

        # Specific generation settings for each layer : coil direction and vias connections
        layers_specs = {
            'top': {
                'anticlockwise': False,
                'outside_connection': CoilConnection.TERMINAL,
                'inside_connection': CoilConnection.INSIDE_OUTER_VIA,
            },
            'in1': {
                'anticlockwise': True,
                'outside_connection': CoilConnection.OUTSIDE_OUTER_RIGHT_VIA,
                'inside_connection': CoilConnection.INSIDE_OUTER_VIA,
            },
            'in2': {
                'anticlockwise': False,
                'outside_connection': CoilConnection.OUTSIDE_OUTER_RIGHT_VIA,
                'inside_connection': CoilConnection.INSIDE_LEFT_VIA,
            },
            'in3': {
                'anticlockwise': True,
                'outside_connection': CoilConnection.OUTSIDE_OUTER_LEFT_VIA,
                'inside_connection': CoilConnection.INSIDE_LEFT_VIA,
            },
            'in4': {
                'anticlockwise': False,
                'outside_connection': CoilConnection.OUTSIDE_OUTER_LEFT_VIA,
                'inside_connection': CoilConnection.INSIDE_INNER_VIA,
            },
            'in5': {
                'anticlockwise': True,
                'outside_connection': CoilConnection.OUTSIDE_INNER_LEFT_VIA,
                'inside_connection': CoilConnection.INSIDE_INNER_VIA,
            },
            'in6': {
                'anticlockwise': False,
                'outside_connection': CoilConnection.OUTSIDE_INNER_LEFT_VIA,
                'inside_connection': CoilConnection.INSIDE_RIGHT_VIA,
            },
            'bottom': {
                'anticlockwise': True,
                'outside_connection': CoilConnection.OUTSIDE_INNER_RIGHT_VIA,
                'inside_connection': CoilConnection.INSIDE_RIGHT_VIA,
            },
        }

        # Generate the coils on all layers
        layers = {}
        for layer_id, specs in layers_specs.items():
            if config.draw_only_layers is None or layer_id in config.draw_only_layers:
                # Generate the base Coil for this layer
                coil_B = Coil.generate(
                    angle = config.coil_angle,
                    outer_radius = config.board_radius - config.board_outer_margin - config.trace_width / 2.0,
                    inner_radius = config.hole_radius + config.board_inner_margin + config.trace_width / 2.0,
                    anticlockwise = specs['anticlockwise'],
                    trace_width = config.trace_width,
                    trace_spacing = config.trace_spacing,
                    outside_vias = outside_vias,
                    inside_vias = inside_vias,
                    outside_connection = specs['outside_connection'],
                    inside_connection = specs['inside_connection'],
                    max_turns = config.max_turns_per_layer,
                    construction_geometry = construction_geometry,
                )

                # Generate the other coils by rotating and mirroring the base coil
                # TODO : adapt for different number of coils
                coil_A = coil_B.rotated(board_center, -config.coil_angle)
                coil_C = coil_B.rotated(board_center, config.coil_angle)
                coil_A2 = coil_C.mirrored_x()
                coil_B2 = coil_B.mirrored_x()
                coil_C2 = coil_A.mirrored_x()

                # Add these coils to the current layer
                layers[layer_id] = [coil_A, coil_B, coil_C, coil_A2, coil_B2, coil_C2]

        # Add the other layers
        layers['outline'] = outline

        # Create the PCB
        return PCB(config, layers, vias, construction_geometry)
    
    def draw_svg(self, drawing: svg.Drawing) -> Self:
        """Draw this PCB on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        # Draw the layers
        for layer, objects in self.layers.items():
            if self.config.draw_only_layers is None or layer in self.config.draw_only_layers:
                for object in objects:
                    object.draw_svg(
                        drawing,
                        color = self.config.layers_color.get(layer),
                        thickness = self.config.trace_width,
                        dashes = "none"
                    )

        # Draw the vias
        if self.config.draw_vias:
            for via in self.vias:
                via.draw_svg(drawing, self.config.via_color, self.config.via_hole_color)

        # Draw the construction geometry
        if self.config.draw_construction_geometry:
            for object in self.construction_geometry:
                match object:
                    case svg.base.BaseElement():
                        drawing.add(object)
                    
                    case DrawableObject():
                        object.draw_svg(drawing)

        return self