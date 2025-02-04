from typing import Self
from enum import Enum
from config import Config, TerminalType
import svgwrite as svg
import math
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
            connection is None and side == CoilSide.OUTER or \
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
        if connection == CoilConnection.TERMINAL:
            return CoilConnection.TERMINAL
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
    
    def is_outer(connection: Self) -> bool:
        return connection in [CoilConnection.OUTSIDE_OUTER_LEFT_VIA, CoilConnection.OUTSIDE_OUTER_RIGHT_VIA]
    
    def is_inner(connection: Self) -> bool:
        return connection in [CoilConnection.OUTSIDE_INNER_LEFT_VIA, CoilConnection.OUTSIDE_INNER_RIGHT_VIA]
    
    def label(connection: Self) -> str:
        match connection:
            case CoilConnection.TERMINAL:
                return "Terminal"
            case CoilConnection.OUTSIDE_OUTER_LEFT_VIA:
                return "Outside_outer_left"
            case CoilConnection.OUTSIDE_OUTER_RIGHT_VIA:
                return "Outside_outer_right"
            case CoilConnection.OUTSIDE_INNER_LEFT_VIA:
                return "Outside_inner_left"
            case CoilConnection.OUTSIDE_INNER_RIGHT_VIA:
                return "Outside_inner_right"
            case CoilConnection.INSIDE_OUTER_VIA:
                return "Inside_outer"
            case CoilConnection.INSIDE_INNER_VIA:
                return "Inside_inner"
            case CoilConnection.INSIDE_LEFT_VIA:
                return "Inside_left"
            case CoilConnection.INSIDE_RIGHT_VIA:
                return "Inside_right"

class Via:
    """A via allowing connections between layers on the board"""

    def __init__(self, center: Point, diameter: float, hole_diameter: float, tag=None, label: str = None):
        self.center: Point = center
        self.diameter: float = diameter
        self.hole_diameter: float = hole_diameter
        self.tag = tag
        self.label = label
    
    def distance(self, object) -> float:
        """Calculate the distance between this Via and the given object

        Supported objects : Via, and any object supported by Point.distance().
        Throws a TypeError if any other object type is given.
        """
        match object:
            case Via():
                return self.center.distance(object.center) - self.diameter

            case _:
                return self.center.distance(object) - self.diameter / 2.0
    
    def rotated(self, rotation_center: Self, angle: float) -> Self:
        """Create a copy of this Via rotated around the given center point by the given angle"""
        return Via(
            center = self.center.rotated(rotation_center, angle),
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            tag = self.tag,
            label = self.label,
        )
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Via mirrored about the X axis"""
        return Via(
            center = self.center.mirrored_x(),
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            tag = self.tag,
            label = self.label,
        )
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Via mirrored about the Y axis"""
        return Via(
            center = self.center.mirrored_y(),
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            tag = self.tag,
            label = self.label,
        )
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement, via_color: str, hole_color: str, opacity: float = 1.0):
        """Draw this Via on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        # Group the two circles together
        group = drawing.g(label = self.label)

        # Draw the pad
        group.add(drawing.circle(
            self.center.to_viewport().as_tuple(),
            self.diameter / 2.0,
            stroke = "none",
            fill = via_color,
            opacity = opacity,
            label = f"{self.label}_pad" if self.label is not None else None,
        ))

        # Draw the hole
        group.add(drawing.circle(
            self.center.to_viewport().as_tuple(),
            self.hole_diameter / 2.0,
            stroke = "none",
            fill = hole_color,
            opacity = opacity,
            label = f"{self.label}_hole" if self.label is not None else None,
        ))

        # Add the group to the parent
        parent.add(group)

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

class Terminal:
    """A terminal for external connection on the board"""

    def __init__(self, center: Point, terminal_type: TerminalType, diameter: float, hole_diameter: float, tag=None, label: str = None):
        self.center: Point = center
        self.terminal_type: TerminalType = terminal_type
        self.diameter: float = diameter
        self.hole_diameter: float = hole_diameter
        self.tag = tag
        self.label = label
    
    def distance(self, object) -> float:
        """Calculate the distance between this Terminal and the given object

        Supported objects : Terminal, and any object supported by Point.distance().
        Throws a TypeError if any other object type is given.
        """
        match object:
            case Terminal():
                return self.center.distance(object.center) - self.diameter

            case _:
                return self.center.distance(object) - self.diameter / 2.0
    
    def rotated(self, rotation_center: Self, angle: float) -> Self:
        """Create a copy of this Terminal rotated around the given center point by the given angle"""
        return Terminal(
            center = self.center.rotated(rotation_center, angle),
            terminal_type = self.terminal_type,
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            tag = self.tag,
            label = self.label,
        )
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Terminal mirrored about the X axis"""
        return Terminal(
            center = self.center.mirrored_x(),
            terminal_type = self.terminal_type,
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            tag = self.tag,
            label = self.label,
        )
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Terminal mirrored about the Y axis"""
        return Terminal(
            center = self.center.mirrored_y(),
            terminal_type = self.terminal_type,
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            tag = self.tag,
            label = self.label,
        )
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement, pad_color: str, hole_color: str, opacity: float = 1.0, clip_path_id = None):
        """Draw this Terminal on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if self.terminal_type == TerminalType.NONE:
            return

        # Group the two circles together
        group = drawing.g(label = self.label)

        # Draw the pad
        if clip_path_id:
            group.add(drawing.circle(
                self.center.to_viewport().as_tuple(),
                self.diameter / 2.0,
                stroke = "none",
                fill = pad_color,
                opacity = opacity,
                clip_path = f"url(#{clip_path_id})" if clip_path_id is not None else None,
                label = f"{self.label}_pad" if self.label is not None else None,
            ))
        else:
            group.add(drawing.circle(
                self.center.to_viewport().as_tuple(),
                self.diameter / 2.0,
                stroke = "none",
                fill = pad_color,
                opacity = opacity,
                label = f"{self.label}_pad" if self.label is not None else None,
            ))

        # Draw the hole
        if self.terminal_type in [TerminalType.THROUGH_HOLE, TerminalType.CASTELLATED]:
            if clip_path_id:
                group.add(drawing.circle(
                    self.center.to_viewport().as_tuple(),
                    self.hole_diameter / 2.0,
                    stroke = "none",
                    fill = hole_color,
                    opacity = opacity,
                    clip_path = f"url(#{clip_path_id})" if clip_path_id is not None else None,
                    label = f"{self.label}_hole" if self.label is not None else None,
                ))
            else:
                group.add(drawing.circle(
                    self.center.to_viewport().as_tuple(),
                    self.hole_diameter / 2.0,
                    stroke = "none",
                    fill = hole_color,
                    opacity = opacity,
                    label = f"{self.label}_hole" if self.label is not None else None,
                ))

        # Add the group to the parent
        parent.add(group)

        return self

class Coil:
    """A single coil on the board"""

    def __init__(self, path: Path, rotation: float, n_turns: int, label: str = None):
        self.path: Path = path
        self.rotation: float = rotation
        self.n_turns: int = n_turns
        self.label = label

    def generate(
            config: Config,
            board_center: Point,
            angle: float,
            outer_radius: float,
            inner_radius: float,
            anticlockwise: bool,
            outside_vias: dict[CoilConnection, Via],
            inside_vias: dict[CoilConnection, Via],
            terminal: Terminal,
            outside_connection: CoilConnection,
            inside_connection: CoilConnection,
            max_outside_connection_length: float,
            mirror_outside_outer_via: bool,
            mirror_outside_inner_via: bool,
            construction_geometry: list,
        ) -> Self:
        """Generate a coil centered around the vertical axis based on the given parameters"""
        
        # If anticlockwise is set, the coil will be mirrored at the end, which means we need to connect it at
        # first to the opposite via relative to the Y axis so it will end up in the right place at the end.
        # If the mirror_outside_outer_via flags is set and the outside_connection of the coil is connected to
        # one of the outer vias, it should be connected to the opposite one instead, and similarly for the
        # inner connection.
        mirror_outside_via = \
            mirror_outside_outer_via and CoilConnection.is_outer(outside_connection) \
            or mirror_outside_inner_via and CoilConnection.is_inner(outside_connection)
        if (anticlockwise and not mirror_outside_via) or (not anticlockwise and mirror_outside_via):
            # Do not mirror if both flags are set as the two mirror operations cancel each other
            target_outside_connection = CoilConnection.mirrored_y(outside_connection)
        else:
            target_outside_connection = outside_connection

        # Distance between the centerline of each adjacent coil turn
        loop_offset = config.trace_spacing + config.trace_width
        
        # If the outer vias are used, calculate the best outer fillet radius of the outermost coil turn to fit the via
        # Twice the trace width is a sane value to start for any track width
        outer_fillet_radius = config.trace_width * 2
        if CoilConnection.OUTSIDE_OUTER_RIGHT_VIA in outside_vias or CoilConnection.OUTSIDE_OUTER_LEFT_VIA in outside_vias:
            outer_fillet_radius, success = Coil._compute_fillet(
                arc_radius = outer_radius,
                opposite_radius = inner_radius,
                angle = angle,
                initial_fillet_radius = outer_fillet_radius,
                via = outside_vias[CoilConnection.OUTSIDE_OUTER_RIGHT_VIA],
                trace_width = config.trace_width,
                trace_spacing = config.trace_spacing,
                construction_geometry = construction_geometry,
            )
            if not success:
                print("Warning : unable to compute a fillet for the outer side of the coil, check for collision with the via and try increasing outer_vias_offset")

        # Same calculation for the inner fillet radius
        inner_fillet_radius, success = Coil._compute_fillet(
            arc_radius = inner_radius,
            opposite_radius = outer_radius,
            angle = angle,
            initial_fillet_radius = config.trace_width * 2,
            via = outside_vias[CoilConnection.OUTSIDE_INNER_RIGHT_VIA],
            trace_width = config.trace_width,
            trace_spacing = config.trace_spacing,
            construction_geometry = construction_geometry,
        )
        if not success:
            print("Warning : unable to compute a fillet for the inner side of the coil, check for collision with the via and try increasing inner_vias_offset")

        # Calculate the spiral of the coil
        outer_radius_initial = outer_radius
        inner_radius_initial = inner_radius
        line_left = Line.from_two_points(board_center, Point.polar(-angle/2.0, outer_radius)).offset(loop_offset * 0.5 + outer_fillet_radius)
        arc_outer = Arc(Point.polar(-angle/2.0, outer_radius), Point.polar(angle/2.0, outer_radius), outer_radius)
        point_start = line_left.intersect(arc_outer)
        path = Path(point_start)
        n_turns = 0
        collision_via = None
        for i in range(config.max_turns_per_layer):
            # Construction geometry
            line_left = Line.from_two_points(board_center, Point.polar(-angle/2.0, outer_radius)).offset(loop_offset * (i + 0.5))
            line_right = Line.from_two_points(board_center, Point.polar(angle/2.0, outer_radius)).offset(-loop_offset * (i + 0.5))
            arc_outer = Arc(Point.polar(-angle/2.0, outer_radius), Point.polar(angle/2.0, outer_radius), outer_radius)
            inner_radius = inner_radius_initial + i * loop_offset
            arc_inner = Arc(Point.polar(-angle/2.0, inner_radius), Point.polar(angle/2.0, inner_radius), inner_radius)

            # Arc to outer right
            point_outer_right = line_right.intersect(arc_outer)
            arc = Arc(path.end_point, point_outer_right, outer_radius)
            closest_via, distance = Via.closest_in_list(inside_vias.values(), arc)
            if distance < config.trace_width / 2.0 + config.trace_spacing: # Collision
                collision_via = closest_via
            path.append_arc(point_outer_right, outer_radius, anticlockwise=False, fillet_radius=outer_fillet_radius, tag=CoilSide.OUTER)
            if collision_via is not None:
                break

            # Segment to inner right
            point_inner_right = line_right.intersect(arc_inner)
            if point_inner_right.x < config.trace_width:
                break
            segment = Segment(path.end_point, point_inner_right)
            closest_via, distance = Via.closest_in_list(inside_vias.values(), segment)
            if distance < config.trace_width / 2.0 + config.trace_spacing: # Collision
                collision_via = closest_via
            path.append_segment(point_inner_right, fillet_radius=outer_fillet_radius, tag=CoilSide.RIGHT)
            if collision_via is not None:
                break

            # Arc to inner left
            point_inner_left = line_left.intersect(arc_inner)
            arc = Arc(path.end_point, point_inner_left, outer_radius, reverse=True)
            closest_via, distance = Via.closest_in_list(inside_vias.values(), arc)
            if distance < config.trace_width / 2.0 + config.trace_spacing: # Collision
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
            if distance < config.trace_width / 2.0 + config.trace_spacing: # Collision
                collision_via = closest_via
            path.append_segment(point_outer_left, fillet_radius=inner_fillet_radius, tag=CoilSide.LEFT)
            if collision_via is not None:
                break

            # Reduce the fillet radius for the next loop
            min_radius = config.trace_width * 2
            outer_fillet_radius -= loop_offset
            if outer_fillet_radius < min_radius:
                outer_fillet_radius = min_radius
            inner_fillet_radius -= loop_offset
            if inner_fillet_radius < min_radius:
                inner_fillet_radius = min_radius
            
            # Count the number of turns in the coil
            n_turns += 1

        # Connect the start of the coil to the requested via or terminal
        while not CoilConnection.match_side(target_outside_connection, path.first().tag):
            # Pop the first elements from the path until we find the one that should be connected to the target via or terminal
            path.pop_first()
            path.pop_first() # Fillet
            if len(path.elements) == 1:
                raise ValueError("Error : unable to connect the outside of the coil to a via")
        if outside_connection == CoilConnection.TERMINAL or outside_connection is None:
            # Replace the first element in the path with a corner connected to the terminal
            if terminal is not None:
                corner = terminal.center.projected(path.first_geometry())
                connection_point = terminal.center
                connection_segment = Segment(corner, connection_point)
                if connection_segment.length() > max_outside_connection_length:
                    connection_point = corner + connection_segment.unit_vector() * max_outside_connection_length
            else:
                # If there is no terminal to connect to, simply add a stub to connect something later
                point = Point(0, outer_radius_initial)
                line = Line.from_two_points(board_center, point)
                corner = path.first_geometry().intersect(line)
                connection_point = point.offset(line, max_outside_connection_length)
            replaced_element = path.pop_first()
            fillet_radius = config.trace_width * 2
            if path.start_point.x < fillet_radius:
                fillet_radius = None # There isn't enough room for the fillet
            match replaced_element:
                case PathSegment():
                    path.prepend_segment(corner)
                    path.prepend_segment(connection_point, fillet_radius=fillet_radius)

                case PathArc():
                    path.prepend_arc(corner, replaced_element.radius, replaced_element.anticlockwise)
                    path.prepend_segment(connection_point, fillet_radius=fillet_radius)

        else:
            # Via to connect to
            try:
                target_via = outside_vias[outside_connection]
                if (anticlockwise and not mirror_outside_via) or (not anticlockwise and mirror_outside_via):
                    # Do not mirror if both flags are set as the two mirror operations cancel each other
                    target_via = target_via.mirrored_y()
            except KeyError:
                raise ValueError(f"Invalid outside_connection : {outside_connection}")
            
            # Connect the outside of the coil to the target via
            if target_outside_connection in [CoilConnection.OUTSIDE_OUTER_RIGHT_VIA, CoilConnection.OUTSIDE_INNER_LEFT_VIA]:
                segment = Segment(path.first().p2, path.start_point)
                tangent_arc = segment.tangent_arc_through_point(target_via.center)
                path.prepend_arc(target_via.center, tangent_arc.radius, anticlockwise = segment.p2 == tangent_arc.p1)
            elif target_outside_connection in [CoilConnection.OUTSIDE_OUTER_LEFT_VIA, CoilConnection.OUTSIDE_INNER_RIGHT_VIA]:
                first = path.first()
                target_via_radius = board_center.distance(target_via.center)
                if target_via_radius <= first.radius:
                    # The via is inside the arc radius, connect it with a tangeant arc
                    arc = Arc(first.p2, path.start_point, first.radius, reverse = not first.anticlockwise)
                    tangent_arc = arc.tangent_arc_through_point(target_via.center, at_start = not first.anticlockwise)
                    path.prepend_arc(target_via.center, tangent_arc.radius, anticlockwise = first.anticlockwise)
                else:
                    # The via is outside the arc radius, connect it with a fillet and a segment
                    path.prepend_segment(target_via.center, fillet_radius = (target_via_radius - first.radius) / 2.0)

        # Connect the inside of the coil to the requested via
        if inside_connection is not None:
            # Via to connect to
            try:
                target_via = inside_vias[inside_connection]
                if (anticlockwise and not mirror_outside_via) or (not anticlockwise and mirror_outside_via):
                    # Do not mirror if both flags are set as the two mirror operations cancel each other
                    target_via = target_via.mirrored_y()
            except KeyError:
                raise ValueError(f"Invalid inside_connection : {inside_connection}")

            # Pop the last elements from the path until we find the one that should be connected to the target via
            while not inside_connection.match_side(path.last().tag):
                path.pop()
                if len(path.elements) == 1:
                    raise ValueError("Error : unable to connect the inside of the coil to a via, try increasing hole_diameter or reducing n_slots_per_phase")
            
            # Replace the last element in the path with a corner connected to the target via
            last_element = path.last_geometry()
            replaced_element = path.pop()
            fillet_radius = config.trace_width * 2
            corner = target_via.center.projected(last_element)
            if not last_element.contains_point(corner):
                corner = last_element.midpoint()
            if replaced_element.p2.distance(corner) < fillet_radius:
                fillet_radius = None # There is no room for the fillet
            match replaced_element:
                case PathSegment():
                    path.append_segment(corner)
                    path.append_arc(target_via.center, board_center.distance(target_via.center), anticlockwise=corner.x >= target_via.center.x, fillet_radius=fillet_radius)

                case PathArc():
                    path.append_arc(corner, replaced_element.radius, replaced_element.anticlockwise)
                    path.append_segment(target_via.center, fillet_radius=fillet_radius)

        # Return the coil based on this path, mirrored relative to the Y axis if anticlockwise
        if anticlockwise:
            path = path.mirrored_y()
        return Coil(path, 0.0, n_turns)

    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Coil rotated around the given center point by the given angle"""
        path = self.path.rotated(center, angle)
        return Coil(
            path = path,
            rotation = self.rotation + angle,
            n_turns = self.n_turns,
            label = self.label,
        )
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Coil mirrored about the X axis"""
        path = self.path.mirrored_x()
        return Coil(
            path = path,
            rotation = self.rotation,
            n_turns = self.n_turns,
            label = self.label,
        )
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Coil mirrored about the Y axis"""
        path = self.path.mirrored_y()
        return Coil(
            path = path,
            rotation = -self.rotation,
            n_turns = self.n_turns,
            label = self.label,
        )
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement, color=None, opacity=None, thickness=None, dashes=None):
        """Draw this Coil on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        self.path.draw_svg(
            drawing = drawing,
            parent = parent,
            color = color,
            opacity = opacity,
            thickness = thickness,
            dashes = dashes,
            label = self.label,
        )
        return self
    
    def _compute_fillet(
            arc_radius: float,
            opposite_radius: float,
            angle: float,
            initial_fillet_radius: float,
            via: Via,
            trace_width: float,
            trace_spacing: float,
            construction_geometry = None
        ) -> tuple[float, bool]:
        radius = initial_fillet_radius
        arc = Arc(Point.polar(-angle/2.0, arc_radius), Point.polar(angle/2.0, arc_radius), arc_radius)
        segment_base = Segment(Point.polar(angle/2.0, opposite_radius), Point.polar(angle/2.0, arc_radius))
        segment = segment_base.offset_closest_to(arc.p1, (trace_spacing + trace_width) * 0.5)
        corner = segment.intersect(arc)
        start_point = arc.midpoint()
        end_point = segment.midpoint()
        while radius < arc_radius:
            # Try a fillet with its radius increased by 10% or the width of a trace, whichever is larger
            try_fillet_radius = max(radius * 0.1, radius + trace_width)

            # Construct a path with half an outer arc and half a right segment, and a fillet of the current radius in the middle
            path = Path(start_point)
            path.append_arc(corner, arc_radius, anticlockwise=False)
            path.append_segment(end_point, fillet_radius=try_fillet_radius, suppress_warning=True)

            # Get the actual fillet radius and construct an arc at the same place as the fillet
            fillet_real_radius = path.elements[-2].radius
            fillet_arc_p1 = path.elements[-3].p2
            fillet_arc_p2 = path.elements[-2].p2
            anticlockwise_fillet = Vector.from_two_points(corner, fillet_arc_p1).cross(Vector.from_two_points(corner, fillet_arc_p2)) < 0
            fillet_arc = Arc(fillet_arc_p1, fillet_arc_p2, radius=fillet_real_radius, reverse=anticlockwise_fillet)

            # Check if the fillet has the same radius as the one we are trying to make
            if not math.isclose(try_fillet_radius, fillet_real_radius, abs_tol=1e-9):
                return radius, False

            # Check if the end of the fillet is too close to the center line
            if fillet_arc.p2.x < trace_width * 2.0:
                return radius, False

            # Check if the fillet is large enough to prevent a collision with the via
            if via.distance(fillet_arc) >= trace_width / 2.0 + trace_spacing:
                return try_fillet_radius, True

            # This radius looks ok, save it and try again with a larger one
            radius = try_fillet_radius

class Link:
    """A trace on the board linking two coils"""

    def __init__(self, path: Path, trace_width: float, layer: str, label: str = None):
        self.path: Path = path
        self.trace_width: float = trace_width
        self.layer: str = layer
        self.label: str = label

    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Link rotated around the given center point by the given angle"""
        path = self.path.rotated(center, angle)
        return Link(
            path = path,
            trace_width = self.trace_width,
            layer = self.layer,
            label = self.label
        )
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Link mirrored about the X axis"""
        path = self.path.mirrored_x()
        return Link(
            path = path,
            trace_width = self.trace_width,
            layer = self.layer,
            label = self.label
        )
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Link mirrored about the Y axis"""
        path = self.path.mirrored_y()
        return Link(
            path = path,
            trace_width = self.trace_width,
            layer = self.layer,
            label = self.label
        )
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement, color=None, opacity=None, thickness=None, dashes=None):
        """Draw this Link on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        self.path.draw_svg(
            drawing = drawing,
            parent = parent,
            color = color,
            opacity = opacity,
            thickness = thickness,
            dashes = dashes,
            label = self.label,
        )
        return self

class Outline:
    """The shape and dimensions of the board"""

    def __init__(self, board_center: Point, board_diameter: float, hole_diameter: float):
        self.board_center = board_center
        self.board_diameter = board_diameter
        self.hole_diameter = hole_diameter
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement, color=None, thickness=None, dashes=None):
        """Draw this Link on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        parent.add(drawing.circle(
            self.board_center.to_viewport().as_tuple(),
            self.board_diameter / 2.0,
            stroke = color,
            stroke_width = thickness,
            stroke_dasharray = dashes,
            fill = "none",
            label = f"Outline_edge",
        ))
        parent.add(drawing.circle(
            self.board_center.to_viewport().as_tuple(),
            self.hole_diameter / 2.0,
            stroke = color,
            stroke_width = thickness,
            stroke_dasharray = dashes,
            fill = "none",
            label = f"Outline_hole",
        ))
        return self

class PCB:
    """A PCB containing layers"""

    def __init__(self,
        config: Config,
        board_center: Point,
        outline: list,
        coils: dict,
        vias: list[Via],
        terminals: list[Terminal],
        links: list[Link],
        construction_geometry: list
    ):
        self.config = config
        self.board_center = board_center
        self.outline = outline
        self.coils = coils
        self.vias = vias
        self.terminals = terminals
        self.links = links
        self.construction_geometry = construction_geometry
    
    def generate(config: Config):
        """Generate a new PCB based on the given config"""
        
        coils = {}
        vias = []
        terminals = []
        links = []
        construction_geometry = []

        # Names of the coils : A1, A2, A3, B1, ...
        coil_names = []
        for slot in range(config.n_slots_per_phase):
            for phase in range(config.n_phases):
                coil_names.append(f"{chr(ord('A') + phase)}{slot + 1}")

        # Center the board on the origin
        board_center = Point.origin()

        # Board outline
        outline = Outline(board_center, config.board_diameter, config.hole_diameter)

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

        # Specific generation settings for each layer : coil direction and vias connections
        layers_specs = {}
        match config.n_layers:
            case 2:
                outside_connections = [
                    CoilConnection.TERMINAL if config.terminal_type != TerminalType.NONE else None,
                    CoilConnection.OUTSIDE_INNER_RIGHT_VIA,
                ]
                inside_connections = [
                    CoilConnection.INSIDE_OUTER_VIA,
                ]
                optimise_outer_vias = False
                optimise_inner_vias = True
            case 4:
                # For 4-layer boards, there is two possibilities : other place the two inner vias to have
                # only traces on the outer side of the coils and maximize useful traces there, or place
                # the vias along the middle line (using the "optimise" flags) to minimize the place taken
                # by the vias. This choice is left to the user in the config.
                if config.four_layers_inside_vias:
                    outside_connections = [
                        CoilConnection.TERMINAL if config.terminal_type != TerminalType.NONE else None,
                        CoilConnection.OUTSIDE_INNER_LEFT_VIA,
                        CoilConnection.OUTSIDE_INNER_RIGHT_VIA,
                    ]
                    inside_connections = [
                        CoilConnection.INSIDE_OUTER_VIA,
                        CoilConnection.INSIDE_INNER_VIA,
                    ]
                    optimise_outer_vias = False
                    optimise_inner_vias = False
                else:
                    outside_connections = [
                        CoilConnection.TERMINAL if config.terminal_type != TerminalType.NONE else None,
                        CoilConnection.OUTSIDE_OUTER_RIGHT_VIA,
                        CoilConnection.OUTSIDE_INNER_RIGHT_VIA,
                    ]
                    inside_connections = [
                        CoilConnection.INSIDE_OUTER_VIA,
                        CoilConnection.INSIDE_INNER_VIA,
                    ]
                    optimise_outer_vias = True
                    optimise_inner_vias = True
            case 6:
                outside_connections = [
                    CoilConnection.TERMINAL if config.terminal_type != TerminalType.NONE else None,
                    CoilConnection.OUTSIDE_OUTER_RIGHT_VIA,
                    CoilConnection.OUTSIDE_INNER_LEFT_VIA,
                    CoilConnection.OUTSIDE_INNER_RIGHT_VIA,
                ]
                inside_connections = [
                    CoilConnection.INSIDE_OUTER_VIA,
                    CoilConnection.INSIDE_RIGHT_VIA,
                    CoilConnection.INSIDE_INNER_VIA,
                ]
                optimise_outer_vias = True
                optimise_inner_vias = False
            case 8:
                outside_connections = [
                    CoilConnection.TERMINAL if config.terminal_type != TerminalType.NONE else None,
                    CoilConnection.OUTSIDE_OUTER_RIGHT_VIA,
                    CoilConnection.OUTSIDE_OUTER_LEFT_VIA,
                    CoilConnection.OUTSIDE_INNER_LEFT_VIA,
                    CoilConnection.OUTSIDE_INNER_RIGHT_VIA,
                ]
                inside_connections = [
                    CoilConnection.INSIDE_OUTER_VIA,
                    CoilConnection.INSIDE_LEFT_VIA,
                    CoilConnection.INSIDE_INNER_VIA,
                    CoilConnection.INSIDE_RIGHT_VIA,
                ]
                optimise_outer_vias = False
                optimise_inner_vias = False
        for i, layer_id in enumerate(config.layers):
            layers_specs[layer_id] = {
                'anticlockwise': i % 2 == 1,
                'outside_connection': outside_connections[int((i + 1) / 2)],
                'inside_connection': inside_connections[int(i / 2)],
            }

        # Inside vias
        inside_vias = {}
        coil_center_radius = ((config.board_radius - config.board_outer_margin) + (config.hole_radius + config.board_inner_margin)) / 2.0
        coil_center = Point(0, coil_center_radius)
        tangential_length = 2 * math.pi * coil_center_radius * config.coil_angle / 360.
        radial_length = (config.board_radius - config.board_outer_margin) - (config.hole_radius + config.board_inner_margin)
        match config.n_layers:
            case 2:
                # A single via placed in the middle
                via = Via(coil_center, config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                inside_vias_list = [via]
            case 4:
                # Two vias placed vertically
                inside_outer_via = Via(coil_center + Vector(0, config.via_diameter_w_spacing / 2.0), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                inside_inner_via = Via(coil_center - Vector(0, config.via_diameter_w_spacing / 2.0), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
                inside_vias_list = [inside_outer_via, inside_inner_via]
                pass
            case 6:
                # Three vias, placed either vertically, or on a triangular shape
                if tangential_length >= radial_length:
                    # The coil is wider than it is tall, place the via in a triangle
                    inside_inner_via = Via(coil_center + Vector(0, -config.via_diameter_w_spacing / 2.0), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
                    inside_outer_left_via = Via(inside_inner_via.center + Vector.polar(-30, config.via_diameter_w_spacing), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                    inside_outer_right_via = Via(inside_inner_via.center + Vector.polar(30, config.via_diameter_w_spacing), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_RIGHT_VIA)
                    inside_vias_list = [inside_outer_left_via, inside_outer_right_via, inside_inner_via]
                else:
                    # The coil is taller than it is wide, place the via vertically
                    inside_outer_via = Via(coil_center + Vector(0, config.via_diameter_w_spacing), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                    inside_middle_via = Via(coil_center, config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_RIGHT_VIA)
                    inside_inner_via = Via(coil_center + Vector(0, -config.via_diameter_w_spacing), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
                    inside_vias_list = [inside_outer_via, inside_middle_via, inside_inner_via]
            case 8:
                # Four vias placed in an optimised diamond shape. If the coil is wider than tall (tangential length > radial length),
                # keep the base diamond shape on its side with the two vertical vias adjacent to each other. Otherwise, optimise
                # the shape to fit as many radial lines as possible.
                step = config.trace_width / 10.0
                sep = 0.0
                if radial_length >= tangential_length:
                    while True:
                        # The loop starts with the two vertical vias closest to each other, and the diamond shape is progressively
                        # stretched vertically until the line connecting the inner via and the left via is parallel to the left left.
                        # This ensures that this leaves as much room as possible for the radial traces.

                        # Compute the shape of the diamond
                        sep_test = sep + step
                        inside_outer_via = Via(coil_center + Vector(0, (config.via_diameter_w_spacing + sep_test) / 2.0), config.via_diameter, config.via_hole_diameter)
                        inside_inner_via = Via(coil_center + Vector(0, -(config.via_diameter_w_spacing + sep_test) / 2.0), config.via_diameter, config.via_hole_diameter)
                        c1 = Circle(inside_outer_via.center, config.via_diameter_w_spacing)
                        c2 = Circle(inside_inner_via.center, config.via_diameter_w_spacing)
                        points = c1.intersect(c2)
                        if points is None:
                            # No intersection, the vertical vias are too far appart : this shouldn't happen
                            break
                        if points[0].x < points[1].x:
                            point_left = points[0]
                        else:
                            point_left = points[1]
                        
                        # Compute the cross product between the left line, and the line connecting the inner via and the left via
                        v1 = Line.from_two_points(board_center, Point.polar(-config.coil_angle/2.0, config.board_radius)).unit_vector()
                        v2 = Vector.from_two_points(inside_inner_via.center, point_left)
                        if v1.cross(v2) < 0:
                            # The cross product switched sign : we just passed the parallel
                            break
                        
                        # Make sure there is enough spacing between the horizontal vias
                        inside_via_3 = Via(points[0], config.via_diameter, config.via_hole_diameter)
                        inside_via_4 = Via(points[1], config.via_diameter, config.via_hole_diameter)
                        if inside_via_3.distance(inside_via_4) < config.trace_spacing:
                            break
                        sep = sep_test
                inside_outer_via = Via(coil_center + Vector(0, (config.via_diameter_w_spacing + sep) / 2.0), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                inside_inner_via = Via(coil_center - Vector(0, (config.via_diameter_w_spacing + sep) / 2.0), config.via_diameter, config.via_hole_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
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
        for via in inside_vias_list:
            inside_vias[via.tag] = via

        # Outside vias
        outside_vias = {}
        left_line = Line.from_two_points(
            board_center,
            Point.polar(-config.coil_angle/2.0, config.board_radius)
        )
        left_line_offset = left_line.offset(config.via_diameter_w_spacing / 2.0)
        right_line = Line.from_two_points(
            board_center,
            Point.polar(config.coil_angle/2.0, config.board_radius)
        )
        right_line_offset = right_line.offset(-config.via_diameter_w_spacing / 2.0)
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
        if optimise_outer_vias:
            # If this flag is set, there is only one outer via used, so place it exactly
            # on the line between two coils to optimize space
            points_outer = [
                (left_line.intersect(outer_arc), CoilConnection.OUTSIDE_OUTER_LEFT_VIA),
                (right_line.intersect(outer_arc), CoilConnection.OUTSIDE_OUTER_RIGHT_VIA),
            ]
        else:
            # Offset the point to fit vias on both sides of the line
            points_outer = [
                (left_line_offset.intersect(outer_arc), CoilConnection.OUTSIDE_OUTER_LEFT_VIA),
                (right_line_offset.intersect(outer_arc), CoilConnection.OUTSIDE_OUTER_RIGHT_VIA),
            ]
        if optimise_inner_vias:
            # If this flag is set, there is only one inner via used, so place it exactly
            # on the line between two coils to optimize space
            points_inner = [
                (left_line.intersect(inner_arc), CoilConnection.OUTSIDE_INNER_LEFT_VIA),
                (right_line.intersect(inner_arc), CoilConnection.OUTSIDE_INNER_RIGHT_VIA),
            ]
        else:
            # Offset the point to fit vias on both sides of the line
            points_inner = [
                (left_line_offset.intersect(inner_arc), CoilConnection.OUTSIDE_INNER_LEFT_VIA),
                (right_line_offset.intersect(inner_arc), CoilConnection.OUTSIDE_INNER_RIGHT_VIA),
            ]
        points = points_outer + points_inner
        for point, connection in points:
            if connection in outside_connections:
                outside_vias[connection] = Via(point, config.via_diameter, config.via_hole_diameter, tag=connection)
        if CoilConnection.OUTSIDE_INNER_LEFT_VIA in outside_vias \
                and outside_vias[CoilConnection.OUTSIDE_INNER_LEFT_VIA].distance(outside_vias[CoilConnection.OUTSIDE_INNER_RIGHT_VIA]) < 0:
            print("Warning : collision between the vias closer to the center of the board")

        # Copy the vias for all coils
        for slot in range(config.n_slots_per_phase):
            for phase in range(config.n_phases):
                coil_idx = slot * config.n_phases + phase
                for connection, via in inside_vias.items():
                    if slot % 2 == 1:
                        via = via.mirrored_y()
                    via = via.rotated(board_center, 360.0 * coil_idx / config.n_coils)
                    via.label = f"Via_{coil_names[coil_idx]}_{CoilConnection.label(connection)}"
                    vias.append(via)
                for connection, via in outside_vias.items():
                    # If the optimise flag is set for this via, it should be mirrored twice, so don't mirror it
                    # at all as both operations would cancel each other
                    if slot % 2 == 1 \
                            and not (optimise_outer_vias and slot % 2 == 1 and CoilConnection.is_outer(connection)) \
                            and not (optimise_inner_vias and slot % 2 == 1 and CoilConnection.is_inner(connection)):
                        via = via.mirrored_y()
                    via = via.rotated(board_center, 360.0 * coil_idx / config.n_coils)
                    via.label = f"Via_{coil_names[coil_idx]}_{CoilConnection.label(connection)}"
                    vias.append(via)
        
        # Terminal
        terminal_base = None
        match config.terminal_type:
            case TerminalType.THROUGH_HOLE | TerminalType.SMD:
                # Align the terminal outside the coil
                y = config.board_radius - config.board_outer_margin + config.trace_spacing + config.terminal_offset + config.terminal_diameter / 2.0
            case TerminalType.CASTELLATED:
                # Align the terminal on the edge of the board
                y = config.board_radius
        terminal_base = Terminal(Point(0.0, y), config.terminal_type, config.terminal_diameter, config.terminal_hole_diameter)

        # Copy the terminals for all phases and the COM point
        if terminal_base:
            if config.link_series_coils:
                n_terminals = config.n_phases
            else:
                n_terminals = config.n_coils
            for i in range(n_terminals):
                terminal = terminal_base.rotated(board_center, 360.0 * i / config.n_coils)
                terminal.label = f"Terminal_{coil_names[i]}"
                terminals.append(terminal)
            if config.link_series_coils:
                if config.link_com:
                    terminal = terminal_base.rotated(board_center, - 360.0 / config.n_coils)
                    terminal.label = f"Terminal_COM"
                    terminals.append(terminal)
                else:
                    for i in range(config.n_phases):
                        terminal = terminal_base.rotated(board_center, 360.0 * -(i + 1) / config.n_coils)
                        terminal.label = f"Terminal_{coil_names[-(i+1)]}"
                        terminals.append(terminal)

        # Generate the coils on all layers
        coil_outside_connection_length = config.trace_width * 2
        for layer_id in config.layers:
            specs = layers_specs[layer_id]

            # Generate the base Coil for this layer
            coil_even = Coil.generate(
                config = config,
                board_center = board_center,
                angle = config.coil_angle,
                outer_radius = config.board_radius - config.board_outer_margin - config.trace_width / 2.0,
                inner_radius = config.hole_radius + config.board_inner_margin + config.trace_width / 2.0,
                anticlockwise = specs['anticlockwise'],
                outside_vias = outside_vias,
                inside_vias = inside_vias,
                terminal = terminal_base,
                outside_connection = specs['outside_connection'],
                inside_connection = specs['inside_connection'],
                max_outside_connection_length = coil_outside_connection_length,
                mirror_outside_outer_via = False,
                mirror_outside_inner_via = False,
                construction_geometry = construction_geometry,
            )
            if optimise_outer_vias or optimise_inner_vias:
                # If one of these flags is set, the odd coils are not exactly the same as the even coils mirrored
                # around the Y axis because of the position of the outside via, so we need to generate a new coil
                # with the `mirror_outside_outer_via` and `mirror_outside_inner_via` parameters set accordingly
                # before mirroring the coil.
                coil_odd = Coil.generate(
                    config = config,
                    board_center = board_center,
                    angle = config.coil_angle,
                    outer_radius = config.board_radius - config.board_outer_margin - config.trace_width / 2.0,
                    inner_radius = config.hole_radius + config.board_inner_margin + config.trace_width / 2.0,
                    anticlockwise = specs['anticlockwise'],
                    outside_vias = outside_vias,
                    inside_vias = inside_vias,
                    terminal = terminal_base,
                    outside_connection = specs['outside_connection'],
                    inside_connection = specs['inside_connection'],
                    max_outside_connection_length = coil_outside_connection_length,
                    mirror_outside_outer_via = optimise_outer_vias,
                    mirror_outside_inner_via = optimise_inner_vias,
                    construction_geometry = construction_geometry,
                ).mirrored_y()
            else:
                # In other cases, the odd coils are simply the mirror of the even coils, there is no need
                # to generate a full new path as this is a relatively expensive operation.
                coil_odd = coil_even.mirrored_y()

            # Generate the other coils by rotating the base coils
            layer_coils = []
            for slot in range(config.n_slots_per_phase):
                for phase in range(config.n_phases):
                    coil_index = slot * config.n_phases + phase
                    angle = 360.0 * coil_index / config.n_coils
                    if slot % 2 == 0:
                        coil = coil_even
                    else:
                        coil = coil_odd
                    coil = coil.rotated(board_center, angle)
                    coil.label = f"Coil_{coil_names[coil_index]}_{layer_id}"
                    layer_coils.append(coil)

            # Add these coils to the current layer
            coils[layer_id] = layer_coils

        # Connect the coils that are in series, such as A_1 with A_2
        link_vias = []
        if config.link_series_coils and config.n_slots_per_phase >= 2:
            if config.n_layers >= config.n_phases:
                # Base path for inner connections (even coils)
                # Draw a path offset toward the inner side of the board that connects these two points for the first coil
                connection_point_1 = coils[config.layers[-1]][0].path.start_point
                connection_point_2 = connection_point_1
                if not optimise_inner_vias:
                    connection_point_2 = connection_point_2.mirrored_y()
                connection_point_2 = connection_point_2.rotated(board_center, config.coil_angle * config.n_phases)
                vias_circle_radius = board_center.distance(connection_point_1)
                trace_center_radius = vias_circle_radius - config.via_diameter / 2.0 - config.trace_spacing - config.series_link_inner_trace_width / 2.0 - config.series_link_inner_offset
                base_path_inner = Path(connection_point_1)
                if math.isclose(trace_center_radius, vias_circle_radius):
                    base_path_inner.append_arc(connection_point_2, trace_center_radius, anticlockwise=False)
                else:
                    circle = Circle(board_center, trace_center_radius)
                    line1 = Line.from_two_points(board_center, connection_point_1)
                    line2 = Line.from_two_points(board_center, connection_point_2)
                    corner1 = connection_point_1.closest(line1.intersect(circle))
                    corner2 = connection_point_2.closest(line2.intersect(circle))
                    fillet_radius = connection_point_1.distance(corner1) / 2.0
                    base_path_inner.append_segment(corner1)
                    base_path_inner.append_arc(corner2, trace_center_radius, anticlockwise=False, fillet_radius=fillet_radius)
                    base_path_inner.append_segment(connection_point_2, fillet_radius=fillet_radius)

                # Base path for outer connections (odd coils)
                # Draw a path offset toward the outer side of the board that connects these two points for the first coil
                if config.n_slots_per_phase >= 4:
                    connection_point_1 = coils[config.layers[0]][0].path.start_point
                    connection_point_2 = connection_point_1.mirrored_y().rotated(board_center, config.coil_angle * config.n_phases)
                    outer_vias_radius = config.board_radius - config.board_outer_margin + config.outer_vias_offset + config.trace_spacing + config.series_link_outer_trace_width / 2.0
                    vias_circle_radius = config.board_radius - config.board_outer_margin + max(config.via_diameter / 2.0, config.series_link_outer_trace_width / 2.0) + config.trace_spacing
                    trace_center_radius = max(vias_circle_radius + config.via_diameter / 2.0 + config.trace_spacing + config.series_link_outer_trace_width / 2.0, outer_vias_radius) + config.series_link_outer_offset
                    circle_vias = Circle(board_center, vias_circle_radius)
                    circle_trace = Circle(board_center, trace_center_radius)
                    line1 = Line.from_two_points(board_center, connection_point_1)
                    line2 = Line.from_two_points(board_center, connection_point_2)
                    via_pos_1 = connection_point_1.closest(line1.intersect(circle_vias))
                    via_pos_2 = connection_point_2.closest(line2.intersect(circle_vias))
                    corner1 = connection_point_1.closest(line1.intersect(circle_trace))
                    corner2 = connection_point_2.closest(line2.intersect(circle_trace))
                    fillet_radius = via_pos_1.distance(corner1) / 2.0
                    base_path_outer = Path(via_pos_1)
                    base_path_outer.append_segment(corner1)
                    base_path_outer.append_arc(corner2, trace_center_radius, anticlockwise=False, fillet_radius=fillet_radius)
                    base_path_outer.append_segment(via_pos_2, fillet_radius=fillet_radius)

                # Create the links on different layers by rotating the base path accordingly for each phase
                for slot_pair in range(config.n_slots_per_phase - 1):
                    for phase in range(config.n_phases):
                        angle = ((slot_pair * config.n_phases) + phase) * config.coil_angle
                        label = f"Link_{coil_names[slot_pair * config.n_phases + phase]}_{coil_names[(slot_pair + 1) * config.n_phases + phase]}"
                        if slot_pair % 2 == 0:
                            path = base_path_inner.rotated(board_center, angle)
                            link = Link(path, config.series_link_inner_trace_width, config.layers[phase], label=label)
                            links.append(link)
                        else:
                            path = base_path_outer.rotated(board_center, angle)
                            link = Link(path, config.series_link_outer_trace_width, config.layers[phase], label=label)
                            links.append(link)
                            via1_label = f"Link_via_{coil_names[slot_pair * config.n_phases + phase]}"
                            via1 = Via(via_pos_1, config.via_diameter, config.via_hole_diameter, label=via1_label).rotated(board_center, angle)
                            via2_label = f"Link_via_{coil_names[(slot_pair + 1) * config.n_phases + phase]}"
                            via2 = Via(via_pos_2, config.via_diameter, config.via_hole_diameter, label=via2_label).rotated(board_center, angle)
                            link_vias.extend([via1, via2])
            else:
                print("Warning : unable to link the coils, not enough layers")
        vias.extend(link_vias)

        # Connect the common point in a wye (star) configuration, on the top layer
        link_com_connection_point = None
        if config.link_com:
            layer_id = config.layers[0]
            connection_point_1 = coils[layer_id][0].path.start_point
            connection_point_2 = connection_point_1.rotated(board_center, -config.coil_angle)
            terminal_circle_radius = board_center.distance(terminal_base.center) if terminal_base is not None else 0
            intermediate_circle_radius = config.board_radius - config.board_outer_margin + max(config.com_link_trace_width / 2.0 + config.trace_spacing, coil_outside_connection_length)
            outer_vias_radius = config.board_radius - config.board_outer_margin + config.outer_vias_offset + config.trace_spacing + config.com_link_trace_width / 2.0
            if config.terminal_type == TerminalType.CASTELLATED:
                trace_center_radius = max(intermediate_circle_radius, outer_vias_radius) + config.com_link_offset
            else:
                trace_center_radius = max(terminal_circle_radius, intermediate_circle_radius, outer_vias_radius) + config.com_link_offset
            circle_intermediate = Circle(board_center, intermediate_circle_radius)
            circle_trace = Circle(board_center, trace_center_radius)
            line1 = Line.from_two_points(board_center, connection_point_1)
            line2 = Line.from_two_points(board_center, connection_point_2)
            intermediate_point_1 = connection_point_1.closest(line1.intersect(circle_intermediate))
            intermediate_point_2 = connection_point_2.closest(line2.intersect(circle_intermediate))
            corner1 = connection_point_1.closest(line1.intersect(circle_trace))
            corner2 = connection_point_2.closest(line2.intersect(circle_trace))
            fillet_radius = intermediate_point_1.distance(corner1) / 2.0
            base_path_com = Path(intermediate_point_1)
            if math.isclose(trace_center_radius, intermediate_circle_radius):
                base_path_com.append_arc(corner2, trace_center_radius, anticlockwise=True)
            else:
                base_path_com.append_segment(corner1)
                base_path_com.append_arc(corner2, trace_center_radius, anticlockwise=True, fillet_radius=fillet_radius)
                base_path_com.append_segment(intermediate_point_2, fillet_radius=fillet_radius)
            for i in range(config.n_phases - 1):
                path = base_path_com.rotated(board_center, -(i + 1) * config.coil_angle)
                label = f"COM_{coil_names[-(i+2)]}_{coil_names[-(i+1)]}"
                link = Link(path, config.com_link_trace_width, layer_id, label=label)
                links.append(link)
            link_com_connection_point = intermediate_point_1

        # Connect the terminals and linking vias on the first layer
        layer_id = config.layers[0]
        if config.draw_only_layers is None or layer_id in config.draw_only_layers:
            for i in range(config.n_coils):
                if terminal_base and (i < config.n_phases or (not config.link_com and i >= config.n_coils - config.n_phases) or (i == config.n_coils - 1) or not config.link_series_coils):
                    coils[layer_id][i].path.prepend_segment(terminal_base.rotated(board_center, 360.0 * i / config.n_coils).center)
                elif config.link_series_coils and link_vias and config.n_slots_per_phase >= 4 and i >= config.n_phases and i < (config.n_slots_per_phase - 1) * config.n_phases:
                    radius = board_center.distance(link_vias[0].center)
                    coils[layer_id][i].path.prepend_segment(Point.polar(0, radius).rotated(board_center, 360.0 * i / config.n_coils))
                elif link_com_connection_point and config.link_com and i >= config.n_coils - config.n_phases:
                    coils[layer_id][i].path.prepend_segment(link_com_connection_point.rotated(board_center, 360.0 * i / config.n_coils))

        # Create the PCB
        return PCB(
            config = config,
            board_center = board_center,
            outline = outline,
            coils = coils,
            vias = vias,
            terminals = terminals,
            links = links,
            construction_geometry = construction_geometry,
        )
    
    def draw_svg(self, drawing: svg.Drawing) -> Self:
        """Draw this PCB on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        # Create the SVG layers to organize the drawing.
        # Layer groups created later will be displayed on top of the previous ones, so the board layers are reversed
        # to have the top layer on top and bottom layer on the bottom, as expected.
        svg_layers_ids = list(reversed(self.config.layers)) + ['outline', 'vias', 'terminals', 'magnets', 'construction']
        svg_layers = {}
        for layer_id in svg_layers_ids:
            svg_layers[layer_id] = drawing.layer(label=layer_id.capitalize())

        # Draw the coils
        for layer_id in self.config.layers:
            if self.config.draw_only_layers is None or layer_id in self.config.draw_only_layers:
                for object in self.coils[layer_id]:
                    object.draw_svg(
                        drawing = drawing,
                        parent = svg_layers[layer_id],
                        color = self.config.layers_color.get(layer_id),
                        thickness = self.config.trace_width,
                        dashes = "none",
                    )

        # Draw the link traces
        for link in self.links:
            if self.config.draw_only_layers is None or link.layer in self.config.draw_only_layers:
                link.draw_svg(
                    drawing = drawing,
                    parent = svg_layers[link.layer],
                    color = self.config.layers_color.get(link.layer),
                    thickness = link.trace_width,
                    dashes = "none",
                )

        # Draw the vias
        if self.config.draw_vias:
            for via in self.vias:
                via.draw_svg(
                    drawing = drawing,
                    parent = svg_layers['vias'],
                    via_color = self.config.via_color,
                    hole_color = self.config.via_hole_color,
                    opacity = self.config.via_opacity,
                )

        # Draw the terminals
        # Castellated-hole terminals are clipped for better rendering (only when exporting with the
        # SVG 'full' profile, as clipping paths are not available on the 'tiny' profile).
        clip_path_id = None
        if self.config.terminal_type == TerminalType.CASTELLATED and self.config.svg_profile == 'full':
            clip_path_id = "clip_castellated_terminals"
            clip_path = drawing.defs.add(drawing.clipPath(id=clip_path_id))
            clip_path.add(drawing.circle(self.board_center.to_viewport().as_tuple(), self.config.board_radius))
        if self.config.draw_terminals:
            for terminal in self.terminals:
                terminal.draw_svg(
                    drawing = drawing,
                    parent = svg_layers['terminals'],
                    pad_color = self.config.terminal_color,
                    hole_color = self.config.terminal_hole_color,
                    opacity = self.config.terminal_opacity,
                    clip_path_id = clip_path_id,
                )

        # Draw the outline
        if self.config.draw_outline:
            self.outline.draw_svg(
                drawing = drawing,
                parent = svg_layers['outline'],
                color = self.config.outline_color,
                thickness = self.config.outline_thickness,
                dashes = "none",
            )

        # Draw the construction geometry
        if self.config.draw_construction_geometry:
            for object in self.construction_geometry:
                match object:
                    case svg.base.BaseElement():
                        svg_layers['construction'].add(object)
                    
                    case DrawableObject():
                        object.draw_svg(drawing, svg_layers['construction'])

        # Draw the magnets
        if self.config.draw_magnets:
            for i in range(self.config.n_magnets):
                center = Point.polar(360.0 * i / self.config.n_magnets, self.config.magnets_position_radius)
                svg_layers['magnets'].add(drawing.circle(
                    center = center.to_viewport().as_tuple(),
                    r = self.config.magnets_diameter / 2.0,
                    stroke = self.config.magnets_color,
                    stroke_width = self.config.magnets_thickness,
                    stroke_opacity = self.config.magnets_opacity,
                    stroke_dasharray = self.config.magnets_dashes,
                    label = f"Magnet_{i + 1}"
                ))

        # Add the groups to the final output file
        for layer_id in svg_layers_ids:
            drawing.add(svg_layers[layer_id])

        return self