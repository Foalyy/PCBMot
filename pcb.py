from typing import Self
from enum import Enum
from .config import Config, BoardShape, TerminalType
import svgwrite as svg
import math, json
from .geometry import sin, cos, tan, asin, acos, atan, atan2
from .geometry import DrawableObject, Vector, Point, Line, Segment, Circle, Arc, PathSegment, PathArc, Path
from .kicad import KicadPCB

def calculate_resistance(config: Config, length: float, width: float, thickness: float):
    """Calculate the resistance of a copper trace of the given geometry based on the physical parameters in the config"""
    return (config.copper_resistivity * length / (width * thickness)) * (1 + config.copper_temperature_coefficient * (config.temperature - 20.0))

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

    def __init__(self, center: Point, diameter: float, hole_diameter: float, phase: int = None, tag = None, label: str = None):
        self.center: Point = center
        self.diameter: float = diameter
        self.hole_diameter: float = hole_diameter
        self.phase: int = phase
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
    
    def copy(self) -> Self:
        """Create a copy of this Via"""
        return Via(
            center = self.center,
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            phase = self.phase,
            tag = self.tag,
            label = self.label,
        )
    
    def rotate(self, rotation_center: Self, angle: float):
        """Rotate this Via around the given center point by the given angle"""
        self.center = self.center.rotated(rotation_center, angle)
    
    def rotated(self, rotation_center: Self, angle: float) -> Self:
        """Create a copy of this Via rotated around the given center point by the given angle"""
        via = self.copy()
        via.center = self.center.rotated(rotation_center, angle)
        return via
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Via mirrored about the X axis"""
        via = self.copy()
        via.center = self.center.mirrored_x()
        return via
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Via mirrored about the Y axis"""
        via = self.copy()
        via.center = self.center.mirrored_y()
        return via
    
    def draw_svg(
            self,
            drawing: svg.Drawing,
            parent: svg.base.BaseElement,
            via_color: str,
            hole_color: str,
            opacity: float = None,
        ):
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
    
    def draw_kicad(self, kicadpcb: KicadPCB) -> Self:
        """Draw this Via on the given Kicad board"""
        kicadpcb.via(
            center = self.center,
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
        )
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

    def __init__(
            self,
            center: Point,
            terminal_type: TerminalType,
            diameter: float,
            hole_diameter: float,
            angle: float = None,
            tag=None,
            label: str = None,
            pad_name: str = None,
        ):
        self.center: Point = center
        self.terminal_type: TerminalType = terminal_type
        self.diameter: float = diameter
        self.hole_diameter: float = hole_diameter
        self.angle: float = angle
        self.tag = tag
        self.label = label
        self.pad_name = pad_name
    
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
    
    def copy(self) -> Self:
        """Create a copy of this Terminal"""
        return Terminal(
            center = self.center,
            terminal_type = self.terminal_type,
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            angle = self.angle,
            tag = self.tag,
            label = self.label,
            pad_name = self.pad_name,
        )

    def rotate(self, rotation_center: Self, angle: float):
        """Rotate this Terminal around the given center point by the given angle"""
        self.center = self.center.rotated(rotation_center, angle)
        self.angle += angle

    def rotated(self, rotation_center: Self, angle: float) -> Self:
        """Create a copy of this Terminal rotated around the given center point by the given angle"""
        terminal = self.copy()
        terminal.center = terminal.center.rotated(rotation_center, angle)
        terminal.angle = angle
        return terminal
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Terminal mirrored about the X axis"""
        terminal = self.copy
        terminal.center = terminal.center.mirrored_x()
        terminal.angle += 180
        return terminal
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Terminal mirrored about the Y axis"""
        terminal = self.copy
        terminal.center = terminal.center.mirrored_y()
        terminal.angle = -terminal.angle
        return terminal
    
    def draw_svg(
            self,
            drawing: svg.Drawing,
            parent: svg.base.BaseElement,
            pad_color: str,
            hole_color: str,
            ref_font_family: str,
            ref_color: str,
            ref_font_size_factor: float,
            opacity: float = None,
            clip_path_id = None,
            only_layer = None,
        ):
        """Draw this Terminal on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if self.terminal_type == TerminalType.NONE:
            return

        if only_layer is None or only_layer == 'terminals':
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

        if only_layer is None or only_layer == 'top_silk':
            # Draw the ref
            ref_pos, ref_angle, ref_size = self._ref_pos(- self.diameter * 0.35)
            parent.add(drawing.text(
                text = self.pad_name,
                insert = ref_pos.to_viewport().as_tuple(),
                stroke = "none",
                fill = ref_color,
                text_anchor = "middle",
                font_size = ref_size * ref_font_size_factor,
                font_family = ref_font_family,
                font_weight = "bolder",
                transform = f"rotate({ref_angle}, {ref_pos.to_svg()})",
                label = f"{self.label}_ref",
            ))

        return self
    
    def draw_kicad(self, kicadpcb: KicadPCB) -> Self:
        """Draw this Terminal on the given Kicad board"""
        ref_pos, ref_angle, ref_size = self._ref_pos(- self.diameter * 0.08)
        kicadpcb.terminal(
            center = self.center,
            diameter = self.diameter,
            hole_diameter = self.hole_diameter,
            ref = self.pad_name,
            ref_offset = ref_pos - self.center,
            ref_angle = ref_angle,
            ref_size = ref_size,
        )
        return self
    
    def _ref_pos(self, offset: float) -> tuple[Vector, float, float]:
        board_center = Point.origin()
        ref_size = self.diameter * 0.7
        angle_offset = 180 * ((self.diameter / 2) * 1.4 + len(self.pad_name) * ref_size / 2) / (math.pi * board_center.distance(self.center))
        ref_pos = self.center.rotated(board_center, angle_offset)
        if ref_pos.y < 0:
            offset = -offset
        ref_pos = ref_pos.offset(Line.from_two_points(board_center, ref_pos), offset)
        ref_angle = self.angle + angle_offset
        return ref_pos, ref_angle, ref_size


class Coil:
    """A single coil on the board"""

    def __init__(self, path: Path, trace_width: float, thickness: float, rotation: float, n_turns: int, label: str = None):
        self.path: Path = path
        self.trace_width: float = trace_width
        self.thickness: float = thickness
        self.rotation: float = rotation
        self.n_turns: int = n_turns
        self.label = label

    def generate(
            config: Config,
            layer_id: str,
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
        outer_arc_radius_offset = 0
        if CoilConnection.OUTSIDE_OUTER_RIGHT_VIA in outside_vias or CoilConnection.OUTSIDE_OUTER_LEFT_VIA in outside_vias:
            outer_fillet_radius, outer_arc_radius_offset, success = Coil._compute_fillet(
                arc_radius = outer_radius,
                opposite_radius = inner_radius,
                angle = angle,
                initial_fillet_radius = outer_fillet_radius,
                via = outside_vias[CoilConnection.OUTSIDE_OUTER_RIGHT_VIA],
                trace_width = config.trace_width,
                trace_spacing = config.trace_spacing,
                config = config,
                construction_geometry = construction_geometry,
            )
            if not success:
                print("Warning : unable to compute a fillet for the outer side of the coil, check for collision with the via and try increasing outer_vias_offset")

        # Same calculation for the inner fillet radius
        inner_fillet_radius, inner_arc_radius_offset, success = Coil._compute_fillet(
            arc_radius = inner_radius,
            opposite_radius = outer_radius,
            angle = angle,
            initial_fillet_radius = config.trace_width * 2,
            via = outside_vias[CoilConnection.OUTSIDE_INNER_RIGHT_VIA],
            trace_width = config.trace_width,
            trace_spacing = config.trace_spacing,
            config = config,
            construction_geometry = construction_geometry,
        )
        if not success:
            print("Warning : unable to compute a fillet for the inner side of the coil, check for collision with the via and try increasing inner_vias_offset")

        # Calculate the spiral of the coil
        outer_radius_initial = outer_radius + outer_arc_radius_offset
        inner_radius_initial = inner_radius + inner_arc_radius_offset
        line_left = Line.from_two_points(board_center, Point.polar(-angle/2.0, outer_radius_initial)).offset(loop_offset * 0.5 + outer_fillet_radius)
        arc_outer = Arc(Point.polar(-angle/2.0, outer_radius_initial), Point.polar(angle/2.0, outer_radius_initial), outer_radius_initial)
        point_start = line_left.intersect(arc_outer)
        path = Path(point_start)
        n_turns = 0
        collision_via = None
        triangular = False
        max_turns = 10000 if config.turns_per_layer == 'auto' else config.turns_per_layer
        for i in range(max_turns + 1): # +1 because a part of the first and last turns will be removed later
            # Construction geometry
            line_left = Line.from_two_points(board_center, Point.polar(-angle/2.0, outer_radius)).offset(loop_offset * (i + 0.5))
            line_right = Line.from_two_points(board_center, Point.polar(angle/2.0, outer_radius)).offset(-loop_offset * (i + 0.5))
            line_center = Line.from_two_points(board_center, Point.polar(0, outer_radius))
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

            if triangular:
                # Segment to inner point
                point_inner = line_right.intersect(line_center)
                segment = Segment(path.end_point, point_inner)
                closest_via, distance = Via.closest_in_list(inside_vias.values(), segment)
                if distance < config.trace_width / 2.0 + config.trace_spacing: # Collision
                    collision_via = closest_via
                path.append_segment(point_inner, fillet_radius=outer_fillet_radius, tag=CoilSide.RIGHT)
                if collision_via is not None:
                    break
            else:
                # Segment to inner right
                point_inner_right = line_right.intersect(arc_inner)
                if point_inner_right.x < config.trace_width:
                    triangular = True
                    point_inner = line_right.intersect(line_center)
                    segment = Segment(path.end_point, point_inner)
                    closest_via, distance = Via.closest_in_list(inside_vias.values(), segment)
                    if distance < config.trace_width / 2.0 + config.trace_spacing: # Collision
                        collision_via = closest_via
                    path.append_segment(point_inner, fillet_radius=outer_fillet_radius, tag=CoilSide.RIGHT)
                    if collision_via is not None:
                        break
                else:
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
            fillet_radius = config.trace_width if triangular else inner_fillet_radius
            path.append_segment(point_outer_left, fillet_radius=fillet_radius, tag=CoilSide.LEFT)
            if collision_via is not None:
                break

            # Reduce the fillet radius for the next loop
            min_radius = config.trace_width
            outer_fillet_radius -= loop_offset
            if outer_fillet_radius < min_radius:
                outer_fillet_radius = min_radius
            inner_fillet_radius -= loop_offset
            if inner_fillet_radius < min_radius:
                inner_fillet_radius = min_radius
            
            # Count the number of turns in the coil
            n_turns += 1

        # Connect the start of the coil to the requested via or terminal
        n_removed = 0
        while not CoilConnection.match_side(target_outside_connection, path.first().tag):
            # Pop the first elements from the path until we find the one that should be connected to the target via or terminal
            path.pop_first()
            path.pop_first() # Fillet
            if len(path.elements) == 1:
                raise ValueError("Error : unable to connect the outside of the coil to a via")
            n_removed += 1
        if n_removed >= 2:
            n_turns -= 1
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
                    inside_connection = inside_connection.mirrored_y()
            except KeyError:
                raise ValueError(f"Invalid inside_connection : {inside_connection}")

            if triangular and inside_connection == CoilConnection.INSIDE_INNER_VIA:
                # Special case for triangular traces connecting to the inner via
                path.pop()
                path.pop()
                while path.last().tag != CoilSide.RIGHT:
                    path.pop()
                    if len(path.elements) == 1:
                        raise ValueError("Error : unable to connect the inside of the coil to a via, try increasing hole_diameter or reducing n_coils_per_phase")
                
                # Create a corner at the point of the triangle that connects to the target via
                last_element = path.last_geometry()
                corner = Point(target_via.center.x, last_element.p2.y)
                path.append_segment(corner)
                path.append_segment(target_via.center)

            else:
                # Pop the last elements from the path until we find the one that should be connected to the target via
                while not inside_connection.match_side(path.last().tag):
                    path.pop()
                    if len(path.elements) == 1:
                        raise ValueError("Error : unable to connect the inside of the coil to a via, try increasing hole_diameter or reducing n_coils_per_phase")
                
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

        # If the coil is anticlockwise, mirror the path relative to the Y axis
        if anticlockwise:
            path = path.mirrored_y()
        
        # Get the trace thickness based on the layer
        is_outer_layer: bool = (layer_id == config.copper_layers[0] or layer_id == config.copper_layers[-1])
        thickness = config.outer_layers_copper_thickness if is_outer_layer else config.inner_layers_copper_thickness

        # Return the coil based on this path
        return Coil(
            path = path,
            trace_width = config.trace_width,
            thickness = thickness,
            rotation = 0.0,
            n_turns = n_turns,
        )
    
    def length(self) -> float:
        """Calculate the total length of the path of this coil"""
        length = 0
        p1 = self.path.start_point
        for element in self.path.elements:
            geometry = element.geometry(p1)
            length += geometry.length()
            p1 = element.p2
        return length
    
    def radial_length(self, center: Point) -> float:
        """Calculate the length of this coil projected on its radial direction
        
        This is currently an approximation for arc elements based only on the start and end points, but it is good enough
        for the usual geometries of coils.
        """
        radial_length = 0
        p1 = self.path.start_point
        for element in self.path.elements:
            geometry = element.geometry(p1)
            d1 = center.distance(geometry.p1)
            d2 = center.distance(geometry.p2)
            radial_length += math.fabs(d2 - d1)
            p1 = element.p2
        return radial_length
    
    def resistance(self, config: Config) -> float:
        """Calculate the total resistance of this coil based on its trace width and thickness"""
        resistance = 0
        p1 = self.path.start_point
        for element in self.path.elements:
            geometry = element.geometry(p1)
            resistance += calculate_resistance(config, geometry.length(), self.trace_width, self.thickness)
            p1 = element.p2
        return resistance

    def copy(self) -> Self:
        """Create a copy of this Coil"""
        return Coil(
            path = self.path,
            trace_width = self.trace_width,
            thickness = self.thickness,
            rotation = self.rotation,
            n_turns = self.n_turns,
            label = self.label
        )

    def rotate(self, center: Self, angle: float):
        """Rotate this Coil around the given center point by the given angle"""
        self.path = self.path.rotated(center, angle)
        self.rotation += angle

    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Coil rotated around the given center point by the given angle"""
        coil = self.copy()
        coil.path = coil.path.rotated(center, angle)
        coil.rotation = self.rotation + angle
        return coil
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Coil mirrored about the X axis"""
        coil = self.copy()
        coil.path = coil.path.mirrored_x()
        return coil
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Coil mirrored about the Y axis"""
        coil = self.copy()
        coil.path = coil.path.mirrored_y()
        return coil
    
    def draw_svg(
            self,
            drawing: svg.Drawing,
            parent: svg.base.BaseElement,
            color: str = None,
            opacity: float = None,
            dashes: str = None,
        ):
        """Draw this Coil on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        self.path.draw_svg(
            drawing = drawing,
            parent = parent,
            color = color,
            opacity = opacity,
            line_width = self.trace_width,
            dashes = dashes,
            label = self.label,
        )
        return self
    
    def draw_kicad(self, kicadpcb: KicadPCB, layer: str) -> Self:
        """Draw this Coil on the given Kicad board"""
        kicadpcb.path(
            path = self.path,
            width = self.trace_width,
            layer = layer,
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
            config: Config,
            construction_geometry = None
        ) -> tuple[float, bool]:
        arc_radius_offset = 0
        while arc_radius_offset < via.diameter + max(-config.inner_vias_offset, 0):
            radius = initial_fillet_radius
            arc = Arc(Point.polar(-angle/2.0, (arc_radius + arc_radius_offset)), Point.polar(angle/2.0, (arc_radius + arc_radius_offset)), (arc_radius + arc_radius_offset))
            segment_base = Segment(Point.polar(angle/2.0, opposite_radius), Point.polar(angle/2.0, (arc_radius + arc_radius_offset)))
            segment = segment_base.offset_closest_to(arc.p1, (trace_spacing + trace_width) * 0.5)
            corner = segment.intersect(arc)
            start_point = arc.midpoint()
            end_point = segment.midpoint()
            while radius < (arc_radius + arc_radius_offset):
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
                    # Try a larger arc offset
                    break

                # Check if the end of the fillet is too close to the center line
                if fillet_arc.p2.x < trace_width * 2.0:
                    # Try a larger arc offset
                    break

                # Check if the fillet is large enough to prevent a collision with the via
                if via.distance(fillet_arc) >= trace_width / 2.0 + trace_spacing:
                    # Success : return this fillet radius and arc radius offset
                    return try_fillet_radius, arc_radius_offset, True

                # This radius looks ok, save it and try again with a larger one
                radius = try_fillet_radius
            
            # Unable to find a fillet radius that doesn't collide with the via for this arc radius,
            # add an offset and try again
            if opposite_radius > arc_radius:
                arc_radius_offset += trace_width
            else:
                arc_radius_offset -= trace_width
        
        # Unable to find a fillet radius that doesn't collide with the via for any reasonable arc radius,
        # give up and return the last values tried
        return radius, arc_radius_offset, False

class Link:
    """A trace on the board linking two coils"""

    def __init__(self, path: Path, trace_width: float, thickness: float, layer: str, phase: int = None, label: str = None):
        self.path: Path = path
        self.trace_width: float = trace_width
        self.thickness: float = thickness
        self.layer: str = layer
        self.phase: str = phase
        self.label: str = label
    
    def length(self) -> float:
        """Calculate the length of this Link trace"""
        length = 0
        p1 = self.path.start_point
        for element in self.path.elements:
            geometry = element.geometry(p1)
            length += geometry.length()
            p1 = element.p2
        return length
    
    def resistance(self, config: Config) -> float:
        """Calculate the resistance of this Link trace"""
        resistance = 0
        p1 = self.path.start_point
        for element in self.path.elements:
            geometry = element.geometry(p1)
            resistance += calculate_resistance(config, geometry.length(), self.trace_width, self.thickness)
            p1 = element.p2
        return resistance
    
    def copy(self) -> Self:
        """Create a copy of this Link"""
        return Link(
            path = self.path,
            trace_width = self.trace_width,
            thickness = self.thickness,
            layer = self.layer,
            phase = self.phase,
            label = self.label,
        )

    def rotate(self, center: Self, angle: float):
        """Rotate this Link around the given center point by the given angle"""
        self.path = self.path.rotated(center, angle)

    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Link rotated around the given center point by the given angle"""
        link = self.copy()
        link.path = self.path.rotated(center, angle)
        return link
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Link mirrored about the X axis"""
        link = self.copy()
        link.path = self.path.mirrored_x()
        return link
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Link mirrored about the Y axis"""
        link = self.copy()
        link.path = self.path.mirrored_y()
        return link
    
    def draw_svg(
            self,
            drawing: svg.Drawing,
            parent: svg.base.BaseElement,
            color: str = None,
            opacity: float = None,
            dashes: str = None,
        ):
        """Draw this Link on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        self.path.draw_svg(
            drawing = drawing,
            parent = parent,
            color = color,
            opacity = opacity,
            line_width = self.trace_width,
            dashes = dashes,
            label = self.label,
        )
        return self
    
    def draw_kicad(self, kicadpcb: KicadPCB) -> Self:
        """Draw this Link on the given Kicad board"""
        kicadpcb.path(
            path = self.path,
            width = self.trace_width,
            layer = self.layer,
        )
        return self

class SilkscreenText:
    """A text printed on the silkscreen layer"""

    def __init__(
            self,
            text: str,
            position: Point,
            size: float,
            rotation: float,
            layer: str,
            font_family: str = None,
            label: str = None
        ):
        self.text: str = text
        self.position: Point = position
        self.size: float = size
        self.rotation: float = rotation
        self.layer: str = layer
        self.font_family: float = font_family
        self.label: float = label

    def rotate(self, center: Self, angle: float):
        """Rotate this text around the given center point by the given angle"""
        self.position = self.position.rotated(center, angle)
        self.rotation += angle
    
    def draw_svg(
        self,
        drawing: svg.Drawing,
        parent: svg.base.BaseElement,
        color : str = None,
    ):
        """Draw this text on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        parent.add(drawing.text(
            text = self.text,
            insert = self.position.to_viewport().as_tuple(),
            stroke = "none",
            fill = color,
            text_anchor = "middle",
            font_size = self.size,
            font_family = self.font_family,
            font_weight = "bolder",
            transform = f"rotate({self.rotation}, {self.position.to_svg()})",
            label = self.label,
        ))
        return self
    
    def draw_kicad(self, kicadpcb: KicadPCB) -> Self:
        """Draw this text on the given Kicad board"""
        kicadpcb.text(
            text = self.text,
            center = self.position,
            angle = self.rotation,
            font_size = self.size,
            layer = self.layer,
        )
        return self

class Outline:
    """The shape and dimensions of the board"""

    def __init__(self, board_center: Point, config: Config):
        self.board_center = board_center
        self.board_shape = config.board_shape
        self.board_diameter = config.board_diameter
        self.hole_diameter = config.hole_diameter
        self.board_chamfer = config.board_chamfer
        self.board_fillet = config.board_fillet
        if self.board_chamfer is None:
            self.board_chamfer = 0
        self.board_chamfer = min(self.board_chamfer, self.board_diameter / 2.0)
        self.generate_mountpoints = config.generate_mountpoints
        self.n_mountpoints = config.n_mountpoints
        self.mountpoints_position_radius = config.mountpoints_position_radius
        self.mountpoints_diameter = config.mountpoints_diameter
        self.mountpoints_marking_diameter = config.mountpoints_marking_diameter

        # Calculate the shape of a square board with chamfers and fillets
        if self.board_shape == BoardShape.SQUARE:
            r = self.board_diameter / 2.0
            self.path = Path(Point(0, r))
            self.path.append_segment(Point(r - self.board_chamfer, r))
            if self.board_chamfer > 0:
                self.path.append_segment(Point(r, r - self.board_chamfer), fillet_radius = self.board_fillet)
            self.path.append_segment(Point(r, -(r - self.board_chamfer)), fillet_radius = self.board_fillet)
            if self.board_chamfer > 0:
                self.path.append_segment(Point(r - self.board_chamfer, -r), fillet_radius = self.board_fillet)
            self.path.append_segment(Point(-(r - self.board_chamfer), -r), fillet_radius = self.board_fillet)
            if self.board_chamfer > 0:
                self.path.append_segment(Point(-r, -(r - self.board_chamfer)), fillet_radius = self.board_fillet)
            self.path.append_segment(Point(-r, r - self.board_chamfer), fillet_radius = self.board_fillet)
            if self.board_chamfer > 0:
                self.path.append_segment(Point(-(r - self.board_chamfer), r), fillet_radius = self.board_fillet)
            self.path.append_segment(Point(0, r), fillet_radius = self.board_fillet)
    
    def draw_svg(
            self,
            drawing: svg.Drawing,
            parent: svg.base.BaseElement,
            color: str = None,
            line_width: float = None,
            dashes: str = None,
            marking_color: str = None,
            marking_line_width: float = None,
        ):
        """Draw this Outline on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        # Draw the outer edge
        match self.board_shape:
            case BoardShape.CIRCLE:
                parent.add(drawing.circle(
                    self.board_center.to_viewport().as_tuple(),
                    self.board_diameter / 2.0,
                    stroke = color,
                    stroke_width = line_width,
                    stroke_dasharray = dashes,
                    fill = "none",
                    label = f"Outline_edge",
                ))
            case BoardShape.SQUARE:
                self.path.draw_svg(
                    drawing,
                    parent,
                    color = color,
                    opacity = 1.0,
                    line_width = line_width,
                    dashes = dashes,
                    label = f"Outline_edge",
                )
        
        # Draw the inner hole
        parent.add(drawing.circle(
            self.board_center.to_viewport().as_tuple(),
            self.hole_diameter / 2.0,
            stroke = color,
            stroke_width = line_width,
            stroke_dasharray = dashes,
            fill = "none",
            label = f"Outline_hole",
        ))

        # Draw the mountpoints
        if self.generate_mountpoints and self.n_mountpoints is not None:
            for i in range(self.n_mountpoints):
                center = Point.polar((i + 0.5) * 360 / self.n_mountpoints, self.mountpoints_position_radius)
                parent.add(drawing.circle(
                    center.to_viewport().as_tuple(),
                    self.mountpoints_diameter / 2.0,
                    stroke = color,
                    stroke_width = line_width,
                    stroke_dasharray = dashes,
                    fill = "none",
                    label = f"Mountpoint_{i}",
                ))
                if self.mountpoints_marking_diameter > 0:
                    parent.add(drawing.circle(
                        center.to_viewport().as_tuple(),
                        self.mountpoints_marking_diameter / 2.0,
                        stroke = marking_color,
                        stroke_width = marking_line_width,
                        fill = "none",
                        label = f"Mountpoint_marking_{i}",
                    ))
        
        return self
    
    def draw_kicad(self, kicadpcb: KicadPCB, width: float, marking_width: float) -> Self:
        """Draw this Outline on the given Kicad board"""

        # Draw the outer edge
        match self.board_shape:
            case BoardShape.CIRCLE:
                kicadpcb.gr_circle(
                    center = self.board_center,
                    radius = self.board_diameter / 2.0,
                    width = width,
                    layer = 'outline',
                )
            case BoardShape.SQUARE:
                kicadpcb.gr_path(
                    path = self.path,
                    width = width,
                    layer = 'outline',
                )

        # Draw the inner hole
        kicadpcb.gr_circle(
            center = self.board_center,
            radius = self.hole_diameter / 2.0,
            width = width,
            layer = 'outline',
        )

        # Draw the mountpoints
        if self.generate_mountpoints and self.n_mountpoints is not None:
            for i in range(self.n_mountpoints):
                center = Point.polar((i + 0.5) * 360 / self.n_mountpoints, self.mountpoints_position_radius)
                kicadpcb.gr_circle(
                    center = center,
                    radius = self.mountpoints_diameter / 2.0,
                    width = width,
                    layer = 'outline',
                )
                if self.mountpoints_marking_diameter > 0:
                    kicadpcb.gr_circle(
                        center = center,
                        radius = self.mountpoints_marking_diameter / 2.0,
                        width = marking_width,
                        layer = 'top_silk',
                    )
                    kicadpcb.gr_circle(
                        center = center,
                        radius = self.mountpoints_marking_diameter / 2.0,
                        width = marking_width,
                        layer = 'bottom_silk',
                    )

        return self

class PCBStats:
    """A storage class for a list of computed stats about a PCB"""
    def __init__(
        self,
        layers: list[str],
        coil_turns: int,
        coil_length: float,
        coil_radial_length: float,
        coil_resistance: float,
        phases_length: dict[int, float],
        phases_resistance: dict[int, float],
        n_magnets: int,
        magnets_diameter: float,
    ):
        self.layers = layers
        self.coil_turns = coil_turns
        self.coil_length = coil_length
        self.coil_radial_length = coil_radial_length
        self.coil_resistance = coil_resistance
        self.phases_length = phases_length
        self.phases_resistance = phases_resistance
        self.n_magnets = n_magnets
        self.magnets_diameter = magnets_diameter

        self.coil_radial_length_ratio = coil_radial_length / coil_length
        average_phase_length = 0
        average_phase_resistance = 0
        for phase in phases_length:
            average_phase_length += phases_length[phase]
        self.average_phase_length = average_phase_length / len(phases_length)
        for phase in phases_resistance:
            average_phase_resistance += phases_resistance[phase]
        self.average_phase_resistance = average_phase_resistance / len(phases_resistance)
    
    def json(self) -> str:
        """Return the stats as a JSON-formatted string"""
        return json.dumps({
            'layers': self.layers,
            'coil_turns': self.coil_turns,
            'coil_length': round(self.coil_length, 1),
            'coil_radial_length': round(self.coil_radial_length, 1),
            'coil_radial_length_ratio': round(self.coil_radial_length_ratio, 3),
            'coil_resistance': round(self.coil_resistance, 3),
            'phases_length': self.phases_length,
            'average_phase_length': round(self.average_phase_length, 1),
            'phases_resistance': self.phases_resistance,
            'average_phase_resistance': round(self.average_phase_resistance, 3),
            'n_magnets': self.n_magnets,
            'magnets_diameter': round(self.magnets_diameter, 1),
        }, indent=4)
    
    def write(self, destination: str):
        """Write the stats to the given file, or the standard output if the destination is '-'"""
        data = self.json()
        if destination == '-':
            print(data)
        else:
            with open(destination, 'w') as f:
                f.write(data)


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
        top_silk: list[SilkscreenText],
        bottom_silk: list[SilkscreenText],
        construction_geometry: list,
        layers: list[str],
        stats: PCBStats,
    ):
        self.config = config
        self.board_center = board_center
        self.outline = outline
        self.coils = coils
        self.vias = vias
        self.terminals = terminals
        self.links = links
        self.top_silk = top_silk
        self.bottom_silk = bottom_silk
        self.construction_geometry = construction_geometry
        self.layers = layers
        self.stats = stats
    
    def generate(config: Config, compute_stats: bool):
        """Generate a new PCB based on the given config"""
        
        coils = {}
        vias = []
        terminals = []
        links = []
        top_silk = []
        bottom_silk = []
        construction_geometry = []

        # All layers in the final drawing
        layers = ['bottom_silk'] + list(reversed(config.copper_layers)) + ['top_silk', 'outline', 'vias', 'terminals', 'magnets', 'construction']

        # Names of the coils : A1, A2, A3, B1, ...
        coil_names = []
        for slot in range(config.n_coils_per_phase):
            for phase in range(config.n_phases):
                coil_names.append(f"{chr(ord('A') + phase)}{slot + 1}")

        # Center the board on the origin
        board_center = Point.origin()

        # Board outline
        outline = Outline(board_center, config)

        # Construction geometry
        construction_geometry = [
            # Base coil center line
            Segment(
                p1 = board_center,
                p2 = Point(0, config.board_radius),
            ),
            # Base coil left line
            Segment(
                p1 = board_center,
                p2 = Point.polar(-config.coil_angle/2.0, config.board_radius),
            ),
            # Base coil right line
            Segment(
                p1 = board_center,
                p2 = Point.polar(config.coil_angle/2.0, config.board_radius),
            ),
            # Outer circle
            Circle(
                center = board_center,
                radius = config.coils_outer_radius,
            ),
            # Inner circle
            Circle(
                center = board_center,
                radius = config.coils_inner_radius,
            ),
        ]
        if config.board_shape != BoardShape.CIRCLE:
            construction_geometry.append(
                Circle(
                    center = board_center,
                    radius = config.board_radius,
                )
            )

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
                if config.four_layers_inner_vias:
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
        for i, layer_id in enumerate(config.copper_layers):
            layers_specs[layer_id] = {
                'anticlockwise': i % 2 == 1,
                'outside_connection': outside_connections[int((i + 1) / 2)],
                'inside_connection': inside_connections[int(i / 2)],
            }

        # Inside vias
        inside_vias = {}
        inside_vias_center = Point(0, config.coils_middle_radius + config.inside_vias_offset)
        tangential_length = 2 * math.pi * config.coils_middle_radius * config.coil_angle / 360.
        radial_length = config.coils_outer_radius - config.coils_inner_radius
        match config.n_layers:
            case 2:
                # A single via placed in the middle
                via = Via(inside_vias_center, config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                inside_vias_list = [via]
            case 4:
                # Two vias placed vertically
                inside_outer_via = Via(inside_vias_center + Vector(0, config.via_diameter_w_spacing / 2.0), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                inside_inner_via = Via(inside_vias_center - Vector(0, config.via_diameter_w_spacing / 2.0), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
                inside_vias_list = [inside_outer_via, inside_inner_via]
                pass
            case 6:
                # Three vias, placed either vertically, or on a triangular shape
                if tangential_length >= radial_length:
                    # The coil is wider than it is tall, place the via in a triangle
                    inside_inner_via = Via(inside_vias_center + Vector(0, -config.via_diameter_w_spacing / 2.0), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
                    inside_outer_left_via = Via(inside_inner_via.center + Vector.polar(-30, config.via_diameter_w_spacing), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                    inside_outer_right_via = Via(inside_inner_via.center + Vector.polar(30, config.via_diameter_w_spacing), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_RIGHT_VIA)
                    inside_vias_list = [inside_outer_left_via, inside_outer_right_via, inside_inner_via]
                else:
                    # The coil is taller than it is wide, place the via vertically
                    inside_outer_via = Via(inside_vias_center + Vector(0, config.via_diameter_w_spacing), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                    inside_middle_via = Via(inside_vias_center, config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_RIGHT_VIA)
                    inside_inner_via = Via(inside_vias_center + Vector(0, -config.via_diameter_w_spacing), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
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
                        inside_outer_via = Via(inside_vias_center + Vector(0, (config.via_diameter_w_spacing + sep_test) / 2.0), config.via_diameter, config.via_drill_diameter)
                        inside_inner_via = Via(inside_vias_center + Vector(0, -(config.via_diameter_w_spacing + sep_test) / 2.0), config.via_diameter, config.via_drill_diameter)
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
                        inside_via_3 = Via(points[0], config.via_diameter, config.via_drill_diameter)
                        inside_via_4 = Via(points[1], config.via_diameter, config.via_drill_diameter)
                        if inside_via_3.distance(inside_via_4) < config.trace_spacing:
                            break
                        sep = sep_test
                inside_outer_via = Via(inside_vias_center + Vector(0, (config.via_diameter_w_spacing + sep) / 2.0), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_OUTER_VIA)
                inside_inner_via = Via(inside_vias_center - Vector(0, (config.via_diameter_w_spacing + sep) / 2.0), config.via_diameter, config.via_drill_diameter, tag=CoilConnection.INSIDE_INNER_VIA)
                c1 = Circle(inside_outer_via.center, config.via_diameter_w_spacing)
                c2 = Circle(inside_inner_via.center, config.via_diameter_w_spacing)
                points = c1.intersect(c2)
                if points[0].x < points[1].x:
                    tag3, tag4 = CoilConnection.INSIDE_LEFT_VIA, CoilConnection.INSIDE_RIGHT_VIA
                else:
                    tag3, tag4 = CoilConnection.INSIDE_RIGHT_VIA, CoilConnection.INSIDE_LEFT_VIA
                inside_via_3 = Via(points[0], config.via_diameter, config.via_drill_diameter, tag=tag3)
                inside_via_4 = Via(points[1], config.via_diameter, config.via_drill_diameter, tag=tag4)
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
            Point.polar(-config.coil_angle/2.0, config.coils_outer_radius),
            Point.polar(config.coil_angle/2.0, config.coils_outer_radius),
            config.coils_outer_radius
        ).offset(-config.via_diameter / 2.0 + config.outer_vias_offset)
        inner_arc = Arc(
            Point.polar(-config.coil_angle/2.0, config.coils_inner_radius),
            Point.polar(config.coil_angle/2.0, config.coils_inner_radius),
            config.coils_inner_radius
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
                outside_vias[connection] = Via(point, config.via_diameter, config.via_drill_diameter, tag=connection)
        if CoilConnection.OUTSIDE_INNER_LEFT_VIA in outside_vias \
                and outside_vias[CoilConnection.OUTSIDE_INNER_LEFT_VIA].distance(outside_vias[CoilConnection.OUTSIDE_INNER_RIGHT_VIA]) < 0:
            print("Warning : collision between the vias closer to the center of the board")

        # Copy the vias for all coils
        for slot in range(config.n_coils_per_phase):
            for phase in range(config.n_phases):
                coil_idx = slot * config.n_phases + phase
                for connection, via in inside_vias.items():
                    if slot % 2 == 1:
                        via = via.mirrored_y()
                    via = via.rotated(board_center, 360.0 * coil_idx / config.n_coils)
                    via.phase = phase
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
                    via.phase = phase
                    via.label = f"Via_{coil_names[coil_idx]}_{CoilConnection.label(connection)}"
                    vias.append(via)
        
        # Terminal
        terminal_base = None
        if config.terminal_type != TerminalType.NONE:
            match config.terminal_type:
                case TerminalType.THROUGH_HOLE | TerminalType.SMD:
                    # Align the terminal outside the coil
                    y = config.coils_outer_radius + config.trace_spacing + config.terminal_offset + config.terminal_diameter / 2.0
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
                pad_name = coil_names[i]
                if config.link_series_coils and config.link_com:
                    # Only keep the letter, not the number
                    pad_name = pad_name[0]
                terminal.pad_name = pad_name
                terminals.append(terminal)
            if config.link_series_coils:
                if config.link_com:
                    if config.generate_com_terminal:
                        terminal = terminal_base.rotated(board_center, - 360.0 / config.n_coils)
                        terminal.label = f"Terminal_COM"
                        terminal.pad_name = f"COM"
                        terminals.append(terminal)
                else:
                    if config.n_coils_per_phase % 2 == 0:
                        for i in range(config.n_phases):
                            terminal = terminal_base.rotated(board_center, 360.0 * -(i + 1) / config.n_coils)
                            terminal.label = f"Terminal_{coil_names[-(i + 1)]}"
                            pad_name = coil_names[-(i + 1)]
                            terminal.pad_name = pad_name
                            terminals.append(terminal)

        # Generate the coils on all layers
        coil_outside_connection_length = config.trace_width * 2
        for layer_id in config.copper_layers:
            specs = layers_specs[layer_id]

            # Generate the base Coil for this layer
            coil_even = Coil.generate(
                config = config,
                layer_id = layer_id,
                board_center = board_center,
                angle = config.coil_angle,
                outer_radius = config.coils_outer_radius - config.trace_width / 2.0,
                inner_radius = config.coils_inner_radius + config.trace_width / 2.0,
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
                    layer_id = layer_id,
                    board_center = board_center,
                    angle = config.coil_angle,
                    outer_radius = config.coils_outer_radius - config.trace_width / 2.0,
                    inner_radius = config.coils_inner_radius + config.trace_width / 2.0,
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
            for slot in range(config.n_coils_per_phase):
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
        if config.link_series_coils and config.n_coils_per_phase >= 2:
            if config.n_layers >= config.n_phases:
                # Base path for inner connections (even coils)
                # Draw a path offset toward the inner side of the board that connects these two points for the first coil
                connection_point_1 = coils[config.copper_layers[-1]][0].path.start_point
                connection_point_2 = connection_point_1
                if not optimise_inner_vias:
                    connection_point_2 = connection_point_2.mirrored_y()
                connection_point_2 = connection_point_2.rotated(board_center, (360.0 / config.n_coils) * config.n_phases)
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
                    if fillet_radius < config.series_link_inner_trace_width:
                        fillet_radius = None
                    base_path_inner.append_segment(corner1)
                    base_path_inner.append_arc(corner2, trace_center_radius, anticlockwise=False, fillet_radius=fillet_radius)
                    base_path_inner.append_segment(connection_point_2, fillet_radius=fillet_radius)

                # Base path for outer connections (odd coils)
                # Draw a path offset toward the outer side of the board that connects these two points for the first coil
                if config.n_coils_per_phase >= 3:
                    connection_point_1 = coils[config.copper_layers[0]][0].path.start_point
                    connection_point_2 = connection_point_1.mirrored_y().rotated(board_center, (360.0 / config.n_coils) * config.n_phases)
                    outer_vias_radius = config.coils_outer_radius + config.outer_vias_offset + config.trace_spacing + config.series_link_outer_trace_width / 2.0
                    vias_circle_radius = config.coils_outer_radius + max(config.via_diameter / 2.0, config.series_link_outer_trace_width / 2.0) + config.trace_spacing
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
                    if fillet_radius < config.series_link_outer_trace_width:
                        fillet_radius = None
                    base_path_outer = Path(via_pos_1)
                    base_path_outer.append_segment(corner1)
                    base_path_outer.append_arc(corner2, trace_center_radius, anticlockwise=False, fillet_radius=fillet_radius)
                    base_path_outer.append_segment(via_pos_2, fillet_radius=fillet_radius)

                # Create the links on different layers by rotating the base path accordingly for each phase
                for slot_pair in range(config.n_coils_per_phase - 1):
                    for phase in range(config.n_phases):
                        angle = ((slot_pair * config.n_phases) + phase) * (360.0 / config.n_coils)
                        label = f"Link_{coil_names[slot_pair * config.n_phases + phase]}_{coil_names[(slot_pair + 1) * config.n_phases + phase]}"
                        layer_id = config.copper_layers[phase]
                        is_outer_layer = (layer_id == config.copper_layers[0] or layer_id == config.copper_layers[-1])
                        thickness = config.outer_layers_copper_thickness if is_outer_layer else config.inner_layers_copper_thickness
                        if slot_pair % 2 == 0:
                            path = base_path_inner.rotated(board_center, angle)
                            link = Link(
                                path = path,
                                trace_width = config.series_link_inner_trace_width,
                                thickness = thickness,
                                layer = layer_id,
                                phase = phase,
                                label = label,
                            )
                            links.append(link)
                        else:
                            path = base_path_outer.rotated(board_center, angle)
                            link = Link(
                                path = path,
                                trace_width = config.series_link_outer_trace_width,
                                thickness = thickness,
                                layer = layer_id,
                                phase = phase,
                                label = label,
                            )
                            links.append(link)
                            via1_label = f"Link_via_{coil_names[slot_pair * config.n_phases + phase]}"
                            via1 = Via(via_pos_1, config.via_diameter, config.via_drill_diameter, phase=phase, label=via1_label).rotated(board_center, angle)
                            via2_label = f"Link_via_{coil_names[(slot_pair + 1) * config.n_phases + phase]}"
                            via2 = Via(via_pos_2, config.via_diameter, config.via_drill_diameter, phase=phase, label=via2_label).rotated(board_center, angle)
                            link_vias.extend([via1, via2])
            else:
                print("Warning : unable to link the coils, not enough layers")
        vias.extend(link_vias)

        # Connect the common point in a wye (star) configuration, on the top layer
        link_com_connection_point = None
        if config.link_com:
            layer_id = config.copper_layers[0]
            is_outer_layer = (layer_id == config.copper_layers[0] or layer_id == config.copper_layers[-1])
            thickness = config.outer_layers_copper_thickness if is_outer_layer else config.inner_layers_copper_thickness
            if config.n_coils_per_phase % 2 == 0:
                # Even number of coils : the COM link is on the outer side of the board
                connection_point_1 = coils[layer_id][0].path.start_point
                connection_point_2 = connection_point_1.rotated(board_center, -(360.0 / config.n_coils))
                terminal_circle_radius = board_center.distance(terminal_base.center) if terminal_base is not None else 0
                intermediate_circle_radius = config.coils_outer_radius + max(config.com_link_trace_width / 2.0 + config.trace_spacing, coil_outside_connection_length)
                outer_vias_radius = config.coils_outer_radius + config.outer_vias_offset + config.trace_spacing + config.com_link_trace_width / 2.0
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
                link_com_connection_point = intermediate_point_1
            else:
                # Odd number of coils : the COM link is on the inner side of the board
                connection_point_1 = coils[config.copper_layers[-1]][0].path.start_point
                connection_point_2 = connection_point_1.rotated(board_center, -(360.0 / config.n_coils))
                vias_circle_radius = board_center.distance(connection_point_1)
                trace_center_radius = vias_circle_radius - config.via_diameter / 2.0 - config.trace_spacing - config.com_link_trace_width / 2.0 - config.com_link_offset
                base_path_com = Path(connection_point_1)
                if math.isclose(trace_center_radius, vias_circle_radius):
                    base_path_com.append_arc(connection_point_2, trace_center_radius, anticlockwise=False)
                else:
                    circle = Circle(board_center, trace_center_radius)
                    line1 = Line.from_two_points(board_center, connection_point_1)
                    line2 = Line.from_two_points(board_center, connection_point_2)
                    corner1 = connection_point_1.closest(line1.intersect(circle))
                    corner2 = connection_point_2.closest(line2.intersect(circle))
                    fillet_radius = connection_point_1.distance(corner1) / 2.0
                    if fillet_radius < config.com_link_trace_width:
                        fillet_radius = None
                    base_path_com.append_segment(corner1)
                    base_path_com.append_arc(corner2, trace_center_radius, anticlockwise=True, fillet_radius=fillet_radius)
                    base_path_com.append_segment(connection_point_2, fillet_radius=fillet_radius)
            for i in range(config.n_phases - 1):
                path = base_path_com.rotated(board_center, -(i + 1) * (360.0 / config.n_coils))
                label = f"COM_{coil_names[-(i+2)]}_{coil_names[-(i+1)]}"
                link = Link(
                    path = path,
                    trace_width = config.com_link_trace_width,
                    thickness = thickness,
                    layer = layer_id,
                    label = label,
                )
                links.append(link)

        # Connect the terminals and linking vias on the first layer
        layer_id = config.copper_layers[0]
        if config.draw_only_layers is None or layer_id in config.draw_only_layers:
            for i in range(config.n_coils):
                # Connect the start of the coil to the terminal in the following cases :
                if terminal_base and (
                        i < config.n_phases # - the first coil of each phase (always)
                        or not config.link_series_coils # - every coil when not linking coils in series
                        or (not config.link_com and i >= config.n_coils - config.n_phases) # - the last coil of each phase when not linking the COM point
                        or (config.n_coils_per_phase % 2 == 0 and config.generate_com_terminal and i == config.n_coils - 1) # - the single last coil when linking the COM point on the outer side of the board
                    ):
                    coils[layer_id][i].path.prepend_segment(terminal_base.rotated(board_center, 360.0 * i / config.n_coils).center)
                elif config.link_series_coils and link_vias and config.n_coils_per_phase >= 3 and i >= config.n_phases and (config.n_coils_per_phase % 2 != 0 or i < (config.n_coils_per_phase - 1) * config.n_phases):
                    # When linking series coils, connect the start of the coil to the outer linking via for every coil starting at the second coil of each phase,
                    # except for the last coil of each phase if the number of coils per phase is odd
                    radius = board_center.distance(link_vias[0].center)
                    coils[layer_id][i].path.prepend_segment(Point.polar(0, radius).rotated(board_center, 360.0 * i / config.n_coils))
                elif link_com_connection_point and config.n_coils_per_phase % 2 == 0 and config.link_com and i >= config.n_coils - config.n_phases:
                    # When linking the COM point on the outside of the board, connect the start of the coil to the COM link
                    coils[layer_id][i].path.prepend_segment(link_com_connection_point.rotated(board_center, 360.0 * i / config.n_coils))

        # Print the names of the coils on the top silkscreen
        if config.draw_coil_names:
            for i in range(config.n_coils):
                angle = 360. * i / config.n_coils
                top_silk.append(SilkscreenText(
                    text = coil_names[i],
                    position = Point.polar(angle, config.coil_names_position_radius),
                    size = config.coil_names_font_size * config.svg_font_size_factor,
                    rotation = angle,
                    layer = 'top_silk',
                    font_family = config.silk_font_family,
                    label = f"Coil_name_{coil_names[i]}",
                ))

        # Apply the final rotation on everything on the board except the outline and mountpoints
        if config.rotation:
            for layer in coils:
                for coil in coils[layer]:
                    coil.rotate(board_center, config.rotation)
            for via in vias:
                via.rotate(board_center, config.rotation)
            for terminal in terminals:
                terminal.rotate(board_center, config.rotation)
            for link in links:
                link.rotate(board_center, config.rotation)
            for element in top_silk:
                element.rotate(board_center, config.rotation)
            for element in bottom_silk:
                element.rotate(board_center, config.rotation)
            construction_geometry_rotated = []
            for element in construction_geometry:
                construction_geometry_rotated.append(element.rotated(board_center, config.rotation))
            construction_geometry = construction_geometry_rotated
        
        # Compute the stats
        if compute_stats:
            coil_turns = 0
            coil_length = 0
            coil_radial_length = 0
            coil_resistance = 0
            for layer_id in config.copper_layers:
                coil = coils[layer_id][0]
                coil_turns += coil.n_turns
                length = coil.length()
                resistance = coil.resistance(config)
                coil_length += length
                coil_radial_length += coil.radial_length(board_center)
                coil_resistance += resistance
            phases_length = {}
            phases_resistance = {}
            if config.link_series_coils:
                for phase in range(config.n_phases):
                    length = coil_length * config.n_coils_per_phase
                    resistance = coil_resistance * config.n_coils_per_phase
                    for via in vias:
                        # This includes links vias
                        if via.phase == phase:
                            resistance += config.via_resistance
                    for link in links:
                        if link.phase == phase:
                            length += link.length()
                            resistance += link.resistance(config)
                    phases_length[phase] = length
                    phases_resistance[phase] = resistance
            else:
                for phase in range(config.n_phases):
                    phases_length[phase] = coil_length
                    phases_resistance[phase] = coil_resistance + config.n_layers * config.via_resistance
            stats = PCBStats(
                layers = layers,
                coil_turns = coil_turns,
                coil_length = coil_length,
                coil_radial_length = coil_radial_length,
                coil_resistance = coil_resistance,
                phases_length = phases_length,
                phases_resistance = phases_resistance,
                n_magnets = config.n_magnets,
                magnets_diameter = config.magnets_diameter,
            )
        else:
            stats = None

        # Create the PCB
        return PCB(
            config = config,
            board_center = board_center,
            outline = outline,
            coils = coils,
            vias = vias,
            terminals = terminals,
            links = links,
            top_silk = top_silk,
            bottom_silk = bottom_silk,
            construction_geometry = construction_geometry,
            layers = layers,
            stats = stats,
        )
    
    def _draw_layer(self, layer_id: str, only_layer: str) -> bool:
        return (self.config.draw_only_layers is None or layer_id in self.config.draw_only_layers) and (only_layer is None or only_layer == layer_id)
    
    def draw_svg(self, drawing: svg.Drawing, only_layer: str = None) -> Self:
        """Draw this PCB on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        # Create the SVG layers to organize the drawing.
        # Layer groups created later will be displayed on top of the previous ones, so the board layers are reversed
        # to have the top layer on top and bottom layer on the bottom, as expected.
        svg_layers_ids = self.layers
        svg_layers = {}
        for layer_id in svg_layers_ids:
            if only_layer is None or only_layer == layer_id:
                svg_layers[layer_id] = drawing.layer(label=layer_id.capitalize())

        # Draw the coils
        for layer_id in self.config.copper_layers:
            if self._draw_layer(layer_id, only_layer):
                for object in self.coils[layer_id]:
                    object.draw_svg(
                        drawing = drawing,
                        parent = svg_layers[layer_id],
                        color = self.config.copper_layers_color.get(layer_id),
                        dashes = "none",
                    )

        # Draw the link traces
        for link in self.links:
            if self._draw_layer(link.layer, only_layer):
                link.draw_svg(
                    drawing = drawing,
                    parent = svg_layers[link.layer],
                    color = self.config.copper_layers_color.get(link.layer),
                    dashes = "none",
                )

        # Draw the vias
        if self.config.draw_vias and (only_layer is None or only_layer == 'vias'):
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
        if self.config.draw_terminals and (only_layer is None or only_layer == 'terminals'):
            clip_path_id = None
            if self.config.terminal_type == TerminalType.CASTELLATED and self.config.svg_profile == 'full':
                clip_path_id = "clip_castellated_terminals"
                clip_path = drawing.defs.add(drawing.clipPath(id=clip_path_id))
                clip_path.add(drawing.circle(self.board_center.to_viewport().as_tuple(), self.config.board_radius))
            for terminal in self.terminals:
                terminal.draw_svg(
                    drawing = drawing,
                    parent = svg_layers['terminals'],
                    pad_color = self.config.terminal_color,
                    hole_color = self.config.terminal_hole_color,
                    ref_font_family = self.config.silk_font_family,
                    ref_color = self.config.top_silk_color,
                    ref_font_size_factor = self.config.svg_font_size_factor,
                    opacity = self.config.terminal_opacity,
                    clip_path_id = clip_path_id,
                    only_layer = 'terminals',
                )
        
        # Draw the top and bottom silkscreens
        if only_layer is None or only_layer == 'top_silk':
            for item in self.top_silk:
                item.draw_svg(
                    drawing = drawing,
                    parent = svg_layers['top_silk'],
                    color = self.config.top_silk_color,
                )
            for terminal in self.terminals:
                terminal.draw_svg(
                    drawing = drawing,
                    parent = svg_layers['top_silk'],
                    pad_color = self.config.terminal_color,
                    hole_color = self.config.terminal_hole_color,
                    ref_font_family = self.config.silk_font_family,
                    ref_color = self.config.top_silk_color,
                    ref_font_size_factor = self.config.svg_font_size_factor,
                    opacity = self.config.terminal_opacity,
                    only_layer = 'top_silk',
                )
        if only_layer is None or only_layer == 'bottom_silk':
            for item in self.bottom_silk:
                item.draw_svg(
                    drawing = drawing,
                    parent = svg_layers['bottom_silk'],
                    color = self.config.bottom_silk_color,
                )

        # Draw the outline
        if self.config.draw_outline and (only_layer is None or only_layer == 'outline'):
            self.outline.draw_svg(
                drawing = drawing,
                parent = svg_layers['outline'],
                color = self.config.outline_color,
                line_width = self.config.outline_line_width,
                dashes = "none",
                marking_color = self.config.top_silk_color,
                marking_line_width = self.config.silk_line_width,
            )

        # Draw the construction geometry
        if self.config.draw_construction_geometry and (only_layer is None or only_layer == 'construction'):
            for object in self.construction_geometry:
                match object:
                    case svg.base.BaseElement():
                        svg_layers['construction'].add(object)
                    
                    case DrawableObject():
                        object.draw_svg(
                            drawing = drawing,
                            parent = svg_layers['construction'],
                            color = self.config.construction_geometry_color,
                            line_width = self.config.construction_geometry_line_width,
                            dashes = self.config.construction_geometry_dashes,
                        )

        # Draw the magnets
        if self.config.draw_magnets and (only_layer is None or only_layer == 'magnets'):
            for i in range(self.config.n_magnets):
                center = Point.polar(360.0 * i / self.config.n_magnets + self.config.rotation, self.config.magnets_position_radius)
                svg_layers['magnets'].add(drawing.circle(
                    center = center.to_viewport().as_tuple(),
                    r = self.config.magnets_diameter / 2.0,
                    stroke = self.config.magnets_color,
                    stroke_width = self.config.magnets_line_width,
                    stroke_opacity = self.config.magnets_opacity,
                    stroke_dasharray = self.config.magnets_dashes,
                    label = f"Magnet_{i + 1}"
                ))

        # Add the groups to the final output file
        for layer_id in svg_layers_ids:
            if only_layer is None or only_layer == layer_id:
                drawing.add(svg_layers[layer_id])

        return self
    
    def draw_kicad(self, kicadpcb: KicadPCB) -> Self:
        """Draw this PCB on the given Kicad board"""

        # Draw the coils
        for layer_id in self.config.copper_layers:
            if self.config.draw_only_layers is None or layer_id in self.config.draw_only_layers:
                for object in self.coils[layer_id]:
                    object.draw_kicad(
                        kicadpcb = kicadpcb,
                        layer = layer_id,
                    )

        # Draw the link traces
        for link in self.links:
            if self.config.draw_only_layers is None or link.layer in self.config.draw_only_layers:
                link.draw_kicad(kicadpcb)

        # Draw the vias
        if self.config.draw_vias:
            for via in self.vias:
                via.draw_kicad(kicadpcb)

        # Draw the terminals
        if self.config.draw_terminals:
            for terminal in self.terminals:
                terminal.draw_kicad(kicadpcb)
        
        # Draw the top and bottom silkscreens
        for item in self.top_silk:
            item.draw_kicad(kicadpcb)
        for item in self.bottom_silk:
            item.draw_kicad(kicadpcb)

        # Draw the outline
        if self.config.draw_outline:
            self.outline.draw_kicad(
                kicadpcb = kicadpcb,
                width = self.config.trace_width,
                marking_width = self.config.silk_line_width,
            )

        # Draw the construction geometry
        if self.config.draw_construction_geometry:
            for object in self.construction_geometry:
                object.draw_kicad(
                    kicadpcb = kicadpcb,
                    width = self.config.construction_geometry_line_width,
                    layer = 'construction',
                    stroke_type = 'dash_dot',
                )

        # Draw the magnets
        if self.config.draw_magnets:
            for i in range(self.config.n_magnets):
                center = Point.polar(360.0 * i / self.config.n_magnets, self.config.magnets_position_radius)
                kicadpcb.gr_circle(
                    center = center,
                    radius = self.config.magnets_diameter / 2.0,
                    width = self.config.magnets_line_width,
                    layer = 'magnets',
                )
