from typing import Self
from config import Config
import svgwrite as svg
from geometry import sin, cos, tan, asin, acos, atan, atan2, Vector, Point, Line, Circle, Arc, Path

class Coil:
    """A single coil on the board"""

    def __init__(self, path: Path):
        self.path = path

    def generate(angle: float, outer_radius: float, inner_radius: float, loop_offset: float) -> Self:
        """Generate a clockwise coil centered around the vertical axis based on the given parameters"""

        outer_radius_initial = outer_radius
        inner_radius_initial = inner_radius
        outer_fillet_radius = 1.2
        inner_fillet_radius = 2.2
        point_start = Point(0, outer_radius)
        path = Path(point_start)
        for i in range(13):
            # Construction geometry
            line_left = Line.from_two_points(Point.origin(), Point.polar(-angle/2.0, outer_radius)).offset(loop_offset * (i + 0.5))
            line_right = Line.from_two_points(Point.origin(), Point.polar(angle/2.0, outer_radius)).offset(-loop_offset * (i + 0.5))
            arc_outer = Arc(Point.polar(-angle/2.0, outer_radius), Point.polar(angle/2.0, outer_radius), outer_radius)
            inner_radius = inner_radius_initial + i * loop_offset
            arc_inner = Arc(Point.polar(-angle/2.0, inner_radius), Point.polar(angle/2.0, inner_radius), inner_radius)

            # Arc to outer right
            point_outer_right = line_right.intersect(arc_outer)
            path.append_arc(point_outer_right, outer_radius, anticlockwise=False, fillet_radius=outer_fillet_radius)

            # Segment to inner right
            point_inner_right = line_right.intersect(arc_inner)
            path.append_segment(point_inner_right, fillet_radius=outer_fillet_radius)

            # Arc to inner left
            point_inner_left = line_left.intersect(arc_inner)
            path.append_arc(point_inner_left, inner_radius, anticlockwise=True, fillet_radius=inner_fillet_radius)

            # Construction outer arc with offset
            outer_radius = outer_radius_initial - (i + 1) * loop_offset
            p1 = Point.polar(-angle/2.0, outer_radius)
            p2 = Point.polar(angle/2.0, outer_radius)
            arc_outer = Arc(p1, p2, outer_radius)

            # Segment to outer left with offset
            point_outer_left = line_left.intersect(arc_outer)
            path.append_segment(point_outer_left, fillet_radius=inner_fillet_radius)

            # Reduce the fillet radius for the next loop
            outer_fillet_radius -= loop_offset
            if outer_fillet_radius < 0.1:
                outer_fillet_radius = 0.1
            inner_fillet_radius -= loop_offset
            if inner_fillet_radius < 0.1:
                inner_fillet_radius = 0.1
        
        return Coil(path)
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Coil rotated around the given center point by the given angle"""
        path = self.path.rotated(center, angle)
        return Coil(path)
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Coil mirrored about the X axis"""
        path = self.path.mirrored_x()
        return Coil(path)
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Coil mirrored about the Y axis"""
        path = self.path.mirrored_y()
        return Coil(path)
    
    def draw_svg(self, drawing: svg.Drawing, color=None, opacity=None, thickness=None, dashes=None):
        """Draw this Path on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        self.path.draw_svg(drawing, color, opacity, thickness, dashes)
        return self

class PCB:
    """A PCB containing layers"""

    def __init__(self, config: Config, layers: dict):
        self.config = config
        self.layers = layers
    
    def generate(config: Config):
        """Generate a new PCB based on the given config"""

        # Board outline
        outline = [
            Circle(Point.origin(), config.board_radius),
            Circle(Point.origin(), config.hole_radius),
        ]

        # Generate the coils
        coil_B = Coil.generate(
            angle=config.coil_angle,
            outer_radius=config.board_radius - config.board_outer_margin,
            inner_radius=config.hole_radius + config.board_inner_margin,
            loop_offset=config.trace_spacing + config.trace_width
        )
        coil_A = coil_B.rotated(Point.origin(), -config.coil_angle)
        coil_C = coil_B.rotated(Point.origin(), config.coil_angle)
        coil_A2 = coil_C.mirrored_x()
        coil_B2 = coil_B.mirrored_x()
        coil_C2 = coil_A.mirrored_x()

        # Stack the layers and create the PCB
        layers = {
            'top': [coil_A, coil_B, coil_C, coil_A2, coil_B2, coil_C2],
            'outline': outline,
        }
        return PCB(config, layers)
    
    def draw_svg(self, drawing: svg.Drawing, construction_geometry: bool = False, only_layers: list[str] = None) -> Self:
        """Draw this PCB on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        # Draw the construction geometry
        if construction_geometry:
            drawing.add(drawing.line(
                Point.origin().as_tuple(), Point(0, self.config.board_radius).to_viewport().as_tuple(),
                stroke = self.config.construction_geometry_color,
                stroke_width = self.config.construction_geometry_thickness,
                stroke_dasharray=self.config.construction_geometry_dashes,
            ))
            drawing.add(drawing.line(
                Point.origin().as_tuple(), Point.polar(self.config.coil_angle/2.0, self.config.board_radius).to_viewport().as_tuple(),
                stroke = self.config.construction_geometry_color,
                stroke_width = self.config.construction_geometry_thickness,
                stroke_dasharray = self.config.construction_geometry_dashes,
            ))
            drawing.add(drawing.line(
                Point.origin().as_tuple(), Point.polar(-self.config.coil_angle/2.0, self.config.board_radius).to_viewport().as_tuple(),
                stroke = self.config.construction_geometry_color,
                stroke_width = self.config.construction_geometry_thickness,
                stroke_dasharray = self.config.construction_geometry_dashes,
            ))

        # Draw the layers
        for layer, objects in self.layers.items():
            if only_layers is None or layer in only_layers:
                for object in objects:
                    object.draw_svg(
                        drawing,
                        color = self.config.layers_color.get(layer),
                        thickness = self.config.trace_width,
                        dashes = "none"
                    )

        return self