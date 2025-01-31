"""2D geometry classes

This module contains utility classes representing a few 2D geometric objects,
such as Vector, Point, Line, Arc, Circle and Path. These classes provide
some useful methods to perform computations on these objects, such as
finding intersection points, creating offset geometry and projections, ...

Note that this module is developped as part of the PCBMot project and
mostly provide functions required for the main project. While it defines
reusable components, this module is not anywhere close to be a generic,
fully-featured geometric solver.

Important note : all polar coordinates use the vertical axis as the origin
and clockwise as the positive direction.
"""

from typing import Self
import svgwrite as svg
import math

## Round the coordinates to this precision to account for float precision issues
PRECISION = 10


## SVG default drawing style
class SVGStyle:
    point_color: str = "#C02020"
    point_opacity: float = 0.5
    point_radius: float = 0.05
    point_thickness: float = 0.05

    line_color: str = "#37C837"
    line_opacity: float = 0.5
    line_thickness: float = 0.05
    line_dashes: str = "0.4 0.15"

    circle_color: str = "#37C837"
    circle_opacity: float = 0.5
    circle_thickness: float = 0.05
    circle_dashes: str = "0.4 0.15"

    arc_color: str = "#37C837"
    arc_opacity: float = 0.5
    arc_thickness: float = 0.05
    arc_dashes = "0.4 0.15"

    path_color: str = "#0044AA"
    path_opacity: float = 0.5
    path_thickness: float = 1.0

style = SVGStyle()


## Trigonometry functions using degrees

def sin(x: float) -> float:
    return math.sin(math.radians(x))

def cos(x: float) -> float:
    return math.cos(math.radians(x))

def tan(x: float) -> float:
    return math.tan(math.radians(x))

def asin(x: float) -> float:
    return math.degrees(math.asin(x))

def acos(x: float) -> float:
    return math.degrees(math.acos(x))

def atan(x: float) -> float:
    return math.degrees(math.atan(x))

def atan2(y: float, x: float) -> float:
    return math.degrees(math.atan2(y, x))


## Geometry classes

class Vector:
    """A 2D vector defined by its X and Y components"""

    def __init__(self, x: float, y: float):
        self.x: float = round(float(x), PRECISION)
        self.y: float = round(float(y), PRECISION)

    def __str__(self) -> str:
        return f"Vector(x={self.x}, y={self.y})"
    
    def __add__(self, other: Self) -> Self:
        match other:
            case Vector():
                return Vector(self.x + other.x, self.y + other.y)
            case _:
                raise TypeError(f"Unsupported object added to a Vector : {type(other)}")
    
    def __sub__(self, other) -> Self:
        return self + other * -1
    
    def __mul__(self, factor: float) -> Self:
        return Vector(self.x * factor, self.y * factor)
    
    def __eq__(self, other: Self) -> bool:
        return math.isclose(self.x, other.x) and math.isclose(self.y, other.y)
    
    def from_two_points(p1: "Point", p2: "Point") -> Self:
        """Create a new Vector based on p1 going to p2"""
        return Vector(p2.x - p1.x, p2.y - p1.y)

    def polar(a: float, r: float) -> Self:
        """Create a new Vector using polar coordinates"""
        x = r * sin(a)
        y = r * cos(a)
        return Vector(x, y)
    
    def to_polar(self) -> tuple[float, float]:
        """Return the polar coordinates of this Vector as a tuple (angle, length)"""
        return (self.angle(), self.length())
    
    def angle(self) -> float:
        """Return the angle of the polar coordinates of this Vector"""
        return 90 - atan2(self.y, self.x)
    
    def length(self) -> float:
        """Return the angle of the polar coordinates of this Vector"""
        return math.hypot(self.x, self.y)
    
    def from_tuple(p: tuple[float, float]) -> Self:
        return Vector(p[0], p[1])
    
    def as_tuple(self) -> tuple[float, float]:
        return (self.x, self.y)
    
    def dot(self, other: Self) -> float:
        """Calculate the dot-product of this vector with the given vector"""
        return self.x * other.x + self.y * other.y
    
    def cross(self, other: Self) -> float:
        """Calculate the cross-product of this vector with the given vector"""
        return self.x * other.y - other.x * self.y


class Point(Vector):
    """A point on the 2D plane, based on the Vector class"""

    def __str__(self) -> str:
        return f"Point(x={self.x}, y={self.y})"
    
    def __add__(self, other) -> Self:
        match other:
            case Vector() | Point():
                return Point(self.x + other.x, self.y + other.y)
            case _:
                raise TypeError(f"Unsupported object added to a Point : {type(other)}")
    
    def origin() -> Self:
        return Point(0.0, 0.0)

    def from_vector(vec: Vector) -> Self:
        return Point(vec.x, vec.y)
    
    def to_vector(self) -> Vector:
        return Vector(self.x, self.y)

    def polar(a: float, r: float) -> Self:
        return Point.from_vector(Vector.polar(a, r))
    
    def from_tuple(p: tuple[float, float]) -> Self:
        return Point.from_vector(Vector.from_tuple(p))
    
    def distance(self, object) -> float:
        """Calculate the distance between this Point and the given object

        Supported objects : Point, Line, Circle, Arc (extrapolated as a circle).
        Throws a TypeError if any other object type is given.
        """
        match object:
            case Point():
                return (self + object * -1).length()

            case Line():
                if math.isinf(object.a):
                    return math.fabs(object.b - self.y)
                else:
                    return math.fabs(object.a * self.y - self.x + object.b) / math.sqrt(1 + object.a**2)
            
            case Circle():
                return math.fabs(self.distance(object.center) - object.radius)
            
            case Arc():
                return math.fabs(self.distance(object.center()) - object.radius)

            case _:
                raise TypeError(f"Trying to compute the distance between a Point and an unsupported object : {type(object)}")
    
    def closest(self, objects: list) -> Self:
        """Return the object in the given list that is closest to this Point.

        See distance() for the list of supported objects.
        The given objects do not need to be of the same type as long as they are
        each of one of the supported types.
        Throws a TypeError if any other object type is given.
        """
        min_distance = 0.0
        closest = None
        for other in objects:
            distance = self.distance(other)
            if closest is None or distance < min_distance:
                min_distance = distance
                closest = other
        return closest
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Point rotated around the given center point by the given angle"""
        x = center.x + (self.x - center.x) * cos(-angle) - (self.y - center.y) * sin(-angle)
        y = center.y + (self.x - center.x) * sin(-angle) + (self.y - center.y) * cos(-angle)
        return Point(x, y)
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Point mirrored about the X axis"""
        return Point(self.x, -self.y)
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Point mirrored about the Y axis"""
        return Point(-self.x, self.y)

    def projected(self, object) -> Self:
        """Create a copy of this Point projected on the given object

        Supported objects : Line, Circle, Arc (extrapolated as a circle).
        Throws a TypeError if any other object type is given.
        """
        match object:
            case Line():
                return object.intersect(object.perpendicular(self))
            
            case Circle():
                projection_line = Line.from_two_points(object.center, self)
                return projection_line.intersect(object)
            
            case Arc():
                projection_line = Line.from_two_points(object.center(), self)
                return projection_line.intersect(object)

            case _:
                raise TypeError(f"Trying to project a Point on an unsupported object : {type(object)}")

    def to_viewport(self) -> Self:
        return Vector(round(self.x, PRECISION), round(-self.y, PRECISION))
    
    def to_svg(self) -> str:
        p = self.to_viewport()
        return f"{round(p.x, PRECISION)},{round(p.y, PRECISION)}"
    
    def draw_svg(self, drawing: svg.Drawing, radius=None, color=None, opacity=None, thickness=None):
        """Draw this Point on the given SVG drawing
        
        This method returns self and can therefore be chained."""
        drawing.add(drawing.circle(self.to_viewport().as_tuple(),
            radius or style.point_radius,
            stroke = color or style.point_color,
            stroke_opacity = opacity or style.point_opacity,
            stroke_width = thickness or style.point_thickness,
        ))
        return self


class Line:
    """A line on the 2D plane represented by the following equation `x = a * y + b`

    Note that this is different to the more common `y = a * x + b`.
    This unusual representation is used make vertical lines easier to represent
    and prevent loss of float precision in lines close to vertical, which are more
    common in PCBMot than horizontal lines.
    Horizontal lines are represented with `a=math.inf` and `y=b`.
    """

    def __init__(self, a: float, b: float):
        self.a: float = a
        self.b: float = b

    def __str__(self) -> str:
        return f"Line(a={self.a}, b={self.b})"
    
    def from_two_points(p1: Point, p2: Point) -> Self:
        """Create a new Line that connects the two given points

        If the given points are coincident, return None.
        """
        if math.isclose(p1.y, p2.y):
            if math.isclose(p1.x, p2.x):
                print("Warning: trying to create a Line using two coincident points")
                return None
            else:
                a = math.inf
                b = p1.y
        else:
            a = (p1.x - p2.x) / (p1.y - p2.y)
            b = p1.x - a * p1.y
        return Line(a, b)
    
    def angle(self) -> float:
        """Return the angle of this Line in the range [-90°:90°] relative to the vertical axis"""
        if math.isinf(self.a):
            return 90.0
        else:
            return atan(self.a)
    
    def is_parallel(self, other: Self) -> bool:
        """Check if this Line is parallel with the given line"""
        return math.isclose(self.a, other.a)
    
    def is_colinear(self, other: Self) -> bool:
        """Check if this Line is colinear with the given line"""
        return math.isclose(self.a, other.a) and math.isclose(self.b, other.b)

    def intersect(self, object) -> Point:
        """Calculate the intersection Point(s) between this Line and the given object

        Supported objects : Line, Circle, Arc.
        Throws a TypeError if any other object type is given.
        Returns either a single Point, a tuple of two Points, or None if there is no
        intersection between the objects.
        """
        match object:
            case Line():
                if math.isclose(self.a, object.a):
                    print("Warning: trying to compute intersection between two parallel or colinear lines")
                    return None
                elif math.isinf(object.a):
                    return object.intersect(self)
                elif math.isinf(self.a):
                    x = object.a * self.b + object.b
                    y = self.b
                    return Point(x, y)
                else:
                    if math.isclose(self.a, 0.0):
                        x = self.b
                        y = (x - object.b) / object.a
                    else:
                        x = (self.a * object.b - object.a * self.b) / (self.a - object.a)
                        y = (x - self.b) / self.a
                    return Point(x, y)

            case Circle():
                if math.isinf(self.a):
                    A = 1
                    B = -2 * object.center.x
                    C = (self.b - object.center.y)**2 + object.center.x**2 - object.radius**2
                else:
                    A = 1 + self.a**2
                    B = 2 * (self.a * (self.b - object.center.x) - object.center.y)
                    C = (self.b - object.center.x)**2 + object.center.y**2 - object.radius**2
                discriminant = B**2 - 4 * A * C
                if discriminant < 0:
                    print("Warning: trying to compute intersection between a non-intersecting circle and line")
                    return None
                sqrt_discriminant = math.sqrt(discriminant)
                if math.isinf(self.a):
                    x1 = (-B + sqrt_discriminant) / (2 * A)
                    x2 = (-B - sqrt_discriminant) / (2 * A)
                    y = self.b
                    return (Point(x1, y), Point(x2, y))
                else:
                    y1 = (-B + sqrt_discriminant) / (2 * A)
                    y2 = (-B - sqrt_discriminant) / (2 * A)
                    x1 = self.a * y1 + self.b
                    x2 = self.a * y2 + self.b
                    return (Point(x1, y1), Point(x2, y2))

            case Arc():
                circle = Circle(object.center(), object.radius)
                intersect = self.intersect(circle)
                if intersect is None:
                    return None
                p1, p2 = intersect
                distance_p1 = min(p1.distance(object.p1), p1.distance(object.p2))
                distance_p2 = min(p2.distance(object.p1), p2.distance(object.p2))
                if distance_p1 < distance_p2:
                    return p1
                else:
                    return p2

            case _:
                raise TypeError(f"Unsupported object given to intersect a Line : {type(object)}")
    
    def offset(self, distance: float) -> Self:
        """Create a new Line parallel to this Line separated by the given distance

        The given distance can be negative to flip the side of the offset.
        """
        if math.isinf(self.a):
            return Line(self.a, self.b + distance)
        else:
            angle = self.angle()
            db = distance / cos(angle)
            return Line(self.a, self.b + db)

    def perpendicular(self, point: Point) -> Self:
        """Create a new Line perpendicular to this Line that passes through the given Point"""
        if math.isinf(self.a):
            return Line(0, point.x)
        elif math.isclose(self.a, 0):
            return Line(math.inf, point.y)
        else:
            a = - 1 / self.a
            b = point.x + point.y / self.a
            return Line(a, b)

    def eval(self, y: float) -> float:
        if math.isinf(self.a):
            print("Warning: calling eval() on an horizontal Line")
        return self.a * y + self.b
    
    def draw_svg(self, drawing: svg.Drawing, color=None, opacity=None, thickness=None, dashes=None):
        """Draw this Line on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if math.isinf(self.a):
            x1 = -VIEWPORT_WIDTH/2.0
            x2 = VIEWPORT_WIDTH/2.0
            y = self.b
            p1 = (round(x1, PRECISION), round(-y, PRECISION))
            p2 = (round(x2, PRECISION), round(-y, PRECISION))
        else:
            y1 = -VIEWPORT_HEIGHT/2.0
            y2 = VIEWPORT_HEIGHT/2.0
            p1 = (round(self.eval(y1), PRECISION), round(-y1, PRECISION))
            p2 = (round(self.eval(y2), PRECISION), round(-y2, PRECISION))
        
        drawing.add(drawing.line(
            p1, p2,
            stroke = color or style.line_color,
            stroke_opacity = opacity or style.line_opacity,
            stroke_width = thickness or style.line_thickness,
            stroke_dasharray = dashes or style.line_dashes,
        ))
        return self


class Circle:
    """A 2D Circle represented by its center and radius"""

    def __init__(self, center: Point, radius: float):
        self.center: Point = center
        self.radius: float = radius

    def __str__(self) -> str:
        return f"Circle(center={self.center}, radius={self.radius})"

    def intersect(self, object) -> Point:
        """Calculate the intersection Point between this Circle and the given object

        Supported objects : Line, Circle.
        Throws a TypeError if any other object type is given.
        Returns either a single Point, a tuple of two Points, or None if there is no
        intersection between the objects.
        """
        match object:
            case Line():
                return object.intersect(self)
            
            case Circle():
                if self.center == object.center and math.isclose(self.radius, object.radius):
                    # Same circles
                    print("Warning: trying to compute intersection between two identical circles")
                    return None
                distance = math.hypot(object.center.x - self.center.x, object.center.y - self.center.y)
                if distance > self.radius + object.radius or distance < math.fabs(self.radius - object.radius):
                    # No intersection
                    print("Warning: trying to compute intersection between non-intersecting circles")
                    return None
                a = (self.radius**2 - object.radius**2 + distance**2) / (2.0 * distance)
                h = math.sqrt(self.radius**2 - a**2)
                x0 = self.center.x + a * (object.center.x - self.center.x) / distance
                y0 = self.center.y + a * (object.center.y - self.center.y) / distance
                x1 = x0 + h * (object.center.y - self.center.y) / distance
                y1 = y0 - h * (object.center.x - self.center.x) / distance
                x2 = x0 - h * (object.center.y - self.center.y) / distance
                y2 = y0 + h * (object.center.x - self.center.x) / distance
                p1 = Point(x1, y1)
                p2 = Point(x2, y2)
                if p1 == p2:
                    return p1
                else:
                    return (p1, p2)

            case _:
                raise TypeError(f"Unsupported object given to intersect a Circle : {type(object)}")
    
    def offset(self, distance: float) -> Self:
        """Create a new Circle concentric with this Circle offset by the given distance

        The given distance can be negative to flip the side of the offset.
        """
        radius = self.radius + distance
        if radius <= 0:
            return None
        return Circle(self.center, radius)
    
    def draw_svg(self, drawing: svg.Drawing, color=None, opacity=None, thickness=None, dashes=None):
        """Draw this Circle on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        drawing.add(drawing.circle(
            self.center.to_viewport().as_tuple(),
            self.radius,
            stroke = color or style.circle_color,
            stroke_opacity = opacity or style.circle_opacity,
            stroke_width = thickness or style.circle_thickness,
            stroke_dasharray = dashes or style.circle_dashes,
        ))
        return self


class Arc:
    """A 2D Arc represented by its two endpoints and its radius

    The arc always go clockwise from p1 to p2.
    """

    def __init__(self, p1: Point, p2: Point, radius: float, reverse: bool = False):
        """Create a new Arc using the to given endpoints and the given radius

        If reverse is set, p1 and p2 and simply swapped before constructing the Arc.
        """
        if p1 == p2:
            raise ValueError("Error : unable to create an arc using two coincident points")
        distance = p1.distance(p2)
        min_radius = distance / 2.0
        if radius < min_radius:
            print(f"Warning : trying to create an Arc with a radius (r={radius}) too small for the distance between the two points (d={distance}), setting radius to {min_radius}")
            radius = min_radius
        if reverse:
            p1, p2 = p2, p1
        self.p1: Point = p1
        self.p2: Point = p2
        self.radius: Point = radius

    def __str__(self) -> str:
        return f"Arc(p1={self.p1}, p2={self.p2}, radius={self.radius})"
    
    def __mul__(self, factor: float) -> Self:
        """Create a new Arc concentric to this Arc with the same p1, and with
        p2 rotated around the center so that the length of the arc is grown or
        shrinked by the given factor.
        """
        return Arc(self.p1, self.p1.rotated(self.center(), self.angle() * factor), self.radius)
    
    def center(self) -> Point:
        """Calculate the center Point of this Arc"""
        circle1 = Circle(self.p1, self.radius)
        circle2 = Circle(self.p2, self.radius)
        c0, c1 = circle1.intersect(circle2)
        v1 = Vector.from_two_points(c0, self.p1)
        v2 = Vector.from_two_points(c0, self.p2)
        if v1.cross(v2) > 0:
            return c1
        else:
            return c0

    def angle(self) -> float:
        """Calculate the angle of this Arc"""
        center = self.center()
        v1 = Vector.from_two_points(center, self.p1)
        v2 = Vector.from_two_points(center, self.p2)
        return acos(v1.dot(v2) / (v1.length() * v2.length()))

    def length(self) -> float:
        """Calculate the curve length of this Arc"""
        return math.pi * self.radius * self.angle() / 180.

    def intersect(self, object) -> Point:
        """Calculate the intersection Point between this Arc and the given object

        Supported objects : Line.
        Throws a TypeError if any other object type is given.
        Returns either a single Point, a tuple of two Points, or None if there is no
        intersection between the objects.
        """
        match object:
            case Line():
                return object.intersect(self)
            case _:
                raise TypeError(f"Unsupported object given to intersect an Arc : {type(object)}")
    
    def radial_line(self, point: Point) -> Line:
        """Create a new Line radial to this Arc that passes through the given Point"""
        return Line.from_two_points(self.center(), point)

    def offset(self, distance: float) -> Self:
        """Create a new Arc concentric with this Arc offset by the given distance

        The given distance can be negative to flip the side of the offset.
        """
        radius_offset = max(self.radius + distance, 0.0)
        center = self.center()
        line1 = Line.from_two_points(center, self.p1)
        line2 = Line.from_two_points(center, self.p2)
        circle = Circle(center, radius_offset)
        p1_offset = self.p1.closest(line1.intersect(circle))
        p2_offset = self.p2.closest(line2.intersect(circle))
        return Arc(p1_offset, p2_offset, radius_offset)
    
    def draw_svg(self, drawing: svg.Drawing, color=None, opacity=None, thickness=None, dashes=None):
        """Draw this Arc on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        drawing.add(drawing.path(
            d=f"M{self.p1.to_svg()} A{round(self.radius, PRECISION)},{round(self.radius, PRECISION)} 0 0,1 {self.p2.to_svg()}",
            stroke = color or style.arc_color,
            stroke_opacity = opacity or style.arc_opacity,
            stroke_width = thickness or style.arc_thickness,
            stroke_dasharray = dashes or style.arc_dashes,
        ))
        return self


class Fillet:
    """A fillet between two objects
    
    Supported types : two segments, or a segment and an arc."""

    def __init__(self, p_arc_prev: Point, p_arc_next: Point, fillet_radius: float, anticlockwise_arc: bool, shrinked: bool):
        self.p_arc_prev: Point = p_arc_prev
        self.p_arc_next: Point = p_arc_next
        self.fillet_radius: float = fillet_radius
        self.anticlockwise_arc: bool = anticlockwise_arc
        self.shrinked: bool = shrinked

    def segment_to_segment(p_prev: Point, p_mid: Point, p_next: Point, fillet_radius: float) -> Self:
        """Compute the Fillet between two connected segments

        The segments are represented respectively by (p_prev, p_mid) and (p_mid, p_next),
        and therefore share the p_mid point. The fillet is applied on this p_mid point.
        If the fillet is too large for the the given geometry, its radius will be
        automatically shrinked to allow it to fit the geometry. Check the `shrinked`
        attribute on the returned Fillet object if necessary.
        
        Return None if the two segments are colinear.
        """

        # Check inputs
        if fillet_radius <= 0:
            raise ValueError("Creating a fillet requires a strictly positive radius")
        
        # Compute the two lines adjacent to the mid point
        line_prev = Line.from_two_points(p_prev, p_mid)
        line_next = Line.from_two_points(p_mid, p_next)

        # If the two lines are colinear, ignore the fillet
        if line_prev.is_colinear(line_next):
            return None
        
        shrinked = False
        for i in range(3):
            # Compute the center of the fillet by offseting each of the two lines on the direction of the opposite point,
            # then taking the intersection between them
            line_prev_offset = p_next.closest([line_prev.offset(fillet_radius), line_prev.offset(-fillet_radius)])
            line_next_offset = p_prev.closest([line_next.offset(fillet_radius), line_next.offset(-fillet_radius)])
            fillet_center = line_prev_offset.intersect(line_next_offset)

            # Project the center point of the fillet on the adjacent lines to compute the endpoints of the arc of the fillet
            p_fillet_prev = fillet_center.projected(line_prev)
            p_fillet_next = fillet_center.projected(line_next)

            # Check if these endpoints are inside the adjacent segments
            v_prev = Vector.from_two_points(p_mid, p_prev)
            v_fillet_prev = Vector.from_two_points(p_mid, p_fillet_prev)
            v_next = Vector.from_two_points(p_mid, p_next)
            v_fillet_next = Vector.from_two_points(p_mid, p_fillet_next)
            cut = False
            if v_fillet_prev.length() >= v_prev.length():
                # Compute a new pair of fillet endpoints that are inside the prev segment
                print("Warning : fillet too big for the previous line, cutting the fillet")
                length_1 = v_fillet_prev.length()
                p_fillet_prev = p_mid + v_prev * 0.9
                v_fillet_prev = Vector.from_two_points(p_mid, p_fillet_prev)
                length_2 = v_fillet_prev.length()
                factor = length_2 / length_1
                v_fillet_next = v_fillet_next * factor
                p_fillet_next = p_mid + v_fillet_next
                cut = True
                shrinked = True
            if v_fillet_next.length() >= v_next.length():
                # Compute a new pair of fillet endpoints that are inside the next segment
                print("Warning : fillet too big for the next line, cutting the fillet")
                length_1 = v_fillet_next.length()
                p_fillet_next = p_mid + v_next * 0.9
                v_fillet_next = Vector.from_two_points(p_mid, p_fillet_next)
                length_2 = v_fillet_next.length()
                factor = length_2 / length_1
                v_fillet_prev = v_fillet_prev * factor
                p_fillet_prev = p_mid + v_fillet_prev
                cut = True
                shrinked = True
            if cut:
                # The fillet was modified : compute its new radius and start again
                l1 = line_prev.perpendicular(p_fillet_prev)
                l2 = line_next.perpendicular(p_fillet_next)
                new_center = l1.intersect(l2)
                fillet_radius = new_center.distance(line_prev)
                continue
            
            # p_fillet_prev.draw_svg(dwg)
            # p_fillet_next.draw_svg(dwg)
            # Arc(p_fillet_prev, p_fillet_next, fillet_radius).draw_svg(dwg)
            
            # The points are kept in the order of the path, so we must return a separate flag
            # to indicate an anticlockwise arc direction
            anticlockwise_fillet_arc = v_prev.cross(v_next) < 0

            # Return the two endpoints of the fillet, its finale radius, and the anticlockwise flag
            return Fillet(p_fillet_prev, p_fillet_next, fillet_radius, anticlockwise_fillet_arc, shrinked)

    def segment_to_arc(p_prev: Point, p_mid: Point, p_next: Point, arc_radius: float, anticlockwise: bool, fillet_radius: float) -> Self:
        """Compute the Fillet between a segment connected to an arc of the given radius

        The segment and the arc are represented respectively by (p_prev, p_mid) and (p_mid, p_next),
        and therefore share the p_mid point. The fillet is applied on this p_mid point.
        Set the `anticlockwise` flag if the arc goes anti-clockwise from p_mid to p_next.
        If the fillet is too large for the the given geometry, its radius will be
        automatically shrinked to allow it to fit the geometry. Check the `shrinked`
        attribute on the returned Fillet object if necessary.
        """

        # Compute the line and arc adjacent to the mid point
        line_prev = Line.from_two_points(p_prev, p_mid)
        arc_next = Arc(p_mid, p_next, arc_radius, reverse=anticlockwise)
        
        shrinked = False
        for i in range(3):
            # Compute the center of the fillet by offseting each of the two lines on the direction of the opposite point,
            # then taking the intersection between them
            line_prev_offset = p_next.closest([line_prev.offset(fillet_radius), line_prev.offset(-fillet_radius)])
            arc_next_offset = p_prev.closest([arc_next.offset(fillet_radius), arc_next.offset(-fillet_radius)])
            fillet_center = line_prev_offset.intersect(arc_next_offset)
            if fillet_center is None:
                # TODO : automatically try with a smaller fillet radius
                print("Warning : fillet too big for the given geometry")
                return None

            # Project the center point of the fillet on the adjacent lines to compute the endpoints of the arc of the fillet
            p_fillet_prev = fillet_center.projected(line_prev)
            p_fillet_next = fillet_center.projected(arc_next)
            # Arc(p_fillet_next, p_fillet_prev, fillet_radius, reverse=anticlockwise).draw_svg(dwg)

            # Check if these endpoints are inside the adjacent segments
            v_prev = Vector.from_two_points(p_mid, p_prev)
            v_fillet_prev = Vector.from_two_points(p_mid, p_fillet_prev)
            p_fillet_next
            arc_fillet_next = Arc(p_mid, p_fillet_next, arc_next.radius, reverse=anticlockwise)
            cut = False
            # TODO : replace the 'factor' calculation with something more robust
            if v_fillet_prev.length() >= v_prev.length():
                # Compute a new pair of fillet endpoints that are inside the prev segment
                print("Warning : fillet too big for the previous line, cutting the fillet")
                length_1 = v_fillet_prev.length()
                p_fillet_prev = (p_mid + v_prev * 0.9)
                v_fillet_prev = Vector.from_two_points(p_mid, p_fillet_prev)
                length_2 = v_fillet_prev.length()
                factor = length_2 / length_1
                arc_fillet_next = arc_fillet_next * factor
                p_fillet_next = arc_fillet_next.p2
                cut = True
                shrinked = True
            if arc_fillet_next.length() >= arc_next.length():
                # Compute a new pair of fillet endpoints that are inside the next arc
                print("Warning : fillet too big for the next line, cutting the fillet")
                length_1 = arc_fillet_next.length()
                arc_fillet_next = arc_fillet_next * 0.7
                p_fillet_next = arc_fillet_next.p2
                length_2 = arc_fillet_next.length()
                factor = length_2 / length_1
                v_fillet_prev = v_fillet_prev * factor
                p_fillet_prev = p_mid + v_fillet_prev
                cut = True
                shrinked = True
            if cut:
                # The fillet was modified : compute its new radius and start again
                l1 = line_prev.perpendicular(p_fillet_prev)
                l2 = arc_next.radial_line(p_fillet_next)
                new_center = l1.intersect(l2)
                fillet_radius = min(new_center.distance(line_prev), new_center.distance(arc_next))
                Arc(p_fillet_prev, p_fillet_next, fillet_radius, reverse=anticlockwise)
                continue
            
            # p_fillet_prev.draw_svg(dwg)
            # p_fillet_next.draw_svg(dwg)
            # Arc(p_fillet_prev, p_fillet_next, fillet_radius).draw_svg(dwg)

            # The points are kept in the order of the path, so we must return a separate flag
            # to indicate an anticlockwise arc direction
            anticlockwise_arc = Vector.from_two_points(p_mid, p_fillet_prev).cross(Vector.from_two_points(p_mid, p_fillet_next)) < 0

            # Return the two endpoints of the fillet, its finale radius, and the anticlockwise flag
            return Fillet(p_fillet_prev, p_fillet_next, fillet_radius, anticlockwise_arc, shrinked)

    def arc_to_segment(p_prev: Point, p_mid: Point, p_next: Point, arc_radius: float, anticlockwise: bool, fillet_radius: float) -> Self:
        """Compute the Fillet between an arc of the given radius connected to a segment

        The arc and the segment are represented respectively by (p_prev, p_mid) and (p_mid, p_next),
        and therefore share the p_mid point. The fillet is applied on this p_mid point.
        Set the `anticlockwise` flag if the arc goes anti-clockwise from p_prev to p_mid.
        If the fillet is too large for the the given geometry, its radius will be
        automatically shrinked to allow it to fit the geometry. Check the `shrinked`
        attribute on the returned Fillet object if necessary.
        """

        # Create the equivalent Fillet from the segment to the arc, and reverse it
        fillet = Fillet.segment_to_arc(p_next, p_mid, p_prev, arc_radius, not anticlockwise, fillet_radius)
        if fillet is None:
            return None
        fillet.p_arc_prev, fillet.p_arc_next = fillet.p_arc_next, fillet.p_arc_prev
        fillet.anticlockwise_arc = not fillet.anticlockwise_arc
        return fillet


class PathElement:
    """Empty base class for specific object types inside a Path"""
    pass

class PathSegment(PathElement):
    """A segment as part of a Path"""

    def __init__(self, end_point: Point):
        self.p2: Point = end_point

    def __str__(self) -> str:
        return f"PathSegment(p2={self.p2})"
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this PathSegment rotated around the given center point by the given angle"""
        return PathSegment(self.p2.rotated(center, angle))
    
    def mirrored_x(self) -> Self:
        """Create a copy of this PathSegment mirrored about the X axis"""
        return PathSegment(self.p2.mirrored_x())
    
    def mirrored_y(self) -> Self:
        """Create a copy of this PathSegment mirrored about the Y axis"""
        return PathSegment(self.p2.mirrored_y())

class PathArc(PathElement):
    """An arc segment as part of a Path"""

    def __init__(self, end_point: Point, radius: float, anticlockwise: bool):
        self.p2: Point = end_point
        self.radius: Point = radius
        self.anticlockwise: Point = anticlockwise

    def __str__(self) -> str:
        return f"PathArc(p2={self.p2}, radius={self.radius}, anticlockwise={self.anticlockwise})"
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this PathArc rotated around the given center point by the given angle"""
        return PathArc(self.p2.rotated(center, angle), self.radius, self.anticlockwise)
    
    def mirrored_x(self) -> Self:
        """Create a copy of this PathArc mirrored about the X axis"""
        return PathArc(self.p2.mirrored_x(), self.radius, not self.anticlockwise)
    
    def mirrored_y(self) -> Self:
        """Create a copy of this PathArc mirrored about the Y axis"""
        return PathArc(self.p2.mirrored_y(), self.radius, not self.anticlockwise)

class Path:
    """A Path consisting of a list of connected segments and arcs"""

    def __init__(self, start: Point):
        self.start: Point = start
        self.elements: list[PathElement] = []
        self.current_point: Point = start
    
    def previous_point(self, n: int) -> Point:
        if n < 0:
            return None
        if n == 0:
            return self.current_point
        elif n < len(self.elements):
            return self.elements[-(n+1)].p2
        elif n == len(self.elements):
            return self.start
        else:
            return None
    
    def append_segment(self, to_point: Point, fillet_radius: float = None):
        """Append a new segment at the end of this Path, with an optional fillet"""

        if to_point == self.current_point:
            return
        
        # Fillet between this segment and the previous element in the path
        if not self.elements:
            fillet_radius = None
        if fillet_radius:
            # Get the three points adjacent to the current point
            p_prev = self.previous_point(1)
            p_mid = self.current_point
            p_next = to_point
            last_element = self.elements[-1]
            match last_element:
                case PathSegment():
                    # Compute the fillet
                    fillet = Fillet.segment_to_segment(p_prev, p_mid, p_next, fillet_radius)

                    if fillet is not None:
                        # Replace the last segment with the fillet
                        self.elements.pop()
                        self.elements.append(PathSegment(fillet.p_arc_prev))
                        self.elements.append(PathArc(fillet.p_arc_next, fillet.fillet_radius, fillet.anticlockwise_arc))
                        self.current_point = fillet.p_arc_next
                
                case PathArc():
                    # Compute the fillet
                    fillet = Fillet.arc_to_segment(p_prev, p_mid, p_next, last_element.radius, last_element.anticlockwise, fillet_radius)

                    if fillet is not None:
                        # Replace the last segment with the fillet
                        self.elements.pop()
                        self.elements.append(PathArc(fillet.p_arc_prev, last_element.radius, last_element.anticlockwise))
                        self.elements.append(PathArc(fillet.p_arc_next, fillet.fillet_radius, fillet.anticlockwise_arc))
                        self.current_point = fillet.p_arc_next

        self.elements.append(PathSegment(to_point))
        self.current_point = to_point
    
    def append_arc(self, to_point: Point, radius: float, anticlockwise: bool, fillet_radius: float = None):
        """Append a new arc of the given radius at the end of this Path, with an optional fillet"""

        if to_point == self.current_point:
            return

        # Fillet between this arc and the previous element in the path
        if not self.elements:
            fillet_radius = None
        if fillet_radius:
            # Get the three points adjacent to the current point
            p_prev = self.previous_point(1)
            p_mid = self.current_point
            p_next = to_point
            last_element = self.elements[-1]
            match last_element:
                case PathSegment():
                    # Compute the fillet
                    fillet = Fillet.segment_to_arc(p_prev, p_mid, p_next, radius, anticlockwise, fillet_radius)

                    if fillet is not None:
                        # Replace the last segment with the fillet
                        self.elements.pop()
                        self.elements.append(PathSegment(fillet.p_arc_prev))
                        self.elements.append(PathArc(fillet.p_arc_next, fillet.fillet_radius, fillet.anticlockwise_arc))
                        self.current_point = fillet.p_arc_next
                
                case PathArc():
                    raise ValueError("Fillet between two arcs is not supported")

        self.elements.append(PathArc(to_point, radius, anticlockwise))
        self.current_point = to_point
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Path rotated around the given center point by the given angle"""
        start = self.start.rotated(center, angle)
        path = Path(start)
        for element in self.elements:
            path.elements.append(element.rotated(center, angle))
        path.current_point = self.current_point.rotated(center, angle)
        return path
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Path mirrored about the X axis"""
        start = self.start.mirrored_x()
        path = Path(start)
        for element in self.elements:
            path.elements.append(element.mirrored_x())
        path.current_point = self.current_point.mirrored_x()
        return path
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Path mirrored about the Y axis"""
        start = self.start.mirrored_y()
        path = Path(start)
        for element in self.elements:
            path.elements.append(element.mirrored_y())
        path.current_point = self.current_point.mirrored_y()
        return path

    def to_svg(self) -> str:
        d = ""
        d += f"M{self.start.to_svg()} "
        current_point = self.start
        for element in self.elements:
            match element:
                case PathSegment():
                    d += f"L{element.p2.to_svg()} "
                    current_point = element.p2
                case PathArc():
                    d += f"A{round(element.radius, PRECISION)},{round(element.radius, PRECISION)} 0 0,{int(not element.anticlockwise)} {element.p2.to_svg()} "
                    current_point = element.p2
        return d
    
    def draw_svg(self, drawing: svg.Drawing, color=None, opacity=None, thickness=None, dashes=None):
        """Draw this Path on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        drawing.add(drawing.path(
            d = self.to_svg(),
            stroke = color or style.path_color,
            stroke_opacity = opacity or style.path_opacity,
            stroke_width = thickness or style.path_thickness,
            stroke_dasharray = dashes or style.circle_dashes,
        ))
        return self