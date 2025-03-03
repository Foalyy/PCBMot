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
from enum import Enum
import svgwrite as svg
import math

_d = None
def set_drawing(d):
    global _d
    _d = d

## Round the coordinates to this precision to account for float precision issues
PRECISION = 10


## SVG default drawing style
class SVGStyle:
    point_color: str = "#C02020"
    point_opacity: float = 0.7
    point_radius: float = 0.1
    point_line_width: float = 0.03

    line_color: str = "#37C837"
    line_opacity: float = 0.7
    line_width: float = 0.05
    line_dashes: str = "0.2 0.12"

    circle_color: str = "#37C837"
    circle_opacity: float = 0.7
    circle_line_width: float = 0.05
    circle_dashes: str = "0.2 0.12"

    arc_color: str = "#37C837"
    arc_opacity: float = 0.7
    arc_line_width: float = 0.05
    arc_dashes = "0.2 0.12"

    path_color: str = "#37C837"
    path_opacity: float = 0.7
    path_line_width: float = 0.05
    path_dashes = "0.2 0.12"

    polygon_color: str = "#379fc8"
    polygon_opacity: float = 0.4
    polygon_stroke_color: str = "#1519db"
    polygon_stroke_opacity: float = 0.7
    polygon_stroke_width: float = 0.05
    polygon_stroke_dashes = "0.2 0.12"

style = SVGStyle()


class CapStyle(Enum):
    """Style to apply on the end points of a Path when creating a stroke polygon"""
    FLAT = 1,
    SQUARE = 2,
    ROUND = 3,

class JoinStyle(Enum):
    """Style to apply on the corners of a Path when creating an offset or a stroke polygon"""
    MITRE = 1,
    BEVEL = 2,
    ROUND = 3,


## Trigonometry functions using degrees

def isclose(value1, value2, abs_tol=1e-9) -> bool:
    return math.isclose(value1, value2, abs_tol=abs_tol)

def sin(x: float) -> float:
    return math.sin(math.radians(x))

def cos(x: float) -> float:
    return math.cos(math.radians(x))

def tan(x: float) -> float:
    return math.tan(math.radians(x))

def asin(x: float) -> float:
    # Prevent out-of-domain errors due to float rounding errors
    if isclose(x, -1):
        return math.degrees(-math.pi / 2.0)
    elif isclose(x, 1):
        return math.degrees(math.pi / 2.0)
    else:
        return math.degrees(math.asin(x))

def acos(x: float) -> float:
    # Prevent out-of-domain errors due to float rounding errors
    if isclose(x, -1):
        return math.degrees(math.pi)
    elif isclose(x, 1):
        return 0
    else:
        return math.degrees(math.acos(x))

def atan(x: float) -> float:
    return math.degrees(math.atan(x))

def atan2(y: float, x: float) -> float:
    return math.degrees(math.atan2(y, x))


## Geometry classes

class DrawableObject:
    pass

class Vector:
    """A 2D vector defined by its X and Y components"""

    def __init__(self, x: float, y: float):
        self.x: float = round(float(x), PRECISION)
        self.y: float = round(float(y), PRECISION)

    def __str__(self) -> str:
        return f"Vector(x={self.x:.2f}, y={self.y:.2f})"
    
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
    
    def __truediv__(self, factor: float) -> Self:
        return Vector(self.x / factor, self.y / factor)
    
    def __eq__(self, other: Self) -> bool:
        return isclose(self.x, other.x) and isclose(self.y, other.y)
    
    def from_two_points(p1: "Point", p2: "Point") -> Self:
        """Create a new Vector based on p1 pointing to p2"""
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
        """Return the length (norm) of this Vector"""
        return math.hypot(self.x, self.y)
    
    def from_tuple(p: tuple[float, float]) -> Self:
        return Vector(p[0], p[1])
    
    def as_tuple(self) -> tuple[float, float]:
        return (self.x, self.y)
    
    def normalized(self) -> Self:
        """Return a copy of this Vector divided by its length, or None for a zero-length Vector"""
        length = self.length()
        if length == 0:
            return None
        return self / length

    def reversed(self) -> Self:
        """Return a copy of this Vector pointing in the opposite direction"""
        return Vector(-self.x, -self.y)

    def rotated(self, angle: float) -> Self:
        """Return a copy of this Vector rotated by the given angle"""
        return Vector(
            x = self.x * cos(-angle) - self.y * sin(-angle),
            y = self.x * sin(-angle) + self.y * cos(-angle),
        )
    
    def perpendicular(self) -> Self:
        """Return a Vector perpendicular to this Vector in the positive direction
        
        This is equivalent to `Vector.rotated(90)`."""
        return Vector(self.y, -self.x)
    
    def bisector(self, other: Self) -> Self:
        """Return the bisector between this Vector and the other given Vector, as a unit vector"""
        if isclose(self.length(), 0) or isclose(other.length(), 0):
            # If one of the vectors is null, we cannot calculate the bisector
            return None
        vectors_sum = self + other
        if vectors_sum.length() == 0:
            # If the two vectors are 180° from each other, there are two possible solutions : return
            # one of them arbitrarily
            return self.normalized().rotated(90)
        return vectors_sum.normalized()
    
    def angle_with(self, other: Self) -> float:
        """Calculate the angle between this Vector and the given one"""
        self_length = self.length()
        other_length = other.length()
        if isclose(self_length, 0) or isclose(other_length, 0):
            # Cannot calculate an angle of one of the vectors is null
            return None
        return acos(self.dot(other) / (self_length * other_length))
    
    def dot(self, other: Self) -> float:
        """Calculate the dot-product of this vector with the given vector"""
        return self.x * other.x + self.y * other.y
    
    def cross(self, other: Self) -> float:
        """Calculate the cross-product of this vector with the given vector"""
        return self.x * other.y - other.x * self.y


class Point(Vector, DrawableObject):
    """A point on the 2D plane, based on the Vector class"""

    def __str__(self) -> str:
        return f"Point(x={self.x:.2f}, y={self.y:.2f})"
    
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

        Supported objects : Point, Line, Segment (extrapolated as a line), Circle, Arc (extrapolated as a circle).
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
            
            case Segment():
                return self.distance(object.line())
            
            case Circle():
                return math.fabs(self.distance(object.center) - object.radius)
            
            case Arc():
                return math.fabs(self.distance(object.center()) - object.radius)

            case _:
                raise TypeError(f"Trying to compute the distance between a Point and an unsupported object : {type(object)}")
    
    def furthest(self, objects: list) -> Self:
        """Return the object in the given list that is furthest to this Point

        See distance() for the list of supported objects.
        The given objects do not need to be of the same type as long as they are
        each of one of the supported types.
        Throws a TypeError if any other object type is given.
        """
        max_distance = 0.0
        furthest = None
        for other in objects:
            distance = self.distance(other)
            if furthest is None or distance > max_distance:
                max_distance = distance
                furthest = other
        return furthest
    
    def closest(self, objects: list) -> Self:
        """Return the object in the given list that is closest to this Point

        See distance() for the list of supported objects.
        The given objects do not need to be of the same type as long as they are
        each of one of the supported types.
        Throws a TypeError if any other object type is given.
        """
        min_distance = 0.0
        closest = None
        if objects is None:
            return None
        if not isinstance(objects, (list, tuple)):
            objects = [objects]
        for other in objects:
            distance = self.distance(other)
            if closest is None or distance < min_distance:
                min_distance = distance
                closest = other
        return closest
    
    def translated(self, vector: Vector) -> Self:
        """Create a copy of this Point translated according to the given vector"""
        return Point(self.x + vector.x, self.y + vector.y)
    
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

        Supported objects : Line, Segment, Circle, Arc (extrapolated as a circle).
        Throws a TypeError if any other object type is given.
        """
        match object:
            case Line():
                return object.intersect(object.perpendicular(self))

            case Segment():
                line = object.line()
                return line.intersect(line.perpendicular(self))
            
            case Circle():
                projection_line = Line.from_two_points(object.center, self)
                return self.closest(projection_line.intersect(object))
            
            case Arc():
                return self.projected(object.as_circle())

            case _:
                raise TypeError(f"Trying to project a Point on an unsupported object : {type(object)}")
    
    def offset(self, along_line: "Line", distance: float) -> Self:
        """Create a new Point along the given Line separated by the given distance

        The given distance can be negative to flip the side of the offset.
        """
        return self + along_line.unit_vector() * distance

    def to_viewport(self) -> Self:
        return Vector(round(self.x, PRECISION), round(-self.y, PRECISION))
    
    def to_svg(self) -> str:
        p = self.to_viewport()
        return f"{round(p.x, PRECISION)},{round(p.y, PRECISION)}"
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement=None, radius=None, color=None, opacity=None, line_width=None, dashes=None) -> Self:
        """Draw this Point on the given SVG drawing
        
        This method returns self and can therefore be chained."""
        if parent is None:
            parent = drawing
        parent.add(drawing.circle(self.to_viewport().as_tuple(),
            radius or style.point_radius,
            stroke = color or style.point_color,
            stroke_opacity = opacity or style.point_opacity,
            stroke_width = line_width or style.point_line_width,
        ))
        return self
    
    def draw_kicad(self, kicadpcb: "KicadPCB", width: float, layer: str, stroke_type: str = 'solid') -> Self:
        """No-op : construction Points are not printed on the Kicad board"""
        return self


class Line(DrawableObject):
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
        if isclose(p1.y, p2.y):
            if isclose(p1.x, p2.x):
                print("Warning: trying to create a Line using two coincident points")
                return None
            else:
                a = math.inf
                b = p1.y
        else:
            a = (p1.x - p2.x) / (p1.y - p2.y)
            b = p1.x - a * p1.y
        return Line(a, b)

    def from_point_and_vector(point: Point, vector: Vector) -> Self:
        """Return a new Line based on a Point and a Vector"""
        return Line.from_two_points(point, point + vector)
    
    def horizontal(y: float) -> Self:
        """Create a new horizontal line at the given Y position"""
        return Line(math.inf, y)
    
    def as_base(self) -> Self:
        """Return the base geometric element of this Line, i.e. itself"""
        return self
    
    def unit_vector(self) -> Vector:
        """Return a new Vector of length 1 parallel to this line
        
        The direction of the vector is not guaranteed.
        """
        
        if math.isinf(self.a):
            # Horizontal line
            return Vector(1.0, 0.0)
        else:
            p1 = Point(self.eval(0), 0)
            p2 = self.intersect(Circle(p1, 1.0))[0]
            return Vector.from_two_points(p1, p2)
    
    def angle(self) -> float:
        """Return the angle of this Line in the range [-90°:90°] relative to the vertical axis"""
        if math.isinf(self.a):
            return 90.0
        else:
            return atan(self.a)
    
    def is_parallel(self, other: Self) -> bool:
        """Check if this Line is parallel with the given Line"""
        return isclose(self.a, other.a)
    
    def is_colinear(self, other: Self) -> bool:
        """Check if this Line is colinear with the given Line"""
        return isclose(self.a, other.a) and isclose(self.b, other.b)

    def contains_point(self, point: Point) -> bool:
        """Check if the given point is on this Line"""
        if math.isinf(self.a):
            return isclose(point.y, self.b)
        return isclose(point.x, self.a * point.y + self.b)
    
    def is_tangeant(self, object) -> bool:
        """Check if this Line is tangeant to the given object

        Supported objects : Circle, Arc (extrapolated as Circle).
        Throws a TypeError if any other object type is given."""
        match object:
            case Circle() | Arc():
                return isclose(self.distance(object), 0)

            case _:
                raise TypeError(f"Unsupported object given to Line.is_tangeant(), must be a Circle or Arc, received : {type(object)}")
    
    def distance(self, object) -> float:
        """Calculate the distance between this Line and the given object

        Supported objects : Point, Line, Circle, Arc.
        Throws a TypeError if any other object type is given.
        """
        match object:
            case Point():
                return object.distance(self)

            case Line():
                if self.is_parallel(object):
                    if math.isinf(self.a):
                        return math.fabs(self.b - object.b)
                    else:
                        p1 = Point(self.eval(0), 0)
                        p2 = p1.projected(object)
                        return p1.distance(p2)
                else:
                    # The lines are intersecting, the distance is considered to be zero
                    # because the lines necessarily intersect at some point
                    return 0
            
            case Circle():
                return math.fabs(self.distance(object.center) - object.radius)
            
            case Arc():
                # Extrapolate the Arc as a Circle
                return self.distance(object.as_circle())

            case _:
                raise TypeError(f"Trying to compute the distance between a Point and an unsupported object : {type(object)}")

    def intersect(self, object, suppress_warning: bool = False) -> Point:
        """Calculate the intersection Point(s) between this Line and the given object

        Supported objects : Line, Circle, Arc.
        Throws a TypeError if any other object type is given.
        Returns either a single Point, a tuple of two Points, or None if there is no
        intersection between the objects.
        """
        match object:
            case Line():
                if isclose(self.a, object.a):
                    if not suppress_warning:
                        print("Warning: trying to compute intersection between two parallel or colinear lines")
                    return None
                elif math.isinf(object.a):
                    return object.intersect(self, suppress_warning)
                elif math.isinf(self.a):
                    x = object.a * self.b + object.b
                    y = self.b
                    return Point(x, y)
                else:
                    if isclose(self.a, 0.0):
                        x = self.b
                        y = (x - object.b) / object.a
                    else:
                        x = (self.a * object.b - object.a * self.b) / (self.a - object.a)
                        y = (x - self.b) / self.a
                    return Point(x, y)
            
            case Segment():
                return self.intersect(object.line(), suppress_warning)

            case Circle():
                if self.is_tangeant(object):
                    # If the line is tangeant to the circle, find the intersection by projecting the center
                    # of the circle on the line to avoid float rounding errors, and return it as if it was
                    # the two solutions of the polynome
                    point = object.center.projected(self)
                    return (point, point)
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
                    if not suppress_warning:
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
                intersect = self.intersect(circle, suppress_warning)
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
        elif isclose(self.a, 0):
            return Line(math.inf, point.y)
        else:
            a = - 1 / self.a
            b = point.x + point.y / self.a
            return Line(a, b)

    def eval(self, y: float) -> float:
        if math.isinf(self.a):
            print("Warning: calling eval() on an horizontal Line")
        return self.a * y + self.b
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement=None, color=None, opacity=None, line_width=None, dashes=None) -> Self:
        """Draw this Line on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if parent is None:
            parent = drawing

        # Assume a sufficiently large viewport to cover most cases
        viewport_width = 1000
        viewport_height = 1000

        if math.isinf(self.a):
            x1 = -viewport_width/2.0
            x2 = viewport_width/2.0
            y = self.b
            p1 = (round(x1, PRECISION), round(-y, PRECISION))
            p2 = (round(x2, PRECISION), round(-y, PRECISION))
        else:
            y1 = -viewport_height/2.0
            y2 = viewport_height/2.0
            p1 = (round(self.eval(y1), PRECISION), round(-y1, PRECISION))
            p2 = (round(self.eval(y2), PRECISION), round(-y2, PRECISION))
        
        parent.add(drawing.line(
            p1, p2,
            stroke = color or style.line_color,
            stroke_opacity = opacity or style.line_opacity,
            stroke_width = line_width or style.line_width,
            stroke_dasharray = dashes or style.line_dashes,
        ))
        return self
    
    def draw_kicad(self, kicadpcb: "KicadPCB", width: float, layer: str, stroke_type: str = 'solid') -> Self:
        """No-op : construction Lines are not printed on the Kicad board"""
        return self


class Segment(DrawableObject):
    """A 2D line segment represented by two endpoints"""

    def __init__(self, p1: Point, p2: Point):
        self.p1: Point = p1
        self.p2: Point = p2

    def __str__(self) -> str:
        return f"Segment(p1={self.p1}, p2={self.p2})"

    def from_point_and_vector(point: Point, vector: Vector) -> Self:
        """Return a new Segment based on a Point and a Vector"""
        return Segment(point, point + vector)
    
    def copy(self) -> Self:
        """Return a copy of this Segment"""
        return Segment(self.p1, self.p2)
    
    def start_point(self) -> Point:
        """Return the start point (p1) of this Segment"""
        return self.p1
    
    def end_point(self) -> Point:
        """Return the end point (p2) of this Segment"""
        return self.p2
    
    def line(self) -> Line:
        """Return the Line colinear to this segment, losing the endpoints information"""
        return Line.from_two_points(self.p1, self.p2)
    
    def as_base(self) -> Line:
        """Return the base geometric element of this Segment, i.e. a Line"""
        return self.line()
    
    def as_vector(self) -> Vector:
        """Return a new Vector starting at p1 and going to p2"""
        return Vector.from_two_points(self.start_point(), self.end_point())
    
    def unit_vector(self) -> Vector:
        """Return a new Vector of unit length starting at p1 and pointing toward p2"""
        length = self.length()
        if not isclose(length, 0):
            return self.as_vector() / length
        else:
            # The start and end points are coincident, cannot calculate a unit vector
            return None
    
    def reversed(self) -> Self:
        """Return a copy of this Segment with its end points reversed"""
        return Segment(self.p2, self.p1)
    
    def length(self) -> float:
        """Return the length of this Segment"""
        return self.p1.distance(self.p2)

    def angle(self) -> float:
        """Return the angle of this Segment in the range [-90°:90°] relative to the vertical axis"""
        return self.line().angle()
    
    def translated(self, vector: Vector) -> Self:
        """Create a copy of this Segment translated according to the given vector"""
        return Segment(self.p1.translated(vector), self.p2.translated(vector))
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Segment rotated around the given center point by the given angle"""
        return Segment(self.p1.rotated(center, angle), self.p2.rotated(center, angle))
    
    def is_parallel(self, other: Self) -> bool:
        """Check if this Segment is parallel with the given Segment or Line"""
        match other:
            case Segment():
                return self.line().is_parallel(other.line())

            case Line():
                return self.line().is_parallel(other)

            case _:
                raise TypeError(f"Trying to check parallelism to a Segment on an unsupported object : {type(object)}")
                
    
    def is_colinear(self, other: Self) -> bool:
        """Check if this Segment is colinear with the given Segment or Line"""
        match other:
            case Segment():
                return self.line().is_colinear(other.line())

            case Line():
                return self.line().is_colinear(other)

            case _:
                raise TypeError(f"Trying to check colinearity to a Segment on an unsupported object : {type(object)}")

    def contains_point(self, point: Point) -> bool:
        """Check if the given point is on this Segment"""
        if point == self.p1 or point == self.p2:
            return True
        v1 = self.as_vector()
        v2 = Vector.from_two_points(self.p1, point)
        on_line = isclose(v1.cross(v2), 0.0)
        return on_line and v1.dot(v2) >= 0 and v2.length() <= v1.length()

    def midpoint(self) -> Point:
        """Return the point at the middle of this Segment"""
        return self.p1 + self.as_vector() * 0.5

    def tangeant_at_start(self) -> Vector:
        """Return the unit vector of this Segment
        
        This is provided as a convenience functions for other functions that handle both Segments and Arcs the same way."""
        return self.unit_vector()

    def tangeant_at_end(self) -> Vector:
        """Return the unit vector of this Segment
        
        This is provided as a convenience functions for other functions that handle both Segments and Arcs the same way."""
        return self.unit_vector()

    def intersect(self, object, suppress_warning: bool = False) -> Point:
        """Calculate the intersection Point(s) between this Segment (extrapolated as a Line) and the given object

        Supported objects : see Line.
        Throws a TypeError if any other object type is given.
        Returns either a single Point, a tuple of two Points, or None if there is no
        intersection between the objects.
        """
        return self.line().intersect(object, suppress_warning)
    
    def cut(self, tool, from_end: bool = False) -> Self:
        """Return a copy of this Segment cut by the given "tool"

        Supported objects for the tool : Line, Segment, Circle, Arc.
        The part of the segment that is kept depends on the `from_end` parameter.
        """
        if tool is None:
            # Nothing to cut with
            return self.copy()
        start_point = self.end_point() if from_end else self.start_point()
        intersection_points = self.line().intersect(tool.as_base(), suppress_warning=True)
        if intersection_points is None:
            return self.copy()
        if not isinstance(intersection_points, (list, tuple)):
            intersection_points = [intersection_points]
        cut_segments = []
        for intersection_point in intersection_points:
            if isinstance(tool, (Segment, Arc)) and not tool.contains_point(intersection_point):
                # The tool is a bounded object (Segment or Arc) and the intersection point is not inside it,
                # ignore this point
                continue
            v_self = self.as_vector()
            if from_end:
                v_self = v_self.reversed()
            v_intersect = Vector.from_two_points(start_point, intersection_point)
            if self.contains_point(intersection_point) and intersection_point != start_point:
                # The line or segment intersects the segment : cut it
                segment = Segment(start_point, intersection_point)
                if from_end:
                    cut_segments.append(segment.reversed())
                else:
                    cut_segments.append(segment)
            elif v_self.dot(v_intersect) <= 0:
                # The intersection point is before the start point : assume the segment is complety cut
                return None
            else:
                # The intersection point is after the end point
                continue
        if cut_segments:
            # At least one cut found : return the smallest one
            return min(cut_segments, key = lambda v: v.length())
        else:
            # No cut found : return the segment as-is
            return self.copy()
    
    def offset(self, distance: float, distance2: float=None) -> Self:
        """Create a new Segment parallel to this Segment separated by the given distance

        The given distance can be negative to flip the side of the offset.
        """
        if self.p1 == self.p2:
            print("Warning : cannot offset a zero-length Segment")
            return None
        if distance2 is None:
            distance2 = distance
        p1_offset = self.start_point() + self.tangeant_at_start().perpendicular() * distance
        p2_offset = self.end_point() + self.tangeant_at_end().perpendicular() * distance2
        return Segment(p1_offset, p2_offset)
    
    def offset_closest_to(self, point: Point, distance: float) -> Self:
        """Create a new Segment parallel to this Segment separated by the given distance and is the closest to the given Point

        The sign of the given offset is ignored, the function will compute the offsets in the two directions and return the
        offset segment closest to the given Point.
        """
        return point.closest([self.offset(distance), self.offset(-distance)])

    def tangent_arc_through_point(self, point: Point, at_start: bool = False) -> "Arc":
        """Return an Arc tangeant to this Segment at one of its endpoint that passes through the given Point

        By default, the new Arc is placed at the end of the current Segment, except if `at_start` is set.
        """
        if at_start:
            p1 = self.p1 # Start the new arc at the start of the segment
            opposite_point = self.p2
        else:
            p1 = self.p2 # Start the new arc at the end of the segment
            opposite_point = self.p1
        p2 = point
        if p1 == p2:
            # The endpoint of the segment and the target point are coincident, the arc would be empty
            return None
        v1 = Vector.from_two_points(opposite_point, p1)
        v2 = Vector.from_two_points(p1, p2)
        dot_product = v1.dot(v2)
        if dot_product < 0:
            # The target point is behind the endpoint of the segment, it is not possible to draw an arc
            # tangeant to the segment that would pass by the given point
            return None
        cross_product = v1.cross(v2)
        if isclose(cross_product, 0.0):
            # The target point is coincident with the segment, the result would be a line
            # TODO : allow the Arc to have a radius=math.inf for these cases?
            return None
        segment = Segment(p1, p2)
        midline = segment.line().perpendicular(segment.midpoint())
        center = self.line().perpendicular(p1).intersect(midline)
        radius = center.distance(p1)
        return Arc(p1, p2, radius, cross_product > 0)
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement=None, color=None, opacity=None, line_width=None, dashes=None) -> Self:
        """Draw this Segment on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if parent is None:
            parent = drawing

        parent.add(drawing.line(
            self.p1.to_viewport().as_tuple(), self.p2.to_viewport().as_tuple(),
            stroke = color or style.line_color,
            stroke_opacity = opacity or style.line_opacity,
            stroke_width = line_width or style.line_width,
            stroke_dasharray = dashes or style.line_dashes,
        ))
        return self
    
    def draw_kicad(self, kicadpcb: "KicadPCB", width: float, layer: str, stroke_type: str = 'solid') -> Self:
        """Draw this Segment on the given Kicad board"""
        kicadpcb.gr_line(
            p1 = self.p1,
            p2 = self.p2,
            width = width,
            layer = layer,
            stroke_type = stroke_type,
        )
        return self

class Circle(DrawableObject):
    """A 2D circle represented by its center and radius"""

    def __init__(self, center: Point, radius: float):
        self.center: Point = center
        self.radius: float = radius

    def __str__(self) -> str:
        return f"Circle(center={self.center}, radius={self.radius:.2f})"
    
    def as_base(self) -> Self:
        """Return the base geometric element of this Circle, i.e. itself"""
        return self
    
    def translated(self, vector: Vector) -> Self:
        """Create a copy of this Circle translated according to the given vector"""
        return Circle(self.center.translated(vector), self.radius)
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Circle rotated around the given center point by the given angle"""
        return Circle(self.center.rotated(center, angle), self.radius)
    
    def perimeter(self) -> float:
        """Return the perimeter of this Circle"""
        return 2 * math.pi * self.radius

    def contains_point(self, point: Point) -> bool:
        """Check if the given point is on this Circle"""
        if point == self.p1 or point == self.p2:
            return True
        return isclose(self.center().distance(point), self.radius)
    
    def is_tangeant(self, object) -> bool:
        """Check if this Circle is tangeant to the given object

        Supported objects : Line, Circle, Arc (extrapolated as Circle).
        Throws a TypeError if any other object type is given."""
        match object:
            case Line():
                return object.is_tangeant(self)

            case Circle():
                if isclose(self.center.distance(object.center), math.fabs(self.radius - object.radius)):
                    # The circles are tangeant to one another, with one inside the other
                    return True
                elif isclose(self.center.distance(object.center), math.fabs(self.radius + object.radius)):
                    # The circles are tangeant to one another, outside of each other
                    return True
                return False

            case Arc():
                return self.is_tangeant(object.as_circle())

            case _:
                raise TypeError(f"Unsupported object given to Line.is_tangeant(), must be a Circle or Arc, received : {type(object)}")

    def intersect(self, object, suppress_warning: bool = False) -> Point:
        """Calculate the intersection Point between this Circle and the given object

        Supported objects : Line, Segment, Circle.
        Throws a TypeError if any other object type is given.
        Returns either a single Point, a tuple of two Points, or None if there is no
        intersection between the objects.
        """
        match object:
            case Line() | Segment():
                return object.intersect(self, suppress_warning)
            
            case Circle():
                if self.center == object.center:
                    # Concentric circles : they may be identical (infinite intersection) or not,
                    # but in both cases return None
                    if not suppress_warning:
                        print("Warning: trying to compute intersection between two concentric circles")
                    return None
                elif self.is_tangeant(object):
                    # The two circles are tangeant : intersect the line connecting the center points with the bigger
                    # circle, and take the intersection that is closest to the center of the smaller circle
                    if object.radius > self.radius:
                        bigger_circle = object
                        smaller_circle = self
                    else:
                        bigger_circle = self
                        smaller_circle = object
                    intersect = smaller_circle.center.closest(Line.from_two_points(self.center, object.center).intersect(bigger_circle, suppress_warning))
                    return intersect
                distance = math.hypot(object.center.x - self.center.x, object.center.y - self.center.y)
                if distance > self.radius + object.radius or distance < math.fabs(self.radius - object.radius):
                    # No intersection
                    if not suppress_warning:
                        print("Warning: trying to compute intersection between non-intersecting circles")
                    return None
                a = (self.radius**2 - object.radius**2 + distance**2) / (2.0 * distance)
                x = self.radius**2 - a**2
                if isclose(x, 0):
                    # Prevent out-of-domain errors due to float rounding errors
                    h = 0
                else:
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
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement=None, color=None, opacity=None, line_width=None, dashes=None) -> Self:
        """Draw this Circle on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if parent is None:
            parent = drawing

        parent.add(drawing.circle(
            self.center.to_viewport().as_tuple(),
            self.radius,
            stroke = color or style.circle_color,
            stroke_opacity = opacity or style.circle_opacity,
            stroke_width = line_width or style.circle_line_width,
            stroke_dasharray = dashes or style.circle_dashes,
        ))
        return self
    
    def draw_kicad(self, kicadpcb: "KicadPCB", width: float, layer: str, stroke_type: str = 'solid') -> Self:
        """Draw this Circle on the given Kicad board"""
        kicadpcb.gr_circle(
            center = self.center,
            radius = self.radius,
            width = width,
            layer = layer,
            stroke_type = stroke_type,
        )
        return self


class Arc(DrawableObject):
    """A 2D Arc represented by its two endpoints and its radius

    The arc always goes clockwise from p1 to p2.
    """

    def __init__(self, p1: Point, p2: Point, radius: float, reverse: bool = False, suppress_warning: bool = False):
        """Create a new Arc using the to given endpoints and the given radius

        If reverse is set, p1 and p2 are swapped before constructing the Arc.
        """
        if p1 == p2:
            raise ValueError("Error : unable to create an arc using two coincident points")
        distance = p1.distance(p2)
        min_radius = distance / 2.0
        if radius < min_radius:
            if not isclose(radius, min_radius) and not suppress_warning:
                print(f"Warning : trying to create an Arc with a radius (r={radius}) too small for the distance between the two points (d={distance}), setting radius to {min_radius}")
            radius = min_radius
        if reverse:
            p1, p2 = p2, p1
        self.p1: Point = p1
        self.p2: Point = p2
        self.radius: Point = radius
        self._reversed = reverse

    def __str__(self) -> str:
        return f"Arc(p1={self.p1}, p2={self.p2}, radius={self.radius:.2f}{', reversed' if self._reversed else ''})"
    
    def __mul__(self, factor: float) -> Self:
        """Create a new Arc concentric to this Arc with the same p1, and with
        p2 rotated around the center so that the length of the arc is grown or
        shrinked by the given factor.
        """
        return Arc(self.p1, self.p1.rotated(self.center(), self.angle() * factor), self.radius)
    
    def copy(self) -> Self:
        """Return a copy of this Arc"""
        arc = Arc(self.p1, self.p2, self.radius)
        arc._reversed = self._reversed
        return arc

    def start_point(self) -> Point:
        """Return the start point of this Arc
        
        Depending on the `reverse` flag provided during the creation of the Arc, this is either p1 (`reverse=False`) or p2 (`reverse=False`)."""
        if self._reversed:
            return self.p2
        else:
            return self.p1
    
    def end_point(self) -> Point:
        """Return the end point of this Arc
        
        Depending on the `reverse` flag provided during the creation of the Arc, this is either p2 (`reverse=False`) or p1 (`reverse=False`)."""
        if self._reversed:
            return self.p1
        else:
            return self.p2
    
    def translated(self, vector: Vector) -> Self:
        """Create a copy of this Arc translated according to the given vector"""
        return Arc(self.p1.translated(vector), self.p2.translated(vector), self.radius, self.reverse)
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Segment rotated around the given center point by the given angle"""
        return Arc(self.p1.rotated(center, angle), self.p2.rotated(center, angle), self.radius, self.reverse)
    
    def center(self) -> Point:
        """Calculate the center Point of this Arc"""

        # Check for edge conditions
        if isclose(self.p1.distance(self.p2), self.radius*2, abs_tol=1e-6):
            # The arc spans 180° and finding the intersection between two circles might not work,
            # depending on float rounding errors. Instead, calculate the center of the arc by taking
            # the midpoint of the two endpoints.
            return Segment(self.p1, self.p2).midpoint()

        # Otherwise, find the center by calculating the point at the same distance d=radius from the two endpoints
        circle1 = Circle(self.p1, self.radius)
        circle2 = Circle(self.p2, self.radius)
        intersect = circle1.intersect(circle2)
        match intersect:
            case None:
                raise ValueError(f"Cannot calculate the center of arc p1={self.p1} p2={self.p2} r={self.radius}")
            
            case Point():
                return intersect
            
            case [c0, c1]:
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
    
    def as_circle(self) -> Circle:
        """Return the Circle that this Arc is part of"""
        return Circle(self.center(), self.radius)
    
    def as_base(self) -> Circle:
        """Return the base geometric element of this Arc, i.e. a Circle"""
        return self.as_circle()

    def contains_point(self, point: Point) -> bool:
        """Check if the given point is on this Arc"""
        if point == self.p1 or point == self.p2:
            return True
        center = self.center()
        if not isclose(center.distance(point), self.radius):
            return False
        v1 = Vector.from_two_points(center, self.p1)
        v2 = Vector.from_two_points(center, self.p2)
        v3 = Vector.from_two_points(center, point)
        return v3.cross(v1) >= 0 and v1.dot(v3) >= v1.dot(v2)

    def midpoint(self) -> Point:
        """Return the point at the middle of this Arc"""
        return self.p1.rotated(self.center(), self.angle() / 2.0)
    
    def distance(self, object) -> float:
        """Calculate the distance between this Arc and the given object

        Supported objects : Point, Line.
        Throws a TypeError if any other object type is given.
        """
        match object:
            case Point():
                return object.distance(self)

            case Line():
                return object.distance(self)

            case _:
                raise TypeError(f"Trying to compute the distance between a Point and an unsupported object : {type(object)}")
    
    def closest(self, objects: list) -> Self:
        """Return the object in the given list that is closest to this Arc

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
    
    def is_tangeant(self, object) -> bool:
        """Check if this Arc (extrapolated as a Circle) is tangeant to the given object

        Supported objects : see Circle.is_tangeant().
        Throws a TypeError if any other object type is given."""
        return self.as_circle().is_tangeant(object)

    def tangeant_at_start(self) -> Vector:
        """Return a unit Vector tangeant to the start point of this Arc, oriented in the same direction as the Arc"""
        return Segment(self.center(), self.start_point()).unit_vector().rotated(-90 if self._reversed else 90)

    def tangeant_at_point(self, point: Point) -> Vector:
        """Return a unit Vector tangeant to this Arc at the given point, oriented in the same direction as the Arc"""
        center = self.center()
        if point == center:
            return None
        if not self.contains_point(point):
            # If the point is not on the arc, project it
            point = point.closest(self.intersect(Line.from_two_points(center, point)))
        return Segment(center, point).unit_vector().rotated(-90 if self._reversed else 90)

    def tangeant_at_end(self) -> Vector:
        """Return a unit Vector tangeant to the end point of this Arc, oriented in the same direction as the Arc"""
        return Segment(self.center(), self.end_point()).unit_vector().rotated(-90 if self._reversed else 90)

    def intersect(self, object, suppress_warning: bool = False) -> Point:
        """Calculate the intersection Point between this Arc (extrapolated as a Circle) and the given object

        Supported objects : Line, Segment, Circle, Arc (extrapolated as a Circle).
        Throws a TypeError if any other object type is given.
        Returns either a single Point, a tuple of two Points, or None if there is no
        intersection between the objects.
        """
        match object:
            case Line() | Segment():
                return object.intersect(self, suppress_warning)

            case Circle():
                return object.intersect(self.as_circle(), suppress_warning)

            case Arc():
                return self.as_circle().intersect(object.as_circle(), suppress_warning)

            case _:
                raise TypeError(f"Unsupported object given to intersect an Arc : {type(object)}")
    
    def cut(self, tool, from_end: bool = False) -> Self:
        """Return a copy of this Arc cut by the given "tool"
        
        Supported objects for the tool : Line, Segment, Circle, Arc.
        The part of the arc that is kept depends on the `from_end` parameter.
        """
        if tool is None:
            # Nothing to cut with
            return self.copy()
        start_point = self.end_point() if from_end else self.start_point()
        # end_point = self.start_point() if from_end else self.end_point()
        intersection_points = self.as_circle().intersect(tool.as_base(), suppress_warning=True)
        if intersection_points is None:
            return self.copy()
        if not isinstance(intersection_points, (list, tuple)):
            intersection_points = [intersection_points]
        center = self.center()
        cut_arcs = []
        # v_start_point = Vector.from_two_points(center, self.start_point())
        # v_end_point = Vector.from_two_points(center, self.end_point())
        for intersection_point in intersection_points:
            if isinstance(tool, (Segment, Arc)) and not tool.contains_point(intersection_point):
                # The tool is a bounded object (Segment or Arc) and the intersection point is not inside it,
                # ignore this point
                continue
            v_intersection_point = Vector.from_two_points(center, intersection_point)
            if self.contains_point(intersection_point) and intersection_point != start_point:
                # The intersection point is on the current arc : create the cut arc and add to the list
                if from_end:
                    arc = Arc(intersection_point, start_point, self.radius, reverse=self._reversed)
                else:
                    arc = Arc(start_point, intersection_point, self.radius, reverse=self._reversed)
                cut_arcs.append(arc)
            else:
                # Unline for segments, on a circle there is no definitive way to tell if a point is "before"
                # or "after" another one. We can either look at which side the intersection point is on
                # relative to the start point, or just never cut the arc to zero if the intersection point
                # is outside it. Right now, we do the latter.
                continue
                # if v_start_point.cross(v_intersection_point) * v_start_point.cross(v_end_point) <= 0:
                #     # The intersection point is "before" the start point
                #     return None
                # else:
                #     # The intersection point is "after" the end point
                #     continue
        if cut_arcs:
            # At least one cut found : return the smallest one
            return min(cut_arcs, key = lambda v: v.length())
        else:
            # No cut found : return the arc as-is
            return self.copy()
    
    def radial_line(self, point: Point) -> Line:
        """Create a new Line radial to this Arc that passes through the given Point"""
        return Line.from_two_points(self.center(), point)

    def offset(self, distance: float, distance2: float = None) -> Self:
        """Create a new Arc concentric with this Arc offset by the given distance

        The given distance can be negative to flip the side of the offset.
        """
        # TODO : when distance2 != distance, replace this approximation with an ellipse
        if self.p1 == self.p2:
            print("Warning : cannot offset a zero-length Arc")
            return None
        if distance2 is None:
            distance2 = distance
        if self._reversed:
            radius_offset_1 = self.radius + distance
            radius_offset_2 = self.radius + distance2
        else:
            radius_offset_1 = self.radius - distance
            radius_offset_2 = self.radius - distance2
        if radius_offset_1 <= 0 or radius_offset_2 <= 0:
            # The distance is too large to compute an offset for this Arc
            return None
        p1_offset = self.start_point() + self.tangeant_at_start().perpendicular() * distance
        p2_offset = self.end_point() + self.tangeant_at_end().perpendicular() * distance2
        if p1_offset == p2_offset:
            # An arc cannot be created with two coincident points
            return None
        return Arc(p1_offset, p2_offset, (radius_offset_1 + radius_offset_2) / 2, self._reversed)

    def tangent_arc_through_point(self, point: Point, at_start: bool = False) -> Self:
        """Return an Arc tangeant to this Arc at one of its endpoint that passes through the given Point

        By default, the new Arc is placed at the end of the current Arc, except if `at_start` is set.
        """
        if at_start:
            p1 = self.p1 # Start the new arc at the start of the current arc
            opposite_point = self.p2
        else:
            p1 = self.p2 # Start the new arc at the end of the current arc
            opposite_point = self.p1
        p2 = point
        if p1 == p2:
            # The endpoint of the current arc and the target point are coincident, the new arc would be empty
            return None
        line = Line.from_two_points(self.center(), p1).perpendicular(p1)
        p_segment = opposite_point.closest(line.intersect(Circle(p1, 1.0)))
        tangent_segment = Segment(p_segment, p1)
        v1 = Vector.from_two_points(p_segment, p1)
        v2 = Vector.from_two_points(p1, p2)
        dot_product = v1.dot(v2)
        if dot_product < 0:
            # The target point is behind the endpoint of the current arc, cannot calculate a new arc tangeant to the end of this arc
            return None
        cross_product = v1.cross(v2)
        if isclose(cross_product, 0.0):
            # The target point is coincident with the line tangeant at the end of the arc, cannot calculate a finite arc
            return None
        segment = Segment(p1, p2)
        midline = segment.line().perpendicular(segment.midpoint())
        center = line.perpendicular(p1).intersect(midline)
        radius = center.distance(p1)
        return Arc(p1, p2, radius, cross_product > 0)
    
    def discretize(self, segment_length: float = 1, max_deviation: float = 5) -> list[Point]:
        """Return an approximation of this Arc as a series of small segments represented by a list of Points"""
        start_point = self.start_point()
        points = [start_point]
        current_point = start_point
        previous_point = start_point + self.tangeant_at_start().reversed()
        
        for i in range(1000):
            # Find the next point on the arc at the segment_length distance from the current point
            next_point_by_length = self.end_point()
            circle = Circle(center=current_point, radius=segment_length)
            intersection_points = self.as_circle().intersect(circle, suppress_warning=True)
            if intersection_points is not None:
                if not isinstance(intersection_points, (list, tuple)):
                    intersection_points = [intersection_points]
                point = previous_point.furthest(intersection_points)
                if self.contains_point(point):
                    next_point_by_length = point

            # Find the next point on the arc at the max_deviation angle from the current point
            next_point_by_angle = self.end_point()
            line = Line.from_point_and_vector(current_point, self.tangeant_at_point(current_point).rotated(-max_deviation if self._reversed else max_deviation))
            intersection_points = self.as_circle().intersect(line, suppress_warning=True)
            if intersection_points is not None:
                if not isinstance(intersection_points, (list, tuple)):
                    intersection_points = [intersection_points]
                point = previous_point.furthest(intersection_points)
                if self.contains_point(point):
                    next_point_by_angle = point
            
            # Take the closest of these two points
            next_point = current_point.closest([next_point_by_length, next_point_by_angle])
            
            # Add this point to the list
            if next_point == self.end_point() or not self.contains_point(next_point):
                break
            points.append(next_point)
            previous_point = current_point
            current_point = next_point
        
        points.append(self.end_point())
        return points
    
    def draw_svg(self, drawing: svg.Drawing, parent: svg.base.BaseElement=None, color=None, opacity=None, line_width=None, dashes=None) -> Self:
        """Draw this Arc on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if parent is None:
            parent = drawing

        parent.add(drawing.path(
            d=f"M{self.p1.to_svg()} A{round(self.radius, PRECISION)},{round(self.radius, PRECISION)} 0 0,1 {self.p2.to_svg()}",
            stroke = color or style.arc_color,
            stroke_opacity = opacity or style.arc_opacity,
            stroke_width = line_width or style.arc_line_width,
            stroke_dasharray = dashes or style.arc_dashes,
        ))
        return self
    
    def draw_kicad(self, kicadpcb: "KicadPCB", width: float, layer: str, stroke_type: str = 'solid') -> Self:
        """Draw this Arc on the given Kicad board"""
        kicadpcb.gr_arc(
            p1 = self.p1,
            p2 = self.p2,
            midpoint = self.midpoint(),
            width = width,
            layer = layer,
            stroke_type = stroke_type,
        )
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

    def segment_to_segment(p_prev: Point, p_mid: Point, p_next: Point, fillet_radius: float, suppress_warning: bool = False) -> Self:
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
        warning_displayed = False
        for i in range(3):
            # Compute the center of the fillet by offsetting each of the two lines in the direction of the opposite point,
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
                if not warning_displayed and not suppress_warning:
                    print("Warning : segment-to-segment fillet too big for the previous line, cutting the fillet")
                    warning_displayed = True
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
                if not warning_displayed and not suppress_warning:
                    print("Warning : segment-to-segment fillet too big for the next line, cutting the fillet")
                    warning_displayed = True
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
            
            # The points are kept in the order of the path, so we must return a separate flag
            # to indicate an anticlockwise arc direction
            anticlockwise_fillet_arc = v_prev.cross(v_next) < 0

            # Return the two endpoints of the fillet, its finale radius, and the anticlockwise flag
            return Fillet(p_fillet_prev, p_fillet_next, fillet_radius, anticlockwise_fillet_arc, shrinked)

    def segment_to_arc(
            p_prev: Point,
            p_mid: Point,
            p_next: Point,
            arc_radius: float,
            anticlockwise: bool,
            fillet_radius: float,
            suppress_warning: bool = False,
        ) -> Self:
        """Compute the Fillet between a segment connected to an arc of the given radius

        The segment and the arc are represented respectively by (p_prev, p_mid) and (p_mid, p_next),
        and therefore share the p_mid point. The fillet is applied on this p_mid point.
        Set the `anticlockwise` flag if the arc goes anti-clockwise from p_mid to p_next.
        If the fillet is too large for the the given geometry, its radius will be
        automatically shrinked to allow it to fit the geometry. Check the `shrinked`
        attribute on the returned Fillet object if necessary.
        """

        warning_displayed = False

        # Shrink the fillet radius a first time if it is larger than the segment
        segment_length = p_prev.distance(p_mid)
        if fillet_radius > segment_length:
            if not warning_displayed and not suppress_warning:
                print("Warning : segment-to-arc fillet too big for the given geometry, cutting the fillet")
                warning_displayed = True
            fillet_radius = segment_length

        # Compute the line and arc adjacent to the mid point
        segment_prev = Segment(p_prev, p_mid)
        arc_next = Arc(p_mid, p_next, arc_radius, reverse=anticlockwise)

        # If the line is tangeant to the arc, ignore the fillet
        if segment_prev.tangeant_at_end() == arc_next.tangeant_at_start():
            return None
        
        shrinked = False
        for i in range(3):
            # Compute the center of the fillet by offsetting each of the two objects, then taking the intersection
            # between them.
            # If no intersection is found, shrink the fillet until one can be successfully computed.
            while fillet_radius > 0.0:
                # The side to apply the fillet on depends on the geometry. Even if the next point is on the
                # right of the current segment, we cannot assume that the fillet will point to the right like
                # we would when filleting between two segments : if the radius of the arc is small enough, the
                # arc might start to go toward the left of the segment, then turn to the right to connect to
                # p_next.
                # We check for this condition by computing the tangeant unit vector at the start point of the
                # arc, to determine whether the arc starts in the direction of the next point.
                segment_prev_offset_1 = segment_prev.offset(fillet_radius)
                segment_prev_offset_2 = segment_prev.offset(-fillet_radius)
                vec_prev = Segment(p_prev, p_mid).as_vector()
                vec_tangeant_start_arc = arc_next.tangeant_at_start()
                vec_to_next = Segment(p_prev, p_next).as_vector()
                if isclose(vec_prev.cross(vec_to_next), 0):
                    # The next point is aligned with the previous segment, we cannot use p_next.closest() : select
                    # the correct offset based on the direction of the tangeant relative to the segment
                    if vec_prev.cross(vec_tangeant_start_arc) >= 0:
                        segment_prev_offset = segment_prev_offset_2
                    else:
                        segment_prev_offset = segment_prev_offset_1
                elif vec_prev.cross(vec_tangeant_start_arc) * vec_prev.cross(vec_to_next) >= 0:
                    # The cross products are of the same sign so the arc starts in the direction of the next
                    # point. Select the line closest to the next point to orient the fillet toward the inside.
                    # We check for >= 0 because vec_prev.cross(vec_tangeant_start_arc) might equal zero : this
                    # means that the arc is tangeant to the line but the path is making a U-turn, and the line
                    # closest to the next point should be selected in this case as well.
                    segment_prev_offset = p_next.closest([segment_prev_offset_1, segment_prev_offset_2])
                else:
                    # Otherwise, select the line furthest from the next point to orient the fillet toward the outside
                    segment_prev_offset = p_next.furthest([segment_prev_offset_1, segment_prev_offset_2])

                # Once we have the correct offset line on the side that the fillet needs to point to, we can
                # calculate the center of the fillet by intersecting with the offset arc.
                arc_next_offset_1 = arc_next.offset(fillet_radius)
                if arc_next_offset_1 is not None:
                    intersect_1 = p_mid.closest(segment_prev_offset.intersect(arc_next_offset_1.as_circle(), suppress_warning=True))
                else:
                    intersect_1 = None
                arc_next_offset_2 = arc_next.offset(-fillet_radius)
                if arc_next_offset_2 is not None:
                    intersect_2 = p_mid.closest(segment_prev_offset.intersect(arc_next_offset_2.as_circle(), suppress_warning=True))
                else:
                    intersect_2 = None
                if intersect_1 is None and intersect_2 is None:
                    fillet_center = None
                elif intersect_1 is None:
                    fillet_center = intersect_2
                elif intersect_2 is None:
                    fillet_center = intersect_1
                else:
                    fillet_center = p_prev.closest([intersect_1, intersect_2])
                if fillet_center is not None:
                    break
                else:
                    if not warning_displayed and not suppress_warning:
                        print("Warning : segment-to-arc fillet too big for the given geometry")
                        warning_displayed = True
                    fillet_radius -= 0.1
            if isclose(fillet_radius, 0.0):
                return None

            # Project the center point of the fillet on the adjacent lines to compute the endpoints of the arc of the fillet
            p_fillet_prev = fillet_center.projected(segment_prev)
            p_fillet_next = fillet_center.projected(arc_next)
            if p_fillet_prev == p_fillet_next:
                return None

            # Check if these endpoints are inside the adjacent segments
            v_prev = Vector.from_two_points(p_mid, p_prev)
            v_fillet_prev = Vector.from_two_points(p_mid, p_fillet_prev)
            arc_fillet_next = Arc(p_mid, p_fillet_next, arc_next.radius, reverse=anticlockwise)
            cut = False
            # TODO : replace the 'factor' calculation with something more robust
            if v_fillet_prev.length() >= v_prev.length() or v_prev.dot(v_fillet_prev) < 0:
                # Compute a new pair of fillet endpoints that are inside the prev segment
                if not warning_displayed and not suppress_warning:
                    print("Warning : segment-to-arc fillet too big for the previous line, cutting the fillet")
                    warning_displayed = True
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
                if not warning_displayed and not suppress_warning:
                    print("Warning : segment-to-arc fillet too big for the next line, cutting the fillet")
                    warning_displayed = True
                length_1 = arc_fillet_next.length()
                if anticlockwise:
                    arc_fillet_next = Arc(arc_next.p2.rotated(arc_next.center(), - arc_next.angle() * 0.9), arc_next.p2, arc_next.radius)
                    p_fillet_next = arc_fillet_next.p1
                else:
                    arc_fillet_next = Arc(arc_next.p1, arc_next.p1.rotated(arc_next.center(), arc_next.angle() * 0.9), arc_next.radius)
                    p_fillet_next = arc_fillet_next.p2
                length_2 = arc_fillet_next.length()
                factor = length_2 / length_1
                v_fillet_prev = v_fillet_prev * factor
                p_fillet_prev = p_mid + v_fillet_prev
                cut = True
                shrinked = True
            if cut:
                # The fillet was modified : compute its new radius and start again
                line_prev = segment_prev.line()
                l1 = line_prev.perpendicular(p_fillet_prev)
                l2 = arc_next.radial_line(p_fillet_next)
                new_center = l1.intersect(l2)
                fillet_radius = min(new_center.distance(line_prev), new_center.distance(arc_next))
                # Arc(p_fillet_prev, p_fillet_next, fillet_radius, reverse=anticlockwise)
                continue

            # The points are kept in the order of the path, so we must return a separate flag
            # to indicate an anticlockwise arc direction
            anticlockwise_arc = Vector.from_two_points(p_mid, p_fillet_prev).cross(Vector.from_two_points(p_mid, p_fillet_next)) < 0

            # Return the two endpoints of the fillet, its finale radius, and the anticlockwise flag
            return Fillet(p_fillet_prev, p_fillet_next, fillet_radius, anticlockwise_arc, shrinked)

    def arc_to_segment(
            p_prev: Point,
            p_mid: Point,
            p_next: Point,
            arc_radius: float,
            anticlockwise: bool,
            fillet_radius: float,
            suppress_warning: bool = False,
        ) -> Self:
        """Compute the Fillet between an arc of the given radius connected to a segment

        The arc and the segment are represented respectively by (p_prev, p_mid) and (p_mid, p_next),
        and therefore share the p_mid point. The fillet is applied on this p_mid point.
        Set the `anticlockwise` flag if the arc goes anti-clockwise from p_prev to p_mid.
        If the fillet is too large for the the given geometry, its radius will be
        automatically shrinked to allow it to fit the geometry. Check the `shrinked`
        attribute on the returned Fillet object if necessary.
        """

        # Create the equivalent Fillet from the segment to the arc, and reverse it
        fillet = Fillet.segment_to_arc(p_next, p_mid, p_prev, arc_radius, not anticlockwise, fillet_radius, suppress_warning)
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

    def __init__(self, end_point: Point, width: float = None, tag = None):
        self.p2: Point = end_point
        self.width = width
        self.tag = tag

    def __str__(self) -> str:
        return f"PathSegment(p2={self.p2}, width={self.width:.2f})"
    
    def translated(self, vector: Vector) -> Self:
        """Create a copy of this PathSegment translated according to the given vector"""
        return PathSegment(self.p2.translated(vector), width=self.width, tag=self.tag)
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this PathSegment rotated around the given center point by the given angle"""
        return PathSegment(self.p2.rotated(center, angle), width=self.width, tag=self.tag)
    
    def mirrored_x(self) -> Self:
        """Create a copy of this PathSegment mirrored about the X axis"""
        return PathSegment(self.p2.mirrored_x(), width=self.width, tag=self.tag)
    
    def mirrored_y(self) -> Self:
        """Create a copy of this PathSegment mirrored about the Y axis"""
        return PathSegment(self.p2.mirrored_y(), width=self.width, tag=self.tag)
    
    def geometry(self, p1: Point) -> Segment:
        return Segment(p1, self.p2)


class PathArc(PathElement):
    """An arc segment as part of a Path"""

    def __init__(self, end_point: Point, radius: float, anticlockwise: bool, width: float = None, tag = None):
        self.p2: Point = end_point
        self.radius: Point = radius
        self.anticlockwise: Point = anticlockwise
        self.width = width
        self.tag = tag

    def __str__(self) -> str:
        return f"PathArc(p2={self.p2}, radius={self.radius:.2f}, anticlockwise={self.anticlockwise}, width={self.width:.2f})"
    
    def translated(self, vector: Vector) -> Self:
        """Create a copy of this PathArc translated according to the given vector"""
        return PathArc(self.p2.translated(vector), self.radius, self.anticlockwise, width=self.width, tag=self.tag)
    
    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this PathArc rotated around the given center point by the given angle"""
        return PathArc(self.p2.rotated(center, angle), self.radius, self.anticlockwise, width=self.width, tag=self.tag)
    
    def mirrored_x(self) -> Self:
        """Create a copy of this PathArc mirrored about the X axis"""
        return PathArc(self.p2.mirrored_x(), self.radius, not self.anticlockwise, width=self.width, tag=self.tag)
    
    def mirrored_y(self) -> Self:
        """Create a copy of this PathArc mirrored about the Y axis"""
        return PathArc(self.p2.mirrored_y(), self.radius, not self.anticlockwise, width=self.width, tag=self.tag)
    
    def geometry(self, p1: Point) -> Arc:
        return Arc(p1, self.p2, self.radius, self.anticlockwise)

class Path(DrawableObject):
    """A Path consisting of a list of connected segments and arcs"""

    def __init__(self, start_point: Point, start_point_width: float = None, width: float = 1):
        self.start_point: Point = start_point
        self.start_point_width: float = start_point_width or width
        self.width: float = width
        self.elements: list[PathElement] = []

    def __str__(self) -> str:
        return f"Path(start_point=({self.start_point}, width={self.start_point_width}), {', '.join([str(e) for e in self.elements])})"

    def __len__(self) -> int:
        return len(self.elements)

    def from_points(points: list[Point], fillet_radius: float = None, width: float = None, tag = None, suppress_warning: bool = False) -> Self:
        """Create a new Path using the given list of points, with the same fillet radius and width for each point"""
        path = Path(points[0], width=width)
        for point in points[1:]:
            path.append_segment(point, fillet_radius=fillet_radius, tag=tag, suppress_warning=suppress_warning)
        return path

    def from_points_with_width(points: list[tuple[Point, float]], fillet_radius: float = None, tag = None, suppress_warning: bool = False) -> Self:
        """Create a new Path using the given list of points each represented as a tuple (Point, width), with the same fillet radius for each point"""
        path = Path(points[0][0], start_point_width=points[0][1])
        for point, width in points[1:]:
            path.append_segment(point, fillet_radius=fillet_radius, width=width, tag=tag, suppress_warning=suppress_warning)
        return path

    class _Extremity(Enum):
        START = 0
        END = 1

    class _ElementType(Enum):
        SEGMENT = 0
        ARC = 1

    def _insert_element(
        self,
        extremity: _Extremity,
        element_type: _ElementType,
        to_point: Point,
        radius: float,
        anticlockwise: bool, 
        width: float,
        fillet_radius: float = None,
        fillet_end_width: float = None,
        tag = None,
        suppress_warning: bool = False,
    ) -> Self:
        """Insert a new segment or arc either at the start or at the end of this Path, with an optional fillet
        
        This method returns self and can therefore be chained."""

        # Check inputs
        if extremity not in [self._Extremity.START, self._Extremity.END]:
            raise ValueError("Invalid extremity")
        if element_type not in [self._ElementType.SEGMENT, self._ElementType.ARC]:
            raise ValueError("Invalid element type")

        # If the target point is the same as the current extremity point, it would result in a zero-length element, ignore it
        if (extremity == self._Extremity.START and to_point == self.start_point) or (extremity == self._Extremity.END and to_point == self.last_point()):
            return
        
        # For arcs, check the radius by computing an Arc object
        if element_type == self._ElementType.ARC:
            if extremity == self._Extremity.START:
                arc = Arc(to_point, self.start_point, radius=radius, suppress_warning=suppress_warning)
            elif extremity == self._Extremity.END:
                arc = Arc(self.last_point(), to_point, radius=radius, suppress_warning=suppress_warning)
        
        # Optional fillet between this element and the next element in the path
        if not self.elements:
            # Cannot add a fillet if this is the first element added to the path because there is nothing to fillet with,
            # ignore the fillet param
            fillet_radius = None
        if fillet_radius:
            # # Check the fillet radius
            # TODO : make sure this is not required
            # if element_type == self._ElementType.ARC:
            #     if fillet_radius >= radius or isclose(fillet_radius, radius):
            #         new_fillet_radius = radius - 0.1
            #         if not suppress_warning:
            #             print(f"Warning : trying to append an arc with a fillet radius greater than the arc radius, reducing the fillet radius from {fillet_radius} to {new_fillet_radius}")
            #         fillet_radius = new_fillet_radius

            # Get the three points adjacent to the current point
            if extremity == self._Extremity.START:
                p_prev, width_prev = to_point, width
                p_mid, width_mid = self.start_point, self.start_point_width
                p_next, width_next = self.elements[0].p2, self.elements[0].width
                connected_element = self.elements[0]
            elif extremity == self._Extremity.END:
                p_prev, width_prev = self.point_with_width(-2)
                p_mid, width_mid = self.point_with_width(-1)
                p_next, width_next = to_point, width
                connected_element = self.elements[-1]
            
            # Compute the fillet
            match (extremity, connected_element, element_type):
                case (_, PathSegment(), self._ElementType.SEGMENT):
                    # Append or prepend segment to segment
                    fillet = Fillet.segment_to_segment(p_prev, p_mid, p_next, fillet_radius, suppress_warning)

                case (self._Extremity.START, PathSegment(), self._ElementType.ARC):
                    # Prepend arc before segment
                    fillet = Fillet.arc_to_segment(p_prev, p_mid, p_next, radius, anticlockwise, fillet_radius, suppress_warning)

                case (self._Extremity.END, PathSegment(), self._ElementType.ARC):
                    # Append arc after segment
                    fillet = Fillet.segment_to_arc(p_prev, p_mid, p_next, radius, anticlockwise, fillet_radius, suppress_warning)

                case (self._Extremity.START, PathArc(), self._ElementType.SEGMENT):
                    # Prepend segment before arc
                    fillet = Fillet.segment_to_arc(p_prev, p_mid, p_next, connected_element.radius, connected_element.anticlockwise, fillet_radius, suppress_warning)

                case (self._Extremity.END, PathArc(), self._ElementType.SEGMENT):
                    # Append segment after arc
                    fillet = Fillet.arc_to_segment(p_prev, p_mid, p_next, connected_element.radius, connected_element.anticlockwise, fillet_radius, suppress_warning)

                case (_, PathArc(), self._ElementType.ARC):
                    # Append or prepend arc to arc
                    if not suppress_warning:
                        print("Warning : fillet between two arcs is not supported, ignoring")
                    fillet = None

            if fillet is not None:
                if fillet_end_width is not None:
                    # Custom width specified for the fillet : the start width of the fillet is the width of the point
                    # that the fillet is applied to, and the end width is the given width
                    fillet_start_width = width_mid
                else:
                    # By defaut, calculate the widths at the intermediate points of the fillet for a smooth transition
                    if isclose(width, connected_element.width):
                        fillet_start_width = fillet_end_width = width
                    else:
                        match connected_element:
                            case PathSegment():
                                prev_length_no_fillet = Segment(p_prev, p_mid).length()
                                prev_length_with_fillet = Segment(p_prev, fillet.p_arc_prev).length()
                            case PathArc():
                                prev_length_no_fillet = Arc(p_prev, p_mid, radius=connected_element.radius).length()
                                prev_length_with_fillet = Arc(p_prev, fillet.p_arc_prev, radius=connected_element.radius).length()
                        if element_type == self._ElementType.SEGMENT:
                            next_length_no_fillet = Segment(p_mid, p_next).length()
                            next_length_fillet_cut = Segment(p_mid, fillet.p_arc_next).length()
                        elif element_type == self._ElementType.ARC:
                            next_length_no_fillet = Arc(p_mid, p_next, radius=arc.radius).length()
                            next_length_fillet_cut = Arc(p_mid, fillet.p_arc_next, radius=arc.radius).length()
                        fillet_start_width = width_prev + (width_mid - width_prev) * prev_length_with_fillet / prev_length_no_fillet
                        fillet_end_width = width_mid + (width_next - width_mid) * next_length_fillet_cut / next_length_no_fillet

                # Cut the last element
                if extremity == self._Extremity.END:
                    self.elements.pop()
                    match connected_element:
                        case PathSegment():
                            connected_element_cut = PathSegment(fillet.p_arc_prev, width=fillet_start_width, tag=connected_element.tag)
                        case PathArc():
                            connected_element_cut = PathArc(fillet.p_arc_prev, connected_element.radius, connected_element.anticlockwise, width=fillet_start_width, tag=connected_element.tag)
                    self.elements.append(connected_element_cut)

                # Add the fillet as an arc
                fillet_arc_element = PathArc(fillet.p_arc_next, fillet.fillet_radius, fillet.anticlockwise_arc, width=fillet_end_width)
                if extremity == self._Extremity.START:
                    self.elements.insert(0, fillet_arc_element)
                    self.start_point = fillet.p_arc_prev
                    self.start_point_width = fillet_start_width
                elif extremity == self._Extremity.END:
                    self.elements.append(fillet_arc_element)

        # Add the element
        if extremity == self._Extremity.START:
            if element_type == self._ElementType.SEGMENT:
                element = PathSegment(self.start_point, width=self.start_point_width, tag=tag)
            elif element_type == self._ElementType.ARC:
                element = PathArc(self.start_point, arc.radius, anticlockwise, width=self.start_point_width, tag=tag)
            self.elements.insert(0, element)
            self.start_point = to_point
            self.start_point_width = width
        elif extremity == self._Extremity.END:
            if element_type == self._ElementType.SEGMENT:
                element = PathSegment(to_point, width=width, tag=tag)
            elif element_type == self._ElementType.ARC:
                element = PathArc(to_point, arc.radius, anticlockwise, width=width, tag=tag)
            self.elements.append(element)
        return self

    def append_segment(
            self,
            to_point: Point,
            fillet_radius: float = None,
            width: float = None,
            fillet_end_width: float = None,
            tag = None,
            suppress_warning: bool = False,
        ) -> Self:
        """Append a new segment at the end of this Path, with an optional fillet
        
        This method returns self and can therefore be chained."""
        return self._insert_element(
            extremity = self._Extremity.END,
            element_type = self._ElementType.SEGMENT,
            to_point = to_point,
            radius = None,
            anticlockwise = None,
            width = width or self.width,
            fillet_radius = fillet_radius,
            fillet_end_width = fillet_end_width,
            tag = tag,
            suppress_warning = suppress_warning,
        )

    def append_arc(
            self,
            to_point: Point,
            radius: float,
            anticlockwise: bool,
            fillet_radius: float = None,
            width: float = None,
            fillet_end_width: float = None,
            tag = None,
            suppress_warning: bool = False,
        ) -> Self:
        """Append a new arc of the given radius at the end of this Path, with an optional fillet
        
        This method returns self and can therefore be chained."""
        return self._insert_element(
            extremity = self._Extremity.END,
            element_type = self._ElementType.ARC,
            to_point = to_point,
            radius = radius,
            anticlockwise = anticlockwise,
            width = width or self.width,
            fillet_radius = fillet_radius,
            fillet_end_width = fillet_end_width,
            tag = tag,
            suppress_warning = suppress_warning,
        )

    def prepend_segment(
            self,
            to_point: Point,
            fillet_radius: float = None,
            width: float = None,
            fillet_end_width: float = None,
            tag = None,
            suppress_warning: bool = False,
        ) -> Self:
        """Insert a new segment at the beginning of this Path, with an optional fillet
        
        This method returns self and can therefore be chained."""
        return self._insert_element(
            extremity = self._Extremity.START,
            element_type = self._ElementType.SEGMENT,
            to_point = to_point,
            radius = None,
            anticlockwise = None,
            width = width or self.width,
            fillet_radius = fillet_radius,
            fillet_end_width = fillet_end_width,
            tag = tag,
            suppress_warning = suppress_warning,
        )

    def prepend_arc(
            self,
            to_point: Point,
            radius: float,
            anticlockwise: bool,
            fillet_radius: float = None,
            width: float = None,
            fillet_end_width: float = None,
            tag = None,
            suppress_warning: bool = False,
        ) -> Self:
        """Insert a new arc of the given radius at the beginning of this Path, with an optional fillet
        
        This method returns self and can therefore be chained."""
        return self._insert_element(
            extremity = self._Extremity.START,
            element_type = self._ElementType.ARC,
            to_point = to_point,
            radius = radius,
            anticlockwise = anticlockwise,
            width = width or self.width,
            fillet_radius = fillet_radius,
            fillet_end_width = fillet_end_width,
            tag = tag,
            suppress_warning = suppress_warning,
        )

    def append_segments(
            self,
            points: list[Point],
            fillet_radius: float = None,
            width: float = None,
            tag = None,
            suppress_warning: bool = False,
        ) -> Self:
        """Append a series of new segments at the end of this Path based on the given list of points, with an optional fillet
        
        This method returns self and can therefore be chained."""
        for point in points:
            self.append_segment(
                to_point = point,
                fillet_radius = fillet_radius,
                width = width,
                tag = tag,
                suppress_warning = suppress_warning,
            )
        return self
    
    def append_path(self, path) -> Self:
        """Append the given Path at the end of this Path after connecting them with a straight segment"""
        self.append_segment(path.start_point, width=path.start_point_width)
        for element in path.elements:
            self.elements.append(element)
        return self

    def first(self) -> PathElement:
        return self.elements[0]
    
    def element_geometry(self, index: int) -> Segment|Arc:
        """Return the element at the given index in this Path converted to a Segment or Arc"""
        try:
            if index == 0 or index == -len(self.elements):
                previous_point = self.start_point
            else:
                previous_point = self.elements[index - 1].p2
            return self.elements[index].geometry(previous_point)
        except IndexError:
            return None
    
    def elements_geometries(self) -> list[Segment|Arc]:
        """Return the elements in this Path converted to a list of individual Segments and Arcs"""
        geometries = []
        for i in range(len(self.elements)):
            geometries.append(self.element_geometry(i))
        return geometries
    
    def pop_first(self) -> PathElement:
        first = self.elements.pop(0)
        self.start_point = first.p2
        return first
    
    def first(self) -> PathElement:
        return self.elements[0]
    
    def last(self) -> PathElement:
        return self.elements[-1]
    
    def first_geometry(self) -> Segment|Arc:
        return self.elements[0].geometry(self.start_point)
    
    def last_geometry(self) -> Segment|Arc:
        if len(self.elements) == 0:
            return None
        elif len(self.elements) == 1:
            previous_point = self.start_point
        else:
            previous_point = self.elements[-2].p2
        return self.elements[-1].geometry(previous_point)
    
    def point(self, index: int) -> Point:
        if index == 0 or -index == len(self.elements) + 1:
            return self.start_point
        elif index < 0 and -index <= len(self.elements):
            return self.elements[index].p2
        elif index > 0 and index <= len(self.elements):
            return self.elements[index - 1].p2
        else:
            return None
    
    def point_with_width(self, index: int) -> Point:
        if index == 0 or -index == len(self.elements) + 1:
            return self.start_point, self.start_point_width
        elif index < 0 and -index <= len(self.elements):
            return self.elements[index].p2, self.elements[index].width
        elif index > 0 and index <= len(self.elements):
            return self.elements[index - 1].p2, self.elements[index - 1].width
        else:
            return None, None
    
    def first_point(self) -> Point:
        return self.start_point
    
    def last_point(self) -> Point:
        if self.elements:
            return self.elements[-1].p2
        else:
            return self.start_point
    
    def first_point_width(self) -> float:
        return self.start_point_width
    
    def last_point(self) -> Point:
        """Return the last Point in this Path"""
        if self.elements:
            return self.elements[-1].p2
        else:
            return self.start_point
    
    def last_point_width(self) -> float:
        if self.elements:
            return self.elements[-1].width
        else:
            return self.start_point_width
    
    def length(self) -> float:
        return sum([element.length() for element in self.elements_geometries()])
    
    def pop(self) -> PathElement:
        return self.elements.pop()
    
    def reversed(self) -> Self:
        """Return a new Path that is an copy of this Path but starting at the end"""
        path = None
        previous_element = None
        reversed_elements = list(reversed(self.elements))
        for element in reversed_elements + [None]:
            if path is None:
                if element is not None:
                    path = Path(element.p2, start_point_width=element.width, width=self.width)
                else:
                    # No element in path
                    path = Path(self.start_point, start_point_width=self.start_point_width, width=self.width)
            else:
                match previous_element:
                    case PathSegment():
                        if element is not None:
                            point = element.p2
                            width = element.width
                        else:
                            # Finish the path after the last element
                            point = self.start_point
                            width = self.start_point_width
                        path.append_segment(point, width=width, tag=previous_element.tag)
                    case PathArc():
                        if element is not None:
                            point = element.p2
                            width = element.width
                        else:
                            # Finish the path after the last element
                            point = self.start_point
                            width = self.start_point_width
                        path.append_arc(point, radius=previous_element.radius, anticlockwise=not previous_element.anticlockwise, width=width, tag=previous_element.tag)
            previous_element = element
        return path
    
    def copy(self, copy_elements: bool = True) -> Self:
        """Create a copy of this Path"""
        path = Path(
            start_point = self.start_point,
            start_point_width = self.start_point_width,
            width = self.width,
        )
        if copy_elements:
            path.elements = self.elements.copy()
        return path
    
    def translated(self, vector: Vector) -> Self:
        """Create a copy of this Path translated according to the given vector"""
        path = self.copy(copy_elements=False)
        path.start_point = path.start_point.translated(vector)
        for element in self.elements:
            path.elements.append(element.translated(vector))
        return path

    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Path rotated around the given center point by the given angle"""
        path = self.copy(copy_elements=False)
        path.start_point = path.start_point.rotated(center, angle)
        for element in self.elements:
            path.elements.append(element.rotated(center, angle))
        return path
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Path mirrored about the X axis"""
        path = self.copy(copy_elements=False)
        path.start_point = path.start_point.mirrored_x()
        for element in self.elements:
            path.elements.append(element.mirrored_x())
        return path
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Path mirrored about the Y axis"""
        path = self.copy(copy_elements=False)
        path.start_point = path.start_point.mirrored_y()
        for element in self.elements:
            path.elements.append(element.mirrored_y())
        return path
    
    def offset(self, distance: float, join_style: JoinStyle = JoinStyle.MITRE, mitre_limit=None) -> "Self":
        """Create a new Path offset from this Path separated by the given distance
        
        The distance can be negative to reverse the side of the offset."""

        # If this path is empty, we cannot calculate an offset
        if not self.elements:
            return None
        
        # If the distance is set to the special 'width' (or '-width') value, use the element's width for each point
        param_distance = distance
        if param_distance == 'width':
            start_distance = self.start_point_width / 2
            current_distance = self.elements[0].width / 2
        elif param_distance == '-width':
            start_distance = -self.start_point_width / 2
            current_distance = -self.elements[0].width / 2
        elif type(param_distance) is int or type(param_distance) is float:
            if isclose(param_distance, 0):
                # If the specified distance is zero, simply return a copy of the current path
                return self.copy()
            start_distance = param_distance
            current_distance = param_distance
        else:
            raise ValueError(f"Invalid distance value : {distance}")
        
        def _cut_offset_path(tool):
            """Cut the current offset path with some geometry element"""
            
            # Check if the tool intersects one of the last few path elements
            intersect_index = None
            for i in range(1, 7):
                element = offset_path.element_geometry(-i)
                if element is None:
                    break
                intersect = element.intersect(tool, suppress_warning=True) if element is not None else None
                if intersect is not None:
                    # If an intersection point is found, make sure it is contained on the element and the tool
                    # and not on their extrapolated basis geometries
                    if not isinstance(intersect, (list, tuple)):
                        intersect = [intersect]
                    for point in intersect:
                        if element.contains_point(point) and tool.contains_point(point):
                            intersect_index = i
                            break
            if intersect_index is not None:
                # Remove elements in the path up to the one that intersects with the tool (not included)
                for i in range(intersect_index - 1):
                    el = offset_path.pop()
            
            while len(offset_path) > 0:
                # Take the last element from the path
                last_offset_cut = offset_path.last_geometry()
                offset_path.pop()
                
                # Cut it with the tool
                last_offset_cut = last_offset_cut.cut(tool)
                
                # If the path element was completely erased with the cut, try again with
                # the next element in the path
                if last_offset_cut is None:
                    continue

                # Otherwise append the last cut element to the path again
                _append_element(last_offset_cut, cut=False)
                return

        def _append_element(offset_element, path_element=None, cut=False):
            """Add an offset element in the path up to the given point, based on the type of the previous element"""

            # Cut the offset path with the edges of the offset shape to add, in order to prevent self-intersecting paths
            if cut:
                if path_element is not None:
                    # Center line of the path
                    _cut_offset_path(path_element)
                    # End cap of the offset element
                    _cut_offset_path(Segment(offset_element.end_point(), path_element.end_point()))
                # Edge of the offset element
                _cut_offset_path(offset_element)
        
                # Cut the offset element to add with the last few elements in the offset path, for the same reasons
                offset_element_cut = offset_element
                # for i in range(5, 0, -1):
                # # for i in range(1, 3):
                #     geometry = offset_path.element_geometry(-i)
                #     if geometry is None:
                #         break
                #     offset_element_cut = offset_element_cut.cut(geometry, from_end=True)
                #     if offset_element_cut is None:
                #         # The element was completely removed by the cut
                #         break
                element_to_add = offset_element_cut
            else:
                element_to_add = offset_element

            # Add the offset element to the offset path
            if element_to_add is not None:
                match element_to_add:
                    case Segment():
                        offset_path.append_segment(
                            to_point = element_to_add.end_point(),
                        )
                    case Arc():
                        offset_path.append_arc(
                            to_point = element_to_add.end_point(),
                            radius = element_to_add.radius,
                            anticlockwise = element_to_add._reversed,
                        )
        
        # Start the offset path
        previous_element_geometry = self.elements[0].geometry(self.start_point)
        previous_offset = previous_element_geometry.offset(start_distance, current_distance)
        offset_path = Path(previous_offset.start_point())
        _append_element(previous_offset)

        for previous_element, next_element in zip(self.elements, self.elements[1:]):
            # The point that we are currently considering is at the intersection between the previous element
            # and the next element
            current_point = previous_element.p2

            # If the distance is set to the special 'width' (or '-width') value, use the element's width for each point
            if param_distance == 'width':
                current_distance = previous_element.width / 2
                next_distance = next_element.width / 2
            elif param_distance == '-width':
                current_distance = -previous_element.width / 2
                next_distance = -next_element.width / 2
            else:
                current_distance = next_distance = param_distance

            # Get the geometry object (Segment or Arc) that coincides with the next path element (PathSegment or PathArc)
            next_element_geometry = next_element.geometry(current_point)

            # Compute the offset of this geometry object
            next_offset = next_element_geometry.offset(current_distance, next_distance)

            # Check of the offset could be calculated. If not, this means the offset would be zero-length, for instance
            # the inside of a small arc with a large offset distance. In this case, ignore this element.
            if next_offset is None:
                continue
            
            # Calculate the tangeant vectors and lines on both sides of the corner
            tangeant_previous = previous_element_geometry.tangeant_at_end()
            tangeant_next = next_element_geometry.tangeant_at_start()
            line_previous = Line.from_point_and_vector(previous_offset.end_point(), tangeant_previous)
            line_next = Line.from_point_and_vector(next_offset.start_point(), tangeant_next)

            # Check the direction at which the next element starts relative to the previous element, by comparing the tangeants.
            # This value will be zero if the tangeants are colinear, negative if the next element in the path turns to the left,
            # or positive if the path turns to the right.
            # We compare this to the sign of the offset distance, which is positive for offsets on the right side and negative
            # for offsets on the left side, to determine if the offset is on the inside or the outside of the corner.
            # This will dictate how the offset path should be computed around this corner.
            tangeants_alignement = tangeant_next.cross(tangeant_previous)
            if isclose(tangeants_alignement, 0) and tangeant_next.dot(tangeant_previous) > 0:
                # The tangeants are colinear and pointing to the same direction : the end point of the previous offset should
                # be coincident to the start point of the next offset. Add the next element to the path.
                if previous_offset.end_point() != next_offset.start_point():
                    print(f"Warning : previous_offset.end_point() != next_offset.start_point(), distance={previous_offset.end_point().distance(next_offset.start_point())}")
                _append_element(next_offset, next_element_geometry)
            
            elif tangeants_alignement > 0 and current_distance > 0 or tangeants_alignement < 0 and current_distance < 0:
                # The path is on the inside of the corner : add the next offset element by cutting the current path
                _append_element(next_offset, next_element_geometry, cut=True)
            
            else:
                # The path is on the outside of the corner : the shape to add to the offset path depends on the join style
                match join_style:
                    case JoinStyle.MITRE | JoinStyle.BEVEL:
                        # Set the bevel distance depending on the join style
                        if join_style == JoinStyle.MITRE:
                            # If the mitre limit is not specified, set it to twice the distance as a reasonable value that doesn't
                            # cut right-angle corners while not extending acute corners too far
                            if mitre_limit is None or mitre_limit < 0:
                                bevel_distance = 2 * math.fabs(current_distance)
                            else:
                                bevel_distance = max(current_distance, mitre_limit)
                        else:
                            bevel_distance = math.fabs(current_distance)

                        # Calculate the bevel line and intersection points
                        bisector = tangeant_previous.reversed().bisector(tangeant_next)
                        bevel_center_point = current_point - bisector * bevel_distance
                        bevel_line = Line.from_point_and_vector(bevel_center_point, bisector.perpendicular())
                        bevel_point_1 = line_previous.intersect(bevel_line)
                        bevel_point_2 = line_next.intersect(bevel_line)

                        # Calculate the intersection between the tangeants
                        intersection_point = line_previous.intersect(line_next)
                        if intersection_point is not None and previous_offset.end_point().distance(intersection_point) < previous_offset.end_point().distance(bevel_point_1):
                            # The intersection point is closer than the bevel line : display the join as a mitre.
                            # Add segments to connect to the intersection point then the start of the next offset.
                            _append_element(Segment(previous_offset.end_point(), intersection_point))
                            _append_element(Segment(intersection_point, next_offset.start_point()))
                        else:
                            # The intersection point is further than the bevel line : display the join as a bevel.
                            # Add segments to connect to the bevel points then to the start of the next offset
                            _append_element(Segment(previous_offset.end_point(), bevel_point_1))
                            _append_element(Segment(bevel_point_1, bevel_point_2))
                            _append_element(Segment(bevel_point_2, next_offset.start_point()))
                        
                        # Add the next offset element
                        _append_element(next_offset, next_element_geometry, cut=True)

                    case JoinStyle.ROUND:
                        # Add an arc between the offset paths
                        if previous_offset.end_point() != next_offset.start_point():
                            _append_element(Arc(previous_offset.end_point(), next_offset.start_point(), radius=math.fabs(current_distance), reverse=current_distance > 0))

                        # Add the next offset element
                        _append_element(next_offset, next_element_geometry, cut=True)

            # Advance one step in the path
            previous_element_geometry = next_element_geometry
            previous_offset = next_offset

        return offset_path

    def stroke(self, cap_style: CapStyle = CapStyle.FLAT, join_style = JoinStyle.MITRE, mitre_limit = None, arcs_discretization_length: float = 1, arcs_discretization_angle: float = 1) -> "Polygon":
        """Return a Polygon that represents the outline of the stroke of this Path"""
        if not self.elements:
            # Path is empty
            return None

        path_left = self.offset('-width', join_style=join_style, mitre_limit=mitre_limit)
        path_right = self.offset('width', join_style=join_style, mitre_limit=mitre_limit)
        match cap_style:
            case CapStyle.FLAT:
                # Append a straight segment at the start and end of the left offset path to connect respectively
                # to the start and end of the right offset path
                path_left.prepend_segment(path_right.first_point())
                path_left.append_segment(path_right.last_point())
            case CapStyle.SQUARE:
                # Append segments to extend the offset paths the same amount as the width of each end point
                for path in [path_left, path_right]:
                    # Start
                    geometry = path.last_geometry()
                    path.append_segment(geometry.end_point() + geometry.tangeant_at_end() * self.last_point_width()/2)

                    # End
                    geometry = path.first_geometry()
                    path.prepend_segment(geometry.start_point() - geometry.tangeant_at_start() * self.first_point_width()/2)
                
                # Finish by capping with straight segments
                path_left.append_segment(path_right.last_point())
                path_left.prepend_segment(path_right.first_point())
            case CapStyle.ROUND:
                # Append 180° arcs at the start and end of the left offset path to connect respectively to the
                # start and end of the right offset path
                path_left.prepend_arc(path_right.first_point(), radius=self.first_point_width()/2, anticlockwise=False)
                path_left.append_arc(path_right.last_point(), radius=self.last_point_width()/2, anticlockwise=False)
        polygonPath = path_left.append_path(path_right.reversed())
        return Polygon(polygonPath, arcs_discretization_length=arcs_discretization_length, arcs_discretization_angle=arcs_discretization_angle)
    
    def _export(self):
        # Debug function that returns the Python code that would create the same Path
        output = ""
        output += f"path = Path(Point({self.start_point.x}, {self.start_point.y}), width={self.start_point_width})\n"
        for element in self.elements:
            match element:
                case PathSegment():
                    output += f"path.append_segment(Point({element.p2.x}, {element.p2.y}), width={element.width})\n"
                case PathArc():
                    output += f"path.append_arc(Point({element.p2.x}, {element.p2.y}), radius={element.radius}, anticlockwise={element.anticlockwise}, width={element.width})\n"
        return output

    def to_svg(self, closed: bool = False) -> str:
        d = ""
        d += f"M{self.start_point.to_svg()} "
        current_point = self.start_point
        for element in self.elements:
            match element:
                case PathSegment():
                    d += f"L{element.p2.to_svg()} "
                    current_point = element.p2
                case PathArc():
                    d += f"A{round(element.radius, PRECISION)},{round(element.radius, PRECISION)} 0 0,{int(not element.anticlockwise)} {element.p2.to_svg()} "
                    current_point = element.p2
        if closed:
            d += "Z"
        return d
    
    def draw_svg(
            self,
            drawing: svg.Drawing,
            parent: svg.base.BaseElement = None,
            color: str = None,
            opacity: float = None,
            line_width: float = None,
            dashes: str = None,
            label: str = None,
        ) -> Self:
        """Draw this Path on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if parent is None:
            parent = drawing

        parent.add(drawing.path(
            d = self.to_svg(),
            stroke = color or style.path_color,
            stroke_opacity = opacity or style.path_opacity,
            stroke_width = line_width or style.path_line_width,
            stroke_dasharray = dashes or style.path_dashes,
            label = label,
        ))
        return self
    
    def draw_kicad(self, kicadpcb: "KicadPCB", width: float, layer: str, stroke_type: str = 'solid') -> Self:
        """No-op : construction Paths are not printed on the Kicad board"""
        return self

class Polygon(DrawableObject):
    """A Polygon represented by a closed Path"""

    def __init__(self, path: Path, arcs_discretization_length: float, arcs_discretization_angle: float):
        self.path = path
        self.discretized_points = None
        if arcs_discretization_length and arcs_discretization_angle:
            points = [self.path.start_point]
            for element in self.path.elements_geometries():
                match element:
                    case Segment():
                        points.append(element.end_point())
                    case Arc():
                        points.extend(element.discretize(segment_length=arcs_discretization_length, max_deviation=arcs_discretization_angle)[1:])
            self.discretized_points = points
    
    def copy(self) -> Self:
        """Create a copy of this Polygon"""
        polygon = Polygon(
            path = self.path,
            arcs_discretization_length = None,
            arcs_discretization_angle = None,
        )
        polygon.discretized_points = self.discretized_points
        return polygon
    
    def translated(self, vector: Vector) -> Self:
        """Create a copy of this Polygon translated according to the given vector"""
        polygon = self.copy()
        polygon.path = polygon.path.translated(vector)
        points = []
        for point in polygon.discretized_points:
            points.append(point.translated(vector))
        polygon.discretized_points = points
        return polygon

    def rotated(self, center: Self, angle: float) -> Self:
        """Create a copy of this Polygon rotated around the given center point by the given angle"""
        polygon = self.copy()
        polygon.path = polygon.path.rotated(center, angle)
        points = []
        for point in polygon.discretized_points:
            points.append(point.rotated(center, angle))
        polygon.discretized_points = points
        return polygon
    
    def mirrored_x(self) -> Self:
        """Create a copy of this Polygon mirrored about the X axis"""
        polygon = self.copy()
        polygon.path = polygon.path.mirrored_x()
        points = []
        for point in polygon.discretized_points:
            points.append(point.mirrored_x())
        polygon.discretized_points = points
        return polygon
    
    def mirrored_y(self) -> Self:
        """Create a copy of this Polygon mirrored about the Y axis"""
        polygon = self.copy()
        polygon.path = polygon.path.mirrored_y()
        points = []
        for point in polygon.discretized_points:
            points.append(point.mirrored_y())
        polygon.discretized_points = points
        return polygon
    
    def draw_svg(
            self,
            drawing: svg.Drawing,
            parent: svg.base.BaseElement = None,
            color: str = None,
            opacity: float = None,
            stroke_color: str = None,
            stroke_width: float = None,
            stroke_opacity: float = None,
            stroke_dashes: str = None,
            label: str = None,
        ) -> Self:
        """Draw this Polygon on the given SVG drawing
        
        This method returns self and can therefore be chained."""

        if parent is None:
            parent = drawing

        parent.add(drawing.path(
            d = self.path.to_svg(closed=True),
            fill = color or style.polygon_color,
            opacity = opacity or style.polygon_opacity,
            stroke = stroke_color or style.polygon_stroke_color,
            stroke_width = stroke_width or style.polygon_stroke_width,
            stroke_opacity = stroke_opacity or style.polygon_stroke_opacity,
            stroke_dasharray = stroke_dashes or style.polygon_stroke_dashes,
            label = label,
        ))
        return self
    
    def draw_kicad(self, kicadpcb: "KicadPCB", layer: str) -> Self:
        """Draw this Polygon on the given Kicad board"""
        kicadpcb.gr_poly(
            points = self.discretized_points,
            layer = layer,
            fill = True,
        )
        return self