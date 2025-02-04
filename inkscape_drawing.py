# Based on 
#   https://github.com/mozman/svgwrite/blob/master/examples/inkscape_drawing.py
# and
#   https://github.com/mozman/svgwrite/blob/master/svgwrite/elementfactory.py#L71
# Original authors : Antonio Ospite <ao2@ao2.it>, Manfred Moitzi <mozman@gmx.at>

import svgwrite as svg
from svgwrite.data.types import SVGAttribute

class InkscapeDrawing(svg.Drawing):
    """An svgwrite.Drawing subclass which supports Inkscape layers"""

    INKSCAPE_NAMESPACE = 'http://www.inkscape.org/namespaces/inkscape'
    SODIPODI_NAMESPACE = 'http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd'

    # SVG objects that will customized to support Inkscape properties
    object_types = {
        'g': svg.container.Group,
        'rect': svg.shapes.Rect,
        'circle': svg.shapes.Circle,
        'path': svg.path.Path,
        'text': svg.text.Text,
    }

    def __init__(self, *args, **extra):
        # Set some default values to the drawing
        extra.setdefault('fill', 'none')
        extra.setdefault('stroke', '#000000')
        extra.setdefault('stroke_width', 0.5)
        extra.setdefault('stroke_linejoin', 'round')
        extra.setdefault('stroke_linecap', 'round')
        background_color = extra.pop('background_color', None)
        coords = extra.pop('coords', None)
        if coords is not None:
            extra['viewBox'] = f"{coords[0]} {coords[1]} {coords[2]} {coords[3]}"

        # Call the parent constructor
        super(InkscapeDrawing, self).__init__(*args, **extra)

        # Add the Inkscape-specific attributes to the validator
        inkscape_attributes = {
            'xmlns:inkscape': SVGAttribute(
                'xmlns:inkscape',
                anim=False,
                types=[],
                const=frozenset([self.INKSCAPE_NAMESPACE])
            ),
            'xmlns:sodipodi': SVGAttribute(
                'xmlns:sodipodi',
                anim=False,
                types=[],
                const=frozenset([self.SODIPODI_NAMESPACE])
            ),
            'inkscape:groupmode': SVGAttribute(
                'inkscape:groupmode',
                anim=False,
                types=[],
                const=frozenset(['layer'])
            ),
            'inkscape:label': SVGAttribute(
                'inkscape:label',
                anim=False,
                types=frozenset(['string']),
                const=[]
            ),
            'sodipodi:insensitive': SVGAttribute(
                'sodipodi:insensitive',
                anim=False,
                types=frozenset(['string']),
                const=[]
            ),
        }
        self.validator.attributes.update(inkscape_attributes)
        elements = self.validator.elements

        # Add the attributes to the svg base element
        svg_attributes = set(elements['svg'].valid_attributes)
        svg_attributes.add('xmlns:inkscape')
        svg_attributes.add('xmlns:sodipodi')
        elements['svg'].valid_attributes = frozenset(svg_attributes)

        # Add the attributes to the graphical object types
        for object_type in self.object_types:
            g_attributes = set(elements[object_type].valid_attributes)
            g_attributes.add('inkscape:label')
            g_attributes.add('sodipodi:insensitive')
            if object_type == 'g':
                g_attributes.add('inkscape:groupmode')
            elements[object_type].valid_attributes = frozenset(g_attributes)

        # Include the Inkscape namespaces
        self['xmlns:inkscape'] = self.INKSCAPE_NAMESPACE
        self['xmlns:sodipodi'] = self.SODIPODI_NAMESPACE

        # Create the background
        if background_color is not None:
            self.add(self.rect(
                (coords[0], coords[1]),
                (coords[2], coords[3]),
                stroke = 'none',
                fill = background_color,
                label = 'Background',
                locked = True,
            ))

    def layer(self, **extra):
        """Create an Inkscape layer with an optional `label` and `locked` flag"""

        # Create a group and set the groupmode attribute so Inkscape recognises it as a layer
        layer_group = self.g(**extra)
        layer_group['inkscape:groupmode'] = 'layer'

        return layer_group
    
    class ObjectBuilder():
        """A custom object builder that supports the `label` and `locked` parameters"""
        def __init__(self, cls, factory):
            self.cls = cls
            self.factory = factory
        
        def __call__(self, *args, **kwargs):
            """Create an SVG object with an optional `label` and `locked` flag"""
            label = kwargs.pop('label', None)
            locked = kwargs.pop('locked', None)
            kwargs['factory'] = self.factory
            object = self.cls(*args, **kwargs)
            if label is not None:
                object['inkscape:label'] = label
            if locked is not None:
                object['sodipodi:insensitive'] = '1' if locked else '0'
            return object

    def __getattr__(self, name):
        if name in self.object_types:
            return self.ObjectBuilder(self.object_types[name], self)
        else:
            raise AttributeError(f"'%s' has no attribute '%s'" % (self.__class__.__name__, name))
