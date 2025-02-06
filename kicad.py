from sexpdata import dump, Symbol as S
from datetime import datetime
import uuid
from geometry import Vector, Point, Segment, Arc, Path, PathSegment, PathArc
from version import PCBMOT_VERSION

class KicadPCB:
    def __init__(self, filename, config):
        self.filename = filename
        self.config = config

        self.board_origin = Point(150, 100)

        # Mapping between internal layer names and Kicad layers
        self.layers_def = {
            'top': (0, "F.Cu", S('signal')),
            'in1': (1, "In1.Cu", S('signal')),
            'in2': (2, "In2.Cu", S('signal')),
            'in3': (3, "In3.Cu", S('signal')),
            'in4': (4, "In4.Cu", S('signal')),
            'in5': (5, "In5.Cu", S('signal')),
            'in6': (6, "In6.Cu", S('signal')),
            'bottom': (31, "B.Cu", S('signal')),
            'bottom_silk': (36, "B.SilkS", S('user'), "B.Silkscreen"),
            'top_silk': (37, "F.SilkS", S('user'), "F.Silkscreen"),
            'bottom_mask': (38, "B.Mask", S('user')),
            'top_mask': (39, "F.Mask", S('user')),
            'construction': (50, "User.1", S('user')),
            'magnets': (51, "User.2", S('user')),
            'outline': (44, "Edge.Cuts", S('user')),
        }
        self.layers = config.copper_layers + ['bottom_silk', 'top_silk', 'bottom_mask', 'top_mask', 'construction', 'magnets', 'outline']

        self.data = \
            [S('kicad_pcb'),
                [S('version'), S("20240108")],
                [S('generator'), 'PCBMot'],
                [S('generator_version'), str(PCBMOT_VERSION)],
                [S('general'),
                    [S('thickness'), config.board_thickness],
                ],
                [S('paper'), 'A4'],
                [S('title_block'),
                    [S('title'), "PCB Motor"],
                    [S('date'), datetime.today().strftime('%d/%m/%Y')],
                    [S('comment'), 1, "Generated using PCBMot"],
                    [S('comment'), 1, "https://pcbmot.silica.io"],
                ],
                [S('layers'),
                    *[self.layers_def[layer_id] for layer_id in self.layers],
                ],
                [S('setup'),
                    [S('pad_to_mask_clearance'), 0],
                    [S('aux_axis_origin'), self.board_origin.x, self.board_origin.y],
                    [S('grid_origin'), self.board_origin.x, self.board_origin.y],
                    [S('pcbplotparams')],
                ],
            ]
    
    def map_coords(self, point: Point) -> tuple[float]:
        return (self.board_origin.x + point.x, self.board_origin.y - point.y)

    def segment(self, p1: Point, p2: Point, width: float, layer: str):
        """Add a track segment to the board"""
        self.data.append(
            [S('segment'),
                [S('start'), *self.map_coords(p1)],
                [S('end'), *self.map_coords(p2)],
                [S('width'), width],
                [S('layer'), self.layers_def[layer][1]],
                [S('uuid'), str(uuid.uuid4())],
            ]
        )
    
    def arc(self, p1: Point, p2: Point, midpoint: Point, width: float, layer: str):
        """Add a track arc to the board"""
        self.data.append(
            [S('arc'),
                [S('start'), *self.map_coords(p1)],
                [S('mid'), *self.map_coords(midpoint)],
                [S('end'), *self.map_coords(p2)],
                [S('width'), width],
                [S('layer'), self.layers_def[layer][1]],
                [S('uuid'), str(uuid.uuid4())],
            ]
        )
    
    def path(self, path: Path, width: float, layer: str):
        """Add a track path to the board"""
        p1 = path.start_point
        for element in path.elements:
            match element:
                case PathSegment():
                    self.segment(p1, element.p2, width, layer)
                    p1 = element.p2
                case PathArc():
                    arc = Arc(p1, element.p2, element.radius, element.anticlockwise)
                    self.arc(arc.p1, arc.p2, arc.midpoint(), width, layer)
                    p1 = element.p2
    
    def via(self, center: Point, diameter: float, hole_diameter: float):
        """Add a via to the board"""
        self.data.append(
            [S('via'),
                [S('at'), *self.map_coords(center)],
                [S('size'), diameter],
                [S('drill'), hole_diameter],
                [S('layers'), self.layers_def[self.config.copper_layers[0]][1], self.layers_def[self.config.copper_layers[-1]][1]],
                [S('uuid'), str(uuid.uuid4())],
            ]
        )
    
    def terminal(
            self,
            center: Point,
            diameter: float,
            hole_diameter: float,
            ref: str,
            ref_offset: Vector,
            ref_angle: float,
            ref_size: float,
        ):
        """Add a terminal to the board"""
        self.data.append(
            [S('footprint'), f"Terminal:Terminal_THTPad_D{round(diameter, 1)}mm_Drill{round(hole_diameter, 1)}mm",
                [S('at'), *self.map_coords(center)],
                [S('layer'), self.layers_def[self.config.copper_layers[0]][1]],
                [S('tags'), "terminal"],
                [S('uuid'), str(uuid.uuid4())],
                [S('attr'), S('exclude_from_pos_files'), S('exclude_from_bom')],
                [S('property'), "Reference", f"{ref}",
                    [S('at'), ref_offset.x, -ref_offset.y, -ref_angle],
                    [S('unlocked'), S('yes')],
                    [S('layer'), self.layers_def['top_silk'][1]],
                    [S('uuid'), str(uuid.uuid4())],
                    [S('effects'),
                        [S('font'),
                            [S('size'), ref_size, ref_size],
                            [S('thickness'), ref_size * 0.2],
                        ],
                    ],
                ],
                [S('fp_circle'),
                    [S('center'), 0, 0],
                    [S('end'), 0, (diameter / 2.0) * 1.1],
                    [S('stroke'),
                        [S('width'), (diameter / 2.0) * 0.1],
                        [S('type'), S('solid')],
                    ],
                    [S('fill'), S('none')],
                    [S('layer'), self.layers_def['top_silk'][1]],
                    [S('uuid'), str(uuid.uuid4())],
                ],
                [S('fp_circle'),
                    [S('center'), 0, 0],
                    [S('end'), 0, (diameter / 2.0) * 1.1],
                    [S('stroke'),
                        [S('width'), (diameter / 2.0) * 0.1],
                        [S('type'), S('solid')],
                    ],
                    [S('fill'), S('none')],
                    [S('layer'), self.layers_def['bottom_silk'][1]],
                    [S('uuid'), str(uuid.uuid4())],
                ],
                [S('pad'), "1", S('thru_hole'), S('circle'),
                    [S('at'), 0, 0],
                    [S('size'), diameter, diameter],
                    [S('drill'), hole_diameter, hole_diameter],
                    [S('layers'), "*.Cu", "*.Mask"],
                    [S('remove_unused_layers'), S('no')],
                    [S('uuid'), str(uuid.uuid4())],
                ],
            ]
        )
    
    def gr_line(self, p1: Point, p2: Point, width: float, layer: str, stroke_type: str = 'solid'):
        """Add a graphical line to the board"""
        self.data.append(
            [S('gr_line'),
                [S('start'), *self.map_coords(p1)],
                [S('end'), *self.map_coords(p2)],
                [S('layer'), self.layers_def[layer][1]],
                [S('stroke'),
                    [S('width'), width],
                    [S('type'), S(stroke_type)],
                ],
                [S('uuid'), str(uuid.uuid4())],
            ]
        )
    
    def gr_circle(self, center: Point, radius: float, width: float, layer: str, stroke_type: str = 'solid'):
        """Add a graphical circle to the board"""
        self.data.append(
            [S('gr_circle'),
                [S('center'), *self.map_coords(center)],
                [S('end'), *self.map_coords(center + Vector(radius, 0.0))],
                [S('layer'), self.layers_def[layer][1]],
                [S('stroke'),
                    [S('width'), width],
                    [S('type'), S(stroke_type)],
                ],
                [S('uuid'), str(uuid.uuid4())],
            ]
        )
    
    def gr_arc(self, p1: Point, p2: Point, midpoint: Point, width: float, layer: str, stroke_type: str = 'solid'):
        """Add a graphical arc to the board"""
        self.data.append(
            [S('gr_arc'),
                [S('start'), *self.map_coords(p1)],
                [S('mid'), *self.map_coords(midpoint)],
                [S('end'), *self.map_coords(p2)],
                [S('stroke'),
                    [S('width'), width],
                    [S('type'), S(stroke_type)],
                ],
                [S('layer'), self.layers_def[layer][1]],
                [S('uuid'), str(uuid.uuid4())],
            ]
        )
    
    def gr_path(self, path: Path, width: float, layer: str, stroke_type: str = 'solid'):
        """Add a graphical path to the board"""
        p1 = path.start_point
        for element in path.elements:
            match element:
                case PathSegment():
                    self.gr_line(p1, element.p2, width, layer, stroke_type)
                    p1 = element.p2
                case PathArc():
                    arc = Arc(p1, element.p2, element.radius, element.anticlockwise)
                    self.gr_arc(arc.p1, arc.p2, arc.midpoint(), width, layer, stroke_type)
                    p1 = element.p2
    
    def text(self, text: str, center: Point, angle: float, font_size: float, layer: str):
        """Add a text to the board"""
        self.data.append(
            [S('gr_text'), text,
                [S('at'), *self.map_coords(center), -angle],
                [S('layer'), self.layers_def[layer][1]],
                [S('effects'),
                    [S('font'),
                        [S('size'), font_size, font_size],
                        [S('thickness'), font_size * 0.2],
                    ],
                ],
                [S('uuid'), str(uuid.uuid4())],
            ]
        )
    
    def save(self):
        """Save the board into a file"""
        with open(self.filename, 'w') as file:
            return dump(self.data, file, pretty_print=True)
