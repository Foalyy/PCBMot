from enum import Enum
from typing import Self
import math, json

class BoardShape(Enum):
    CIRCLE = 1
    SQUARE = 2

class TerminalType(Enum):
    NONE = 1
    THROUGH_HOLE = 2
    SMD = 3
    CASTELLATED = 4

## Board config
class Config:
    options = [
        # General
        {
            'name': 'temperature',
            'json': 'general.temperature',
            'type': float,
            'default': 20.0, # °C
        },
        {
            'name': 'copper_resistivity',
            'json': 'general.copper_resistivity',
            'type': float,
            'default': 1.68e-5, # Ohm.mm at 20°C
        },
        {
            'name': 'copper_temperature_coefficient',
            'json': 'general.copper_temperature_coefficient',
            'type': float,
            'default': 4.04e-3, # K^-1
        },

        # Board
        {
            'name': 'board_shape',
            'json': 'board.board_shape',
            'type': BoardShape,
            'decoder': lambda s : {
                'circle': BoardShape.CIRCLE,
                'square': BoardShape.SQUARE,
            }[s],
            'encoder': lambda s : {
                BoardShape.CIRCLE: 'circle',
                BoardShape.SQUARE: 'square',
            }[s],
            'enum': ['circle', 'square'],
            'default': 'circle',
        },
        {
            'name': 'board_diameter',
            'json': 'board.board_diameter',
            'type': float,
            'required': True,
        },
        {
            'name': 'hole_diameter',
            'json': 'board.hole_diameter',
            'type': float,
            'required': True,
        },
        {
            'name': 'board_chamfer',
            'json': 'board.board_chamfer',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'board_fillet',
            'json': 'board.board_fillet',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'board_thickness',
            'json': 'board.board_thickness',
            'type': float,
            'required': True,
        },
        {
            'name': 'n_layers',
            'json': 'board.n_layers',
            'type': int,
            'enum': [2, 4, 6, 8],
            'required': True,
        },
        {
            'name': 'board_outer_margin',
            'json': 'board.board_outer_margin',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'board_inner_margin',
            'json': 'board.board_inner_margin',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'rotation',
            'json': 'board.rotation',
            'type': float,
            'default': 0.0,
        },
        {
            'name': 'outer_layers_copper_thickness',
            'json': 'board.outer_layers_copper_thickness',
            'type': float,
            'default': 0.0348, # 1 oz
        },
        {
            'name': 'inner_layers_copper_thickness',
            'json': 'board.inner_layers_copper_thickness',
            'type': float,
            'default': 0.0174, # 0.5 oz
        },

        # Board / mountpoints
        {
            'name': 'mountpoints',
            'json': 'board.mountpoints.mountpoints',
            'type': bool,
            'default': False,
        },
        {
            'name': 'n_mountpoints',
            'json': 'board.mountpoints.n_mountpoints',
            'type': int,
            'default': 0,
        },
        {
            'name': 'mountpoints_position_radius',
            'json': 'board.mountpoints.mountpoints_position_radius',
            'type': float,
            'default': 0.0,
        },
        {
            'name': 'mountpoints_diameter',
            'json': 'board.mountpoints.mountpoints_diameter',
            'type': float,
            'default': 0.0,
        },
        {
            'name': 'mountpoints_marking_diameter',
            'json': 'board.mountpoints.mountpoints_marking_diameter',
            'type': float,
            'default': 0.0,
        },

        # Coils
        {
            'name': 'trace_width',
            'json': 'coils.trace_width',
            'type': float,
            'required': True,
        },
        {
            'name': 'trace_spacing',
            'json': 'coils.trace_spacing',
            'type': float,
            'required': True,
        },
        {
            'name': 'n_phases',
            'json': 'coils.n_phases',
            'type': int,
            'required': True,
        },
        {
            'name': 'n_slots_per_phase',
            'json': 'coils.n_slots_per_phase',
            'type': int,
            'required': True,
        },
        {
            'name': 'coil_angle',
            'json': 'coils.coil_angle',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'four_layers_inner_vias',
            'json': 'coils.four_layers_inner_vias',
            'type': bool,
            'default': False,
        },
        {
            'name': 'max_turns_per_layer',
            'json': 'coils.max_turns_per_layer',
            'type': int,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'coil_names_font_size',
            'json': 'coils.coil_names_font_size',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'coil_names_offset',
            'json': 'coils.coil_names_offset',
            'type': float,
            'nullable': True,
            'default': None,
        },

        # Coils / series link
        {
            'name': 'link_series_coils',
            'json': 'coils.series_link.link_series_coils',
            'type': bool,
            'default': True,
        },
        {
            'name': 'series_link_inner_trace_width',
            'json': 'coils.series_link.series_link_inner_trace_width',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'series_link_outer_trace_width',
            'json': 'coils.series_link.series_link_outer_trace_width',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'series_link_inner_offset',
            'json': 'coils.series_link.series_link_inner_offset',
            'type': float,
            'default': 0.0,
        },
        {
            'name': 'series_link_outer_offset',
            'json': 'coils.series_link.series_link_outer_offset',
            'type': float,
            'default': 0.0,
        },

        # Coils / COM link
        {
            'name': 'link_com',
            'json': 'coils.com_link.link_com',
            'type': bool,
            'default': True,
        },
        {
            'name': 'com_link_trace_width',
            'json': 'coils.com_link.com_link_trace_width',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'com_link_offset',
            'json': 'coils.com_link.com_link_offset',
            'type': float,
            'default': 0.0,
        },

        # Vias
        {
            'name': 'via_diameter',
            'json': 'vias.via_diameter',
            'type': float,
            'required': True,
        },
        {
            'name': 'via_drill_diameter',
            'json': 'vias.via_drill_diameter',
            'type': float,
            'required': True,
        },
        {
            'name': 'outer_vias_offset',
            'json': 'vias.outer_vias_offset',
            'type': float,
            'default': 0.0,
        },
        {
            'name': 'inner_vias_offset',
            'json': 'vias.inner_vias_offset',
            'type': float,
            'default': 0.0,
        },
        {
            'name': 'via_resistance',
            'json': 'vias.via_resistance',
            'type': float,
            'default': 0.0, # Ohm
        },

        # Terminals
        {
            'name': 'terminal_type',
            'json': 'terminals.terminal_type',
            'type': TerminalType,
            'decoder': lambda s : {
                'none': TerminalType.NONE,
                'through_hole': TerminalType.THROUGH_HOLE,
                'smd': TerminalType.SMD,
                'castellated': TerminalType.CASTELLATED,
            }[s],
            'encoder': lambda s : {
                TerminalType.NONE: 'none',
                TerminalType.THROUGH_HOLE: 'through_hole',
                TerminalType.SMD: 'smd',
                TerminalType.CASTELLATED: 'castellated',
            }[s],
            'enum': ['none', 'through_hole', 'smd', 'castellated'],
            'default': 'through_hole',
        },
        {
            'name': 'terminal_diameter',
            'json': 'terminals.terminal_diameter',
            'type': float,
            'default': 1.6,
        },
        {
            'name': 'terminal_hole_diameter',
            'json': 'terminals.terminal_hole_diameter',
            'type': float,
            'default': 0.8,
        },
        {
            'name': 'terminal_offset',
            'json': 'terminals.terminal_offset',
            'type': float,
            'default': 0.0,
        },

        # Magnets
        {
            'name': 'magnets_diameter',
            'json': 'magnets.magnets_diameter',
            'type': float,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'magnets_position_radius',
            'json': 'magnets.magnets_position_radius',
            'type': float,
            'nullable': True,
            'default': None,
        },

        # Display
        {
            'name': 'draw_vias',
            'json': 'display.draw_vias',
            'type': bool,
            'default': True,
        },
        {
            'name': 'draw_terminals',
            'json': 'display.draw_terminals',
            'type': bool,
            'default': True,
        },
        {
            'name': 'draw_outline',
            'json': 'display.draw_outline',
            'type': bool,
            'default': True,
        },
        {
            'name': 'draw_construction_geometry',
            'json': 'display.draw_construction_geometry',
            'type': bool,
            'default': True,
        },
        {
            'name': 'draw_only_layers',
            'json': 'display.draw_only_layers',
            'type': list,
            'nullable': True,
            'default': None,
        },
        {
            'name': 'draw_coil_names',
            'json': 'display.draw_coil_names',
            'type': bool,
            'default': True,
        },
        {
            'name': 'draw_magnets',
            'json': 'display.draw_magnets',
            'type': bool,
            'default': True,
        },

        # SVG
        {
            'name': 'svg_profile',
            'json': 'svg.svg_profile',
            'type': str,
            'enum': ['tiny', 'full'],
            'default': 'tiny',
        },
        {
            'name': 'svg_scale',
            'json': 'svg.svg_scale',
            'type': float,
            'default': 100,
        },
    ]

    def __init__(
        self,
        temperature: float,
        copper_resistivity: float,
        copper_temperature_coefficient: float,
        board_shape: BoardShape,
        board_diameter: float,
        hole_diameter: float,
        board_chamfer: float,
        board_fillet: float,
        mountpoints: bool,
        n_mountpoints: int,
        mountpoints_position_radius: float,
        mountpoints_diameter: float,
        mountpoints_marking_diameter: float,
        board_thickness: float,
        board_outer_margin: float,
        board_inner_margin: float,
        rotation: float,
        inner_layers_copper_thickness: float,
        outer_layers_copper_thickness: float,
        n_phases: int,
        four_layers_inner_vias: bool,
        n_slots_per_phase: int,
        coil_angle: float,
        n_layers: int,
        max_turns_per_layer: int,
        trace_width: float,
        trace_spacing: float,
        via_diameter: float,
        via_drill_diameter: float,
        outer_vias_offset: float,
        inner_vias_offset: float,
        via_resistance: float,
        terminal_type: TerminalType,
        terminal_diameter: float,
        terminal_hole_diameter: float,
        terminal_offset: float,
        link_series_coils: bool,
        series_link_inner_trace_width: float,
        series_link_outer_trace_width: float,
        series_link_inner_offset: float,
        series_link_outer_offset: float,
        link_com: bool,
        com_link_trace_width: float,
        com_link_offset: float,
        coil_names_font_size: float,
        coil_names_offset: float,
        magnets_diameter: float,
        magnets_position_radius: float,
        draw_vias: bool,
        draw_terminals: bool,
        draw_outline: bool,
        draw_construction_geometry: bool,
        draw_only_layers: list[str],
        draw_magnets: bool,
        draw_coil_names: bool,
        svg_profile: str,
        svg_scale: int,
    ):
        # Check parameters
        if board_diameter <= 0:
            raise ValueError("The board diameter must be positive")
        if hole_diameter <= 0:
            raise ValueError("The hole diameter must be positive")
        if hole_diameter >= board_diameter:
            raise ValueError("The board diameter is smaller than the hole diameter")
        if not (n_slots_per_phase >= 1 and (n_slots_per_phase == 1 or n_slots_per_phase % 2 == 0)):
            raise ValueError("The number of slots per phase must be 1 or an even number")
        if link_series_coils and n_phases > n_layers:
            raise ValueError("Unable to generate the connections between the coils because the number of layers should be greater than or equal to the number of phases, either increase the number of layers or disable link_series_coils")
        if terminal_type not in [TerminalType.NONE, TerminalType.THROUGH_HOLE, TerminalType.SMD, TerminalType.CASTELLATED]:
            raise ValueError("Invalid terminal type")
        if terminal_type != TerminalType.NONE and terminal_diameter <= 0:
            raise ValueError("Please provide a pad diameter for the terminals")
        if terminal_type != TerminalType.NONE and terminal_hole_diameter <= 0:
            raise ValueError("Please provide a hole diameter for the terminals")
        if terminal_type == TerminalType.CASTELLATED and board_shape != BoardShape.CIRCLE:
            raise ValueError("Castellated terminals are only possible for circle boards")
        if n_mountpoints > 0 and mountpoints_position_radius <= 0:
            raise ValueError("Please provide a position radius for the mountpoints")
        if n_mountpoints > 0 and mountpoints_diameter <= 0:
            raise ValueError("Please provide a diameter for the mountpoints")


        # Save the parameters
        self.temperature: float = temperature
        self.copper_resistivity: float = copper_resistivity
        self.copper_temperature_coefficient: float = copper_temperature_coefficient
        self.board_shape: BoardShape = board_shape
        self.board_diameter: float = board_diameter
        self.hole_diameter: float = hole_diameter
        self.board_chamfer: float = board_chamfer
        self.board_fillet: float = board_fillet
        self.mountpoints: bool = mountpoints
        self.n_mountpoints: int = n_mountpoints
        self.mountpoints_position_radius: float = mountpoints_position_radius
        self.mountpoints_diameter: float = mountpoints_diameter
        self.mountpoints_marking_diameter: float = mountpoints_marking_diameter
        self.board_thickness: float = board_thickness
        self.board_outer_margin: float = board_outer_margin
        self.board_inner_margin: float = board_inner_margin
        self.rotation: float = rotation
        self.inner_layers_copper_thickness: float = inner_layers_copper_thickness
        self.outer_layers_copper_thickness: float = outer_layers_copper_thickness
        self.n_phases: int = n_phases
        self.n_slots_per_phase: int = n_slots_per_phase
        self.coil_angle: float = coil_angle
        self.n_layers: int = n_layers
        self.four_layers_inner_vias: bool = four_layers_inner_vias
        self.max_turns_per_layer: int = max_turns_per_layer
        self.trace_width: float = trace_width
        self.trace_spacing: float = trace_spacing
        self.via_diameter: float = via_diameter
        self.via_drill_diameter: float = via_drill_diameter
        self.outer_vias_offset: float = outer_vias_offset
        self.inner_vias_offset: float = inner_vias_offset
        self.via_resistance: float = via_resistance
        self.terminal_type: TerminalType = terminal_type
        self.terminal_diameter: float = terminal_diameter
        self.terminal_hole_diameter: float = terminal_hole_diameter
        self.terminal_offset: float = terminal_offset
        self.link_series_coils: bool = link_series_coils
        self.series_link_inner_trace_width: float = series_link_inner_trace_width
        self.series_link_outer_trace_width: float = series_link_outer_trace_width
        self.series_link_inner_offset: float = series_link_inner_offset
        self.series_link_outer_offset: float = series_link_outer_offset
        self.link_com: bool = link_com
        self.com_link_trace_width: float = com_link_trace_width
        self.com_link_offset: float = com_link_offset
        self.coil_names_font_size: float = coil_names_font_size
        self.coil_names_offset: float = coil_names_offset
        self.magnets_diameter: float = magnets_diameter
        self.magnets_position_radius: float = magnets_position_radius
        self.draw_vias: bool = draw_vias
        self.draw_terminals: bool = draw_terminals
        self.draw_outline: bool = draw_outline
        self.draw_construction_geometry: bool = draw_construction_geometry
        self.draw_only_layers: list[str] = draw_only_layers
        self.draw_magnets: bool = draw_magnets
        self.draw_coil_names: bool = draw_coil_names
        self.svg_profile: str = svg_profile
        self.svg_scale: int = svg_scale

        # Computed parameters
        self.viewport_width: float = self.board_diameter * 1.1
        self.viewport_height: float = self.board_diameter * 1.1
        self.board_radius: float = self.board_diameter/2
        self.hole_radius: float = self.hole_diameter/2
        self.via_diameter_w_spacing: float = self.via_diameter + self.trace_spacing
        self.n_coils: int = self.n_phases * self.n_slots_per_phase
        if self.max_turns_per_layer is None:
            self.max_turns_per_layer: float = 1000
        max_coil_angle = 360.0 / self.n_coils
        if self.coil_angle is None:
            self.coil_angle: float = max_coil_angle
        else:
            self.coil_angle: float = min(max_coil_angle, coil_angle)
        match n_layers:
            case 2:
                self.copper_layers = ['top', 'bottom']
            case 4:
                self.copper_layers = ['top', 'in1', 'in2', 'bottom']
            case 6:
                self.copper_layers = ['top', 'in1', 'in2', 'in3', 'in4', 'bottom']
            case 8:
                self.copper_layers = ['top', 'in1', 'in2', 'in3', 'in4', 'in5', 'in6', 'bottom']
            case _:
                raise ValueError("The number of layers must be 2, 4, 6 or 8")
        if self.series_link_inner_trace_width is None:
            self.series_link_inner_trace_width = self.trace_width * 4
        if self.series_link_outer_trace_width is None:
            self.series_link_outer_trace_width = self.series_link_inner_trace_width
        if self.com_link_trace_width is None:
            self.com_link_trace_width = self.series_link_outer_trace_width
        if self.board_outer_margin is None:
            self.board_outer_margin = self.terminal_diameter + self.terminal_offset + self.trace_spacing * 2
        if self.board_inner_margin is None:
            self.board_inner_margin = self.inner_vias_offset + self.series_link_inner_trace_width + self.trace_spacing * 2
        coil_names_position_radius = ((self.board_radius - self.board_outer_margin) + (self.hole_radius + self.board_inner_margin)) / 2.0
        if coil_names_font_size is None:
            suggested_size_by_width = (2 * math.pi * coil_names_position_radius / self.n_coils) / 6.0
            suggested_size_by_height = ((self.board_radius - self.board_outer_margin) - (self.hole_radius + self.board_inner_margin)) / 10.0
            self.coil_names_font_size = round(min(suggested_size_by_width, suggested_size_by_height), 1)
        if coil_names_offset is None:
            self.coil_names_offset = self.via_diameter * 4
        self.coil_names_position_radius = coil_names_position_radius + self.coil_names_offset
        if self.magnets_position_radius is None:
            self.magnets_position_radius = ((self.board_radius - self.board_outer_margin) + (self.hole_radius + self.board_inner_margin)) / 2.0
        self.n_magnets = 2 * self.n_slots_per_phase
        if magnets_diameter is None:
            self.magnets_diameter = round((self.board_radius - self.board_outer_margin) - (self.hole_radius + self.board_inner_margin), 1)
        max_magnets_diameter = round((2 * math.pi * self.magnets_position_radius / self.n_magnets) * 0.9, 1)
        if self.magnets_diameter >= max_magnets_diameter:
            self.magnets_diameter = max_magnets_diameter
            if magnets_diameter is not None:
                print(f"Warning : the specified magnet diameter of {magnets_diameter}mm is too large, reducing to {self.magnets_diameter}mm")
        if self.rotation is None:
            self.rotation = 0.0

        # SVG style
        self.background_color = "#001023"
        self.construction_geometry_color: str = "#848484"
        self.construction_geometry_line_width: float = board_diameter / 500.
        t = self.construction_geometry_line_width
        self.construction_geometry_dashes: str = f"{t * 4} {t * 2} {t} {t * 2}"
        self.outline_color = "#D0D2CD"
        self.outline_line_width = board_diameter / 300
        self.copper_layers_color = {
            'top': "#C83434",
            'bottom': "#4D7FC4",
            'in1': "#7FC87F",
            'in2': "#CE7D2C",
            'in3': "#4FCBCB",
            'in4': "#DB628B",
            'in5': "#A7A5C6",
            'in6': "#28CCD9",
        }
        self.top_silk_color = "#EEEFA4"
        self.bottom_silk_color = "#E6B2A4"
        self.silk_line_width = board_diameter / 300
        self.silk_font_family = "Arial, sans-serif"
        self.via_color = "#ECECEC"
        self.via_hole_color = "#E3B72E"
        self.via_opacity = 1.0
        self.terminal_color = "#E6B631"
        self.terminal_hole_color = self.background_color
        self.terminal_opacity = 1.0
        self.magnets_color = "#F0F0F0"
        self.magnets_line_width = board_diameter / 200
        self.magnets_opacity = 1.0
        self.magnets_dashes = "none"
        self.svg_font_size_factor = 1.5
    
    def read(filename: str) -> Self:
        """Read and parse the given file as JSON, and return a Config"""
        with open(filename, 'r') as f:
            cfg = json.load(f)
        return Config.from_dict(cfg)
    
    def from_dict(cfg: dict) -> Self:
        """Parse the given dict and return it as a Config"""
        args = {}
        for option in Config.options:
            value = cfg
            for item in option['json'].split('.'):
                try:
                    value = value[item]
                except KeyError:
                    if option.get('required', False):
                        raise ValueError(f"Config option '{option['json']}' is required (type : {option['type'].__name__})")
                    else:
                        try:
                            value = option['default']
                            break
                        except KeyError:
                            raise ValueError(f"Config option '{option['json']}' is not required but doesn't have a default value")
            if not option.get('nullable', False) and value is None:
                raise ValueError(f"Invalid value 'null' for config option '{option['json']}' (type : {option['type'].__name__})")
            if 'enum' in option and value not in option['enum']:
                available_values_str = f", valid values are : {', '.join(map(lambda x: f"'{x}'", option['enum']))}"
                raise ValueError(f"Invalid value '{value}' for config option '{option['json']}' (type : {option['type'].__name__}){available_values_str}")
            if 'decoder' in option:
                try:
                    value = option['decoder'](value)
                except KeyError:
                    available_values_str = ''
                    if 'enum' in option:
                        available_values_str = f", valid values are : {', '.join(map(lambda x: f"'{x}'", option['enum']))}"
                    raise ValueError(f"Invalid value '{value}' for config option '{option['json']}' (type : {option['type'].__name__}){available_values_str}")
            if value is not None:
                try:
                    if option['type'] is float:
                        value = float(value)
                    elif option['type'] is int:
                        value = int(value)
                    elif option['type'] is bool:
                        value = bool(value)
                except (ValueError, TypeError):
                    raise ValueError(f"Invalid value '{value}' for config option '{option['json']}' (type : {option['type'].__name__})")
            args[option['name']] = value
        return Config(**args)
    
    def as_dict(self) -> dict:
        cfg = {}
        for option in Config.options:
            cfg_parent = cfg
            path = option['json'].split('.')
            for i in range(len(path) - 1):
                if not path[i] in cfg_parent:
                    cfg_parent[path[i]] = {}
                cfg_parent = cfg_parent[path[i]]
            value = getattr(self, option['name'])
            if 'encoder' in option:
                value = option['encoder'](value)
            cfg_parent[option['name']] = value
        return cfg
    
    def to_json(self) -> str:
        return json.dumps(self.as_dict(), indent=4)