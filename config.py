from enum import Enum
from typing import Self
import math, json
from .geometry import Point, Circle

class BoardShape(Enum):
    CIRCLE = 1
    SQUARE = 2

class RadialTraces(Enum):
    PARALLEL = 1
    RADIAL = 2

class TracesGeometry(Enum):
    LINES = 1
    POLYGONS = 2

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
            'signed': True,
            'zero': True,
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
            'zero': True,
            'default': 0,
        },
        {
            'name': 'board_fillet',
            'json': 'board.board_fillet',
            'type': float,
            'zero': True,
            'default': 0,
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
            'signed': True,
            'zero': True,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'board_inner_margin',
            'json': 'board.board_inner_margin',
            'type': float,
            'signed': True,
            'zero': True,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'rotation',
            'json': 'board.rotation',
            'type': float,
            'signed': True,
            'zero': True,
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
            'name': 'generate_mountpoints',
            'json': 'board.mountpoints.generate_mountpoints',
            'type': bool,
            'default': False,
        },
        {
            'name': 'n_mountpoints',
            'json': 'board.mountpoints.n_mountpoints',
            'type': int,
            'default': 4,
        },
        {
            'name': 'mountpoints_position_radius',
            'json': 'board.mountpoints.mountpoints_position_radius',
            'type': float,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'mountpoints_diameter',
            'json': 'board.mountpoints.mountpoints_diameter',
            'type': float,
            'default': 3.0,
        },
        {
            'name': 'mountpoints_marking_diameter',
            'json': 'board.mountpoints.mountpoints_marking_diameter',
            'type': float,
            'default': 5.0,
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
            'name': 'n_coils_per_phase',
            'json': 'coils.n_coils_per_phase',
            'type': int,
            'required': True,
        },
        {
            'name': 'coil_angle',
            'json': 'coils.coil_angle',
            'type': float,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'turns_per_layer',
            'json': 'coils.turns_per_layer',
            'type': int,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'radial_traces',
            'json': 'coils.radial_traces',
            'type': RadialTraces,
            'decoder': lambda s : {
                'parallel': RadialTraces.PARALLEL,
                'radial': RadialTraces.RADIAL,
            }[s],
            'encoder': lambda s : {
                RadialTraces.PARALLEL: 'parallel',
                RadialTraces.RADIAL: 'radial',
            }[s],
            'enum': ['parallel', 'radial'],
            'default': 'parallel',
        },
        {
            'name': 'traces_geometry',
            'json': 'coils.traces_geometry',
            'type': TracesGeometry,
            'decoder': lambda s : {
                'lines': TracesGeometry.LINES,
                'polygons': TracesGeometry.POLYGONS,
            }[s],
            'encoder': lambda s : {
                TracesGeometry.LINES: 'lines',
                TracesGeometry.POLYGONS: 'polygons',
            }[s],
            'enum': ['lines', 'polygons'],
            'default': 'lines',
        },
        {
            'name': 'coil_names_font_size',
            'json': 'coils.coil_names_font_size',
            'type': float,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'coil_names_offset',
            'json': 'coils.coil_names_offset',
            'type': float,
            'zero': True,
            'auto': True,
            'default': 'auto',
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
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'series_link_outer_trace_width',
            'json': 'coils.series_link.series_link_outer_trace_width',
            'type': float,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'series_link_inner_offset',
            'json': 'coils.series_link.series_link_inner_offset',
            'type': float,
            'signed': True,
            'zero': True,
            'default': 0.0,
        },
        {
            'name': 'series_link_outer_offset',
            'json': 'coils.series_link.series_link_outer_offset',
            'type': float,
            'signed': True,
            'zero': True,
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
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'com_link_offset',
            'json': 'coils.com_link.com_link_offset',
            'type': float,
            'signed': True,
            'zero': True,
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
            'signed': True,
            'zero': True,
            'default': 0.0,
        },
        {
            'name': 'inner_vias_offset',
            'json': 'vias.inner_vias_offset',
            'type': float,
            'signed': True,
            'zero': True,
            'default': 0.0,
        },
        {
            'name': 'inside_vias_offset',
            'json': 'vias.inside_vias_offset',
            'type': float,
            'signed': True,
            'zero': True,
            'default': 0.0,
        },
        {
            'name': 'four_layers_inner_vias',
            'json': 'vias.four_layers_inner_vias',
            'type': bool,
            'default': False,
        },
        {
            'name': 'via_resistance',
            'json': 'vias.via_resistance',
            'type': float,
            'zero': True,
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
            'zero': True,
            'default': 0.0,
        },
        {
            'name': 'generate_com_terminal',
            'json': 'terminals.generate_com_terminal',
            'type': bool,
            'default': True,
        },

        # Magnets
        {
            'name': 'magnets_diameter',
            'json': 'magnets.magnets_diameter',
            'type': float,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'magnets_position_radius',
            'json': 'magnets.magnets_position_radius',
            'type': float,
            'auto': True,
            'default': 'auto',
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
            'null': True,
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

        # Kicad
        {
            'name': 'arcs_discretization_length',
            'json': 'kicad.arcs_discretization_length',
            'type': float,
            'auto': True,
            'default': 'auto',
        },
        {
            'name': 'arcs_discretization_angle',
            'json': 'kicad.arcs_discretization_angle',
            'type': float,
            'default': 2, # degrees
        },
    ]

    def __init__(self, cfg: dict):
        """Parse the given dict containing settings and return it as a Config"""

        self.auto_parameters = []
        for option in Config.options:
            # Get the value for this option in the dict structure
            value = cfg
            path = []
            for item in option['json'].split('.'):
                path.append(item)
                try:
                    value = value[item]
                except TypeError:
                    raise ValueError(f"Invalid value '{value}' for config option {'.'.join(path[:-1])}")
                except KeyError:
                    if option.get('required', False):
                        raise ValueError(f"Config option '{option['json']}' is required (type : {option['type'].__name__})")
                    else:
                        try:
                            value = option['default']
                            break
                        except KeyError:
                            raise ValueError(f"Config option '{option['json']}' is not required but doesn't have a default value")
            
            # Check if the value is allowed to be 'null'
            if not option.get('null', False) and value is None:
                raise ValueError(f"Invalid value 'null' for config option '{option['json']}' (type : {option['type'].__name__})")
            
            # Check if the value is allowed to be 'auto'
            if value == 'auto':
                if option.get('auto', False):
                    # Keep track of parameters set to 'auto'
                    self.auto_parameters.append(option['name'])
                else:
                    raise ValueError(f"Invalid value 'auto' for config option '{option['json']}' (type : {option['type'].__name__})")

            # Check value for enum types
            if 'enum' in option and value not in option['enum']:
                raise ValueError(f"Invalid value '{value}' for config option '{option['json']}' (type : {option['type'].__name__}), valid values are : {', '.join(map(lambda x: f"'{x}'", option['enum']))}")

            # If a decoder handler is specified for this option, execute it
            if 'decoder' in option:
                try:
                    value = option['decoder'](value)
                except KeyError:
                    available_values_str = ''
                    if 'enum' in option:
                        available_values_str = f", valid values are : {', '.join(map(lambda x: f"'{x}'", option['enum']))}"
                    raise ValueError(f"Invalid value '{value}' for config option '{option['json']}' (type : {option['type'].__name__}){available_values_str}")
            
            # If a real value is given (not 'null' or 'auto'), check it against its type
            if value is not None and value != 'auto':
                try:
                    # Float
                    if option['type'] is float:
                        value = float(value)

                        # Negative and zero values are not allowed by default
                        if not option.get('signed', False) and value < 0:
                            raise ValueError(f"Invalid value negative '{value}' for config option '{option['json']}' (type : {option['type'].__name__})")
                        if not option.get('zero', False) and math.isclose(value, 0, abs_tol=1e-9):
                            raise ValueError(f"Invalid null value '{value}' for config option '{option['json']}' (type : {option['type'].__name__})")
                    
                    # Int
                    elif option['type'] is int:
                        value = int(value)

                        # Negative and zero values are not allowed by default
                        if not option.get('signed', False) and value < 0:
                            raise ValueError(f"Invalid value negative '{value}' for config option '{option['json']}' (type : {option['type'].__name__})")
                        if not option.get('zero', False) and value == 0:
                            raise ValueError(f"Invalid null value '{value}' for config option '{option['json']}' (type : {option['type'].__name__})")
                    
                    # Bool
                    elif option['type'] is bool:
                        value = bool(value)

                except (ValueError, TypeError):
                    raise ValueError(f"Invalid value '{value}' for config option '{option['json']}' (type : {option['type'].__name__})")
            
            # Save this option directly in the Config object
            setattr(self, option['name'], value)


        ## Other specific checks on some options
        
        if self.hole_diameter >= self.board_diameter:
            raise ValueError("The board diameter must be larger than the center hole diameter.")
        
        if self.generate_com_terminal and not (self.n_coils_per_phase >= 1 and (self.n_coils_per_phase == 1 or self.n_coils_per_phase % 2 == 0)):
            raise ValueError("In order to generate a COM terminal, the number of coils per phase must be 1 or an even number, either change the number of coils per phase or disable the COM terminal.")
        
        if self.link_series_coils and self.n_phases > self.n_layers:
            raise ValueError("Unable to generate the connections between the coils because the number of layers should be greater than or equal to the number of phases, either increase the number of layers or disable link_series_coils.")
        
        if self.terminal_type == TerminalType.CASTELLATED and self.board_shape != BoardShape.CIRCLE:
            raise ValueError("Castellated terminals are only possible for circle boards.")

        # Names of the copper layers based on the number of layers
        match self.n_layers:
            case 2:
                self.copper_layers = ['top', 'bottom']
            case 4:
                self.copper_layers = ['top', 'in1', 'in2', 'bottom']
            case 6:
                self.copper_layers = ['top', 'in1', 'in2', 'in3', 'in4', 'bottom']
            case 8:
                self.copper_layers = ['top', 'in1', 'in2', 'in3', 'in4', 'in5', 'in6', 'bottom']
            case _:
                raise ValueError("The number of layers must be 2, 4, 6 or 8.")
        

        ## Computed parameters, mainly parameters set to 'auto'

        # Set the viewport size slightly larger than the board
        self.viewport_width: float = self.board_diameter * 1.1
        self.viewport_height: float = self.board_diameter * 1.1

        # Compute radiuses based on diameters to make computations easier later
        self.board_radius: float = self.board_diameter/2
        self.hole_radius: float = self.hole_diameter/2
        self.via_diameter_w_spacing: float = self.via_diameter + self.trace_spacing

        # Mountpoints in the corners by default, half distance between the chamfer (or the corner) and the board outer radius
        if self.generate_mountpoints and self.mountpoints_position_radius == 'auto':
            if self.board_shape == BoardShape.SQUARE:
                chamfer_radius = self.board_radius * math.sqrt(2) - self.board_chamfer * math.sqrt(2) / 2.0
                self.mountpoints_position_radius = round((chamfer_radius + self.board_radius) / 2.0, 1)
            else:
                self.mountpoints_position_radius = self.board_radius

        # Total number of coils
        self.n_coils: int = self.n_phases * self.n_coils_per_phase
        
        # Coil angle limited by the the number of coils, and set to the max by default
        max_coil_angle = 360.0 / self.n_coils
        if self.coil_angle == 'auto':
            self.coil_angle: float = max_coil_angle
        else:
            self.coil_angle: float = min(max_coil_angle, self.coil_angle)
        
        # Link and COM trace width set to 4 times the coils trace width by default, but capped by via diameter
        # to prevent interference with the coils around the corner vias
        if self.series_link_inner_trace_width == 'auto':
            self.series_link_inner_trace_width = min(self.trace_width * 4, self.via_diameter)
        if self.series_link_outer_trace_width == 'auto':
            self.series_link_outer_trace_width = self.series_link_inner_trace_width
        if self.com_link_trace_width == 'auto':
            if self.n_coils_per_phase % 2 == 0:
                self.com_link_trace_width = self.series_link_outer_trace_width
            else:
                self.com_link_trace_width = self.series_link_inner_trace_width
        
        # Default board margins are set to take into account the terminals and link traces
        if self.board_outer_margin == 'auto':
            if self.terminal_type != TerminalType.NONE:
                margin_terminals = self.trace_spacing + self.terminal_offset + self.terminal_diameter + self.trace_spacing
            else:
                margin_terminals = 0
            if self.link_series_coils:
                margin_series_link = self.trace_spacing + self.via_diameter + self.trace_spacing + self.series_link_outer_offset + self.series_link_outer_trace_width + self.trace_spacing
            else:
                margin_series_link = 0
            if self.link_com:
                margin_com_link = self.trace_spacing + self.com_link_offset + self.com_link_trace_width + self.trace_spacing
            else:
                margin_com_link = 0
            self.board_outer_margin = round(max(margin_terminals, margin_series_link, margin_com_link), 2)
        if self.board_inner_margin == 'auto':
            self.board_inner_margin = round(self.trace_spacing + self.inner_vias_offset + self.trace_spacing + self.series_link_inner_trace_width + self.trace_spacing, 2)
        self.coils_outer_radius = self.board_radius - self.board_outer_margin
        self.coils_inner_radius = self.hole_radius + self.board_inner_margin
        self.coils_middle_radius = (self.coils_outer_radius + self.coils_inner_radius) / 2.0

        # Find a font size that fits both in height and width for the the coils names by default
        if self.coil_names_font_size == 'auto':
            suggested_size_by_width = (Circle(Point.origin(), self.coils_middle_radius).perimeter() / self.n_coils) / 6.0
            suggested_size_by_height = (self.coils_outer_radius - self.coils_inner_radius) / 10.0
            self.coil_names_font_size = round(min(suggested_size_by_width, suggested_size_by_height), 1)
        
        # Place the coils names between the inner vias and the outside edge of the coils by default
        if self.coil_names_offset == 'auto':
            self.coil_names_offset = (self.coils_outer_radius - self.coils_middle_radius) / 2.0
        self.coil_names_position_radius = self.coils_middle_radius + self.coil_names_offset
        
        # Center the magnets with the coils by default
        if self.magnets_position_radius == 'auto':
            self.magnets_position_radius = self.coils_middle_radius
        
        # Put 2 magnets for each coil in series by default, as required by the most common motor configuration
        self.n_magnets = 2 * self.n_coils_per_phase

        # Set the magnets diameter to the height of the coil by default, while making sure that they are not overlapping, and not larger than the coils
        max_magnets_diameter = round((Circle(Point.origin(), self.magnets_position_radius).perimeter() / self.n_magnets) * 0.9, 1)
        coil_width = Circle(Point.origin(), self.coils_middle_radius).perimeter() / self.n_coils
        max_magnets_diameter = min(coil_width, max_magnets_diameter)
        if self.magnets_diameter == 'auto':
            self.magnets_diameter = min(round(self.coils_outer_radius - self.coils_inner_radius, 1), max_magnets_diameter)
        elif self.magnets_diameter > max_magnets_diameter:
            print(f"Warning : the specified magnet diameter of {magnets_diameter}mm is too large")
        
        # No board rotation by default
        if self.rotation == 'auto':
            self.rotation = 0.0
        
        # Set arcs discretization to the trace width as a reasonable default value
        if self.arcs_discretization_length == 'auto':
            self.arcs_discretization_length = self.trace_width


        ## SVG style
        self.background_color = "#001023"
        self.construction_geometry_color: str = "#848484"
        self.construction_geometry_line_width: float = self.board_diameter / 500.
        t = self.construction_geometry_line_width
        self.construction_geometry_dashes: str = f"{t * 4} {t * 2} {t} {t * 2}"
        self.outline_color = "#D0D2CD"
        self.outline_line_width = self.board_diameter / 300
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
        self.silk_line_width = self.board_diameter / 300
        self.silk_font_family = "Arial, sans-serif"
        self.via_color = "#ECECEC"
        self.via_hole_color = "#E3B72E"
        self.via_opacity = 1.0
        self.terminal_color = "#E6B631"
        self.terminal_hole_color = self.background_color
        self.terminal_opacity = 1.0
        self.magnets_color = "#F0F0F0"
        self.magnets_line_width = self.board_diameter / 200
        self.magnets_opacity = 1.0
        self.magnets_dashes = "none"
        self.svg_font_size_factor = 1.5
    
    def read(filename: str) -> Self:
        """Read and parse the given file as JSON, and return a Config"""
        with open(filename, 'r') as f:
            cfg = json.load(f)
        return Config(cfg)
    
    def as_dict(self) -> dict:
        cfg = {}
        for option in Config.options:
            cfg_parent = cfg
            path = option['json'].split('.')
            for i in range(len(path) - 1):
                if not path[i] in cfg_parent:
                    cfg_parent[path[i]] = {}
                cfg_parent = cfg_parent[path[i]]
            if option.get('auto', False) and option['name'] in self.auto_parameters:
                cfg_parent[option['name'] + '_auto'] = True
            value = getattr(self, option['name'])
            if 'encoder' in option:
                value = option['encoder'](value)
            cfg_parent[option['name']] = value
        return cfg
    
    def to_json(self) -> str:
        return json.dumps(self.as_dict(), indent=4)