from enum import Enum
import math

class TerminalType(Enum):
    NONE = 1
    THROUGH_HOLE = 2
    SMD = 3
    CASTELLATED = 4

## Board config
class Config:
    def __init__(
        self,
        board_diameter: float,
        hole_diameter: float,
        board_thickness: float,
        board_outer_margin: float,
        board_inner_margin: float,
        n_phases: int,
        four_layers_inside_vias: bool,
        n_slots_per_phase: int,
        coil_angle: float,
        n_layers: int,
        max_turns_per_layer: int,
        trace_width: float,
        trace_spacing: float,
        via_diameter: float,
        via_hole_diameter: float,
        outer_vias_offset: float,
        inner_vias_offset: float,
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
        if hole_diameter + board_inner_margin >= board_diameter - board_outer_margin:
            raise ValueError("The board diameter is smaller than the hole diameter")
        if not (n_slots_per_phase >= 1 and (n_slots_per_phase == 1 or n_slots_per_phase % 2 == 0)):
            raise ValueError("The number of slots per phase must be 1 or an even number")
        if link_series_coils and n_phases > n_layers:
            raise ValueError("Unable to generate the connections between the coils because the number of layers should be greater than or equal to the number of phases, either increase the number of layers or disable link_series_coils")
        if terminal_type not in [TerminalType.NONE, TerminalType.THROUGH_HOLE, TerminalType.SMD, TerminalType.CASTELLATED]:
            raise ValueError("Invalid terminal type")
        
        # Save parameters
        self.board_diameter: float = board_diameter
        self.hole_diameter: float = hole_diameter
        self.board_thickness: float = board_thickness
        self.board_outer_margin: float = board_outer_margin
        self.board_inner_margin: float = board_inner_margin
        self.n_phases: int = n_phases
        self.n_slots_per_phase: int = n_slots_per_phase
        self.coil_angle: float = coil_angle
        self.n_layers: int = n_layers
        self.four_layers_inside_vias: bool = four_layers_inside_vias
        self.max_turns_per_layer: int = max_turns_per_layer
        self.trace_width: float = trace_width
        self.trace_spacing: float = trace_spacing
        self.via_diameter: float = via_diameter
        self.via_hole_diameter: float = via_hole_diameter
        self.outer_vias_offset: float = outer_vias_offset
        self.inner_vias_offset: float = inner_vias_offset
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
        if self.coil_angle is None:
            self.coil_angle: float = 360.0 / self.n_coils
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

        # SVG style
        self.background_color = "#001023"
        self.construction_geometry_color: str = "#848484"
        self.construction_geometry_thickness: float = board_diameter / 500.
        t = self.construction_geometry_thickness
        self.construction_geometry_dashes: str = f"{t * 4} {t * 2} {t} {t * 2}"
        self.outline_color = "#D0D2CD"
        self.outline_thickness = board_diameter / 300
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
        self.silk_font_family = "Arial, sans-serif"
        self.via_color = "#ECECEC"
        self.via_hole_color = "#E3B72E"
        self.via_opacity = 1.0
        self.terminal_color = "#E6B631"
        self.terminal_hole_color = self.background_color
        self.terminal_opacity = 1.0
        self.magnets_color = "#F0F0F0"
        self.magnets_thickness = board_diameter / 200
        self.magnets_opacity = 1.0
        self.magnets_dashes = "none"
        self.svg_font_size_factor = 1.5