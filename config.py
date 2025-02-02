from enum import Enum

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
        n_phases: int,
        n_slots_per_phase: int,
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
        draw_vias: bool,
        draw_terminals: bool,
        draw_construction_geometry: bool,
        draw_only_layers: list[str],
        svg_profile: str,
    ):
        # Check parameters
        if n_layers not in [2, 4, 6, 8]:
            raise ValueError("The number of layers must be 2, 4, 6 or 8")
        
        # Save parameters
        self.board_diameter: float = board_diameter
        self.hole_diameter: float = hole_diameter
        self.n_phases: int = n_phases
        self.n_slots_per_phase: int = n_slots_per_phase
        self.n_layers: int = n_layers
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
        self.draw_vias: bool = draw_vias
        self.draw_terminals: bool = draw_terminals
        self.draw_construction_geometry: bool = draw_construction_geometry
        self.draw_only_layers: list[str] = draw_only_layers
        self.svg_profile = svg_profile

        # Computed parameters
        self.viewport_width: float = self.board_diameter * 1.1
        self.viewport_height: float = self.board_diameter * 1.1
        self.board_radius: float = self.board_diameter/2
        self.hole_radius: float = self.hole_diameter/2
        self.board_outer_margin: float = 1.8
        self.board_inner_margin: float = 1.0
        self.via_diameter_w_spacing: float = self.via_diameter + self.trace_spacing
        self.n_coils: int = self.n_phases * self.n_slots_per_phase
        self.coil_angle: float = 360.0 / self.n_coils

        # SVG style
        self.background_color = "#001023"
        self.construction_geometry_color: str = "#848484"
        self.construction_geometry_thickness: float = board_diameter / 1000.
        self.construction_geometry_dashes: str = "0.2 0.12"
        self.outline_thickness = 0.1
        self.layers_color = {
            'top': "#C83434",
            'bottom': "#4D7FC4",
            'in1': "#7FC87F",
            'in2': "#CE7D2C",
            'in3': "#4FCBCB",
            'in4': "#DB628B",
            'in5': "#A7A5C6",
            'in6': "#28CCD9",
            'outline': "#D0D2CD",
        }
        self.via_color = "#ECECEC"
        self.via_hole_color = "#E3B72E"
        self.terminal_color = "#E6B631"
        self.terminal_hole_color = self.background_color