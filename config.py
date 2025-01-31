## Board config
class Config:
    def __init__(
        self,
        board_diameter: float,
        hole_diameter: float,
        trace_width: float,
        trace_spacing: float,
        n_phases: int,
        n_slots_per_phase: int,
    ):
        self.board_diameter: float = board_diameter
        self.hole_diameter: float = hole_diameter
        self.trace_width: float = trace_width
        self.trace_spacing: float = trace_spacing
        self.n_phases: int = n_phases
        self.n_slots_per_phase: int = n_slots_per_phase

        # Computed parameters
        self.viewport_width: float = self.board_diameter * 1.1
        self.viewport_height: float = self.board_diameter * 1.1
        self.board_radius: float = self.board_diameter/2
        self.hole_radius: float = self.hole_diameter/2
        self.board_outer_margin: float = 1.8
        self.board_inner_margin: float = 1.0
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