from typing import Self
from inkscape_drawing import InkscapeDrawing
from config import Config, BoardShape, TerminalType
from pcb import PCB
from kicad import KicadPCB


# Config
KICADPCB_FILENAME = 'example.kicad_pcb'
SVG_FILENAME = 'example.svg'
config = Config(
    board_shape = BoardShape.SQUARE,
    board_diameter = 42, # mm
    hole_diameter = 18, # mm
    board_chamfer = 4, # mm
    board_fillet = 3, # mm
    n_mountpoints = 4,
    mountpoints_position_radius = 23, # mm
    mountpoints_diameter = 3, # mm
    mountpoints_marking_diameter = 5.5, # mm
    board_thickness = 1.0, #mm
    board_outer_margin = 3, # mm
    board_inner_margin = 1.0, # mm
    n_phases = 3,
    n_slots_per_phase = 4,
    coil_angle = None,
    n_layers = 8,
    four_layers_inside_vias = False,
    max_turns_per_layer = None,
    trace_width = 0.1, # mm
    trace_spacing = 0.1, # mm
    via_diameter = 0.45, # mm
    via_hole_diameter = 0.3, # mm
    outer_vias_offset = 0.0, # mm
    inner_vias_offset = 0.4, # mm
    terminal_type = TerminalType.THROUGH_HOLE,
    terminal_diameter = 1.5, # mm
    terminal_hole_diameter = 0.8, # mm
    terminal_offset = 0.0, # mm
    link_series_coils = True,
    series_link_inner_trace_width = 0.3, # mm
    series_link_outer_trace_width = 0.6, # mm
    series_link_inner_offset = 0.0, # mm
    series_link_outer_offset = 0.0, # mm
    link_com = True,
    com_link_trace_width = 0.6, # mm
    com_link_offset = 0.0, # mm
    coil_names_font_size = None,
    coil_names_offset = None,
    rotation = -15,
    magnets_diameter = None,
    magnets_position_radius = None,
    draw_vias = True,
    draw_terminals = True,
    draw_outline = True,
    draw_construction_geometry = True,
    draw_only_layers = None,
    draw_magnets = True,
    draw_coil_names = True,
    svg_profile = 'tiny',
    svg_scale = 100, # px/mm
)

# Generate the PCB based on the given config
pcb = PCB.generate(config)

# Create the SVG drawing
dwg = InkscapeDrawing(
    SVG_FILENAME,
    size = (f"{int(config.viewport_width * config.svg_scale)}px", f"{int(config.viewport_height * config.svg_scale)}px"),
    coords = (-config.viewport_width/2.0, -config.viewport_height/2.0, config.viewport_width, config.viewport_height),
    profile = config.svg_profile,
    background_color = config.background_color,
)
pcb.draw_svg(dwg)
dwg.save()

# Create the Kicad board
kicadpcb = KicadPCB(KICADPCB_FILENAME, config)
pcb.draw_kicad(kicadpcb)
kicadpcb.save()