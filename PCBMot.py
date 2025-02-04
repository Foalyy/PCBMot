from typing import Self
from inkscape_drawing import InkscapeDrawing
from config import Config, TerminalType
from pcb import PCB


# Config
SVG_FILENAME = 'example.svg'
config = Config(
    board_diameter = 90, # mm
    hole_diameter = 55, # mm
    board_outer_margin = 1.8, # mm
    board_inner_margin = 1.0, # mm
    n_phases = 3,
    n_slots_per_phase = 10,
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
    magnets_diameter = None,
    magnets_position_radius = None,
    draw_vias = True,
    draw_terminals = True,
    draw_outline = True,
    draw_construction_geometry = True,
    draw_only_layers = None,
    draw_magnets = True,
    draw_coil_names = True,
    svg_profile = 'tiny'
)

# Create the PCB
pcb = PCB.generate(config)

# Create the SVG drawing
dwg = InkscapeDrawing(
    SVG_FILENAME,
    size = ("1000px", "1000px"),
    coords = (-config.viewport_width/2.0, -config.viewport_height/2.0, config.viewport_width, config.viewport_height),
    profile = config.svg_profile,
    background_color = config.background_color,
)

# Draw the PCB
pcb.draw_svg(dwg)

# Save the final SVG file
dwg.save()