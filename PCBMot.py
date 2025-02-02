from typing import Self
import svgwrite as svg
from config import Config, TerminalType
from pcb import PCB


# Config
SVG_FILENAME = 'example.svg'
config = Config(
    board_diameter = 23, # mm
    hole_diameter = 4, # mm
    board_outer_margin = 1.8, # mm
    board_inner_margin = 1.0, # mm
    n_phases = 3,
    n_slots_per_phase = 2,
    n_layers = 8,
    max_turns_per_layer = None,
    trace_width = 0.1, # mm
    trace_spacing = 0.1, # mm
    via_diameter = 0.45, # mm
    via_hole_diameter = 0.3, # mm
    outer_vias_offset = 0.0, # mm
    inner_vias_offset = 0.2, # mm
    terminal_type = TerminalType.THROUGH_HOLE,
    terminal_diameter = 1.5, # mm
    terminal_hole_diameter = 0.8, # mm
    terminal_offset = 0.0, # mm
    link_series_coils = True,
    series_link_trace_width = 0.3, # mm
    series_link_offset = 0.0, # mm
    link_com = True,
    com_link_trace_width = 0.8, # mm
    com_link_offset = 0.0, # mm
    draw_vias = True,
    draw_terminals = True,
    draw_construction_geometry = True,
    draw_only_layers = None,
    svg_profile = 'tiny'
)

# Create the PCB
pcb = PCB.generate(config)

# Create an SVG drawing
dwg = svg.Drawing(
    SVG_FILENAME,
    size=("1000px", "1000px"),
    viewBox=f"{-config.viewport_width/2.0} {-config.viewport_height/2.0} {config.viewport_width} {config.viewport_height}",
    profile=config.svg_profile,
    fill="none",
    stroke="#000000",
    stroke_width=0.5,
    stroke_linejoin="round",
    stroke_linecap="round",
)
dwg.add(dwg.rect( # Background
    (-config.viewport_width/2.0, -config.viewport_height/2.0),
    (config.viewport_width, config.viewport_height),
    stroke = "none",
    fill = config.background_color,
))

# Draw the PCB
pcb.draw_svg(dwg)

# Save the final SVG file
dwg.save()