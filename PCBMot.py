from typing import Self
import svgwrite as svg
from config import Config
from pcb import PCB


# Config
SVG_FILENAME = 'example.svg'
config = Config(
    board_diameter = 23, # mm
    hole_diameter = 4, # mm
    trace_width = 0.1, # mm
    trace_spacing = 0.1, # mm
    n_phases = 3,
    n_slots_per_phase = 2,
)

# Create the PCB
pcb = PCB.generate(config)

# Create an SVG drawing
dwg = svg.Drawing(
    SVG_FILENAME,
    size=("800px", "800px"),
    viewBox=f"{-config.viewport_width/2.0} {-config.viewport_height/2.0} {config.viewport_width} {config.viewport_height}",
    profile='tiny',
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
pcb.draw_svg(dwg, construction_geometry=True)

# Save the final SVG file
dwg.save()