from typing import Self
from inkscape_drawing import InkscapeDrawing
from config import Config, BoardShape, TerminalType
from pcb import PCB
from kicad import KicadPCB


# Config
KICADPCB_FILENAME = 'examples/example.kicad_pcb'
SVG_FILENAME = 'examples/example.svg'
config = Config.read('examples/example_full.config')

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