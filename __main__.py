from typing import Self
from .inkscape_drawing import InkscapeDrawing
from .config import Config, BoardShape, TerminalType
from .pcb import PCB, PCBStats
from .kicad import KicadPCB
import argparse, sys

if __name__ == "__main__":
    # Read arguments from the command line
    parser = argparse.ArgumentParser(
        prog = 'PCBMot',
        description = "Generate PCB motors easily",
    )
    parser.add_argument('-c', '--config', required=True, help="Config file that describes the design to generate")
    parser.add_argument('-o', '--output', help="SVG output file")
    parser.add_argument('-k', '--kicad', help="Kicad PCB output file")
    parser.add_argument('-s', '--stats', help="Stats output file ('-' to print on the standard output)")
    args = parser.parse_args()
    if args.output is None and args.kicad is None:
        print("Please provide at least either the --output or the --kicad argument", file=sys.stderr)
        sys.exit(1)
    compute_stats = args.stats is not None

    # Read the config file
    config = Config.read(args.config)

    # Generate the PCB based on the given config
    pcb = PCB.generate(config, compute_stats)

    # Create the SVG drawing
    if args.output:
        dwg = InkscapeDrawing(
            args.output,
            size = (f"{int(config.viewport_width * config.svg_scale)}px", f"{int(config.viewport_height * config.svg_scale)}px"),
            coords = (-config.viewport_width/2.0, -config.viewport_height/2.0, config.viewport_width, config.viewport_height),
            profile = config.svg_profile,
            background_color = config.background_color,
        )
        pcb.draw_svg(dwg)
        dwg.save()

    # Create the Kicad board
    if args.kicad:
        kicadpcb = KicadPCB(args.kicad, config)
        pcb.draw_kicad(kicadpcb)
        kicadpcb.save()

    # Write the stats either in a file or on the standard output
    if compute_stats:
        pcb.stats.write(args.stats)