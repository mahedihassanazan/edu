# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///

import os
import sys

# Set environment variable for headless rendering
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol # pytype: disable=import-error
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd # pytype: disable=import-error

cmd.load("/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/AF-O70370-F1-model_v6.cif", "protein")

if cmd.count_atoms("all") == 0:
    print("Error: No atoms loaded.")
    cmd.quit()
    sys.exit(1)

cmd.show("cartoon")

# Define AlphaFold's canonical pLDDT colors
cmd.set_color("af_very_low", [0xFF, 0x7D, 0x45])   # orange, pLDDT < 50
cmd.set_color("af_low",      [0xFF, 0xDB, 0x13])   # yellow, 50 <= pLDDT < 70
cmd.set_color("af_confident", [0x65, 0xCB, 0xF3])  # light blue, 70 <= pLDDT < 90
cmd.set_color("af_very_high", [0x00, 0x53, 0xD6])  # dark blue, pLDDT >= 90

# Apply from lowest to highest threshold
cmd.color("af_very_low", "polymer.protein")
cmd.color("af_low", "polymer.protein and b > 50")
cmd.color("af_confident", "polymer.protein and b > 70")
cmd.color("af_very_high", "polymer.protein and b > 90")

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.orient()

cmd.png("/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/CTSS_AlphaFold.png", width=1200, height=900, dpi=150)
cmd.save("/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/CTSS_session.pse")
cmd.quit()
