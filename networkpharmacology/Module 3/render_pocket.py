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

# Load full protein
cmd.load("/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/CTSS.pdb", "protein")

# Load the pocket alpha spheres
cmd.load("/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/CTSS_out/pockets/pocket8_vert.pqr", "pocket8")

if cmd.count_atoms("all") == 0:
    print("Error: No atoms loaded.")
    cmd.quit()
    sys.exit(1)

# Render protein surface
cmd.remove("solvent")
cmd.show("surface", "protein")
cmd.set("surface_color", "white", "protein")
cmd.set("transparency", 0.15, "protein")

# Show the pocket center (alpha spheres) as spheres
cmd.show("spheres", "pocket8")
cmd.color("red", "pocket8")
cmd.set("sphere_scale", 1.0, "pocket8")

# Orient and render to keep pocket in the center
cmd.orient("pocket8")
cmd.zoom("pocket8", buffer=20.0)
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)

output_png = "/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/CTSS_Pocket8.png"
output_pse = "/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/CTSS_Pocket8.pse"

cmd.png(output_png, width=1200, height=900, dpi=150)
cmd.save(output_pse)
cmd.quit()
