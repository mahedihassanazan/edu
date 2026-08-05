# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///

import os
import sys

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

cmd.save("/home/hmalik342/Working/Python-BioPython/Demo Project/Module 3/CTSS.pdb", "protein")
print("Converted successfully to PDB.")
cmd.quit()
