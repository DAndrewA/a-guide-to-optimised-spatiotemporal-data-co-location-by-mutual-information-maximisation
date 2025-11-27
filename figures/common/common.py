"""Author: Andrew Martin
Creation date: 27/11/25

Script holding common variables for all figure creation, analysis, etc.
"""
import os

from matplotlib.pyplot import subplots

DIR_ROOT = os.path.join(
    os.environ["SCRATCH"]
)

# the size of inidividual panels in the results sections plots
PANEL_SIZE = (4,3)
get_figure_panel = lambda : subplots(1,1, figsize=PANEL_SIZE, layout="constrained")