import cmcrameri.cm as cm
import numpy as np
from matplotlib.colors import LogNorm, Normalize, ListedColormap, CenteredNorm
from matplotlib.cm import ScalarMappable


CMAP_probability = ListedColormap(
    cm.acton(np.linspace(0,0.8,128)),
    name="reduced_acton"
)

CMAP_N_events = cm.navia
CMAP_N_profiles = cm.imola

CMAP_MI = ListedColormap(
    cm.lipari(np.linspace(0, 0.85, 128)),
    name="reduced_lipari"
)

CMAP_copula = cm.vik

NORM_probability = Normalize(vmin=0, vmax=1)
NORM_copula = CenteredNorm(vcenter=1, halfrange=0.75)