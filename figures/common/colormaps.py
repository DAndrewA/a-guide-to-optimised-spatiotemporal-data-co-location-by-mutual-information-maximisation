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
NORM_bias = LogNorm(vmax=np.power(10,1.5), vmin=np.power(10,-2.5), clip=True)

COLOR_ATL09 = cm.roma(0.65)
COLOR_Cloudnet = cm.roma(0.2)

MAPPABLE_copula = ScalarMappable(cmap=CMAP_copula, norm=NORM_copula)
MAPPABLE_probability = ScalarMappable(cmap=CMAP_probability, norm=NORM_probability)
MAPPABLE_bias = ScalarMappable(cmap=CMAP_probability, norm=NORM_bias)

#CMAP_plabels = cm.hawaii_r
CMAP_plabels = cm.batlowS
