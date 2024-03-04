"""This is the colormap module.

The colormap module contains lists of perceptually uniform, scientific 
colormaps and a method for returning a specified number of colors
from a specified colormap.
The maximum number of returnet colors is 255.

Methods:
------
returnColors (int, str):
    Return evenly distributed number of colors from a given colormap.
"""

__version__ = '1.0'
__author__ = 'Christiane Rohse'

import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm     

color_maps = ['viridis', 'viridis_r', 'plasma', 'plasma_r', 'batlow', 'batlow_r', 'hawaii', 'hawaii_r', 'cividis', 'cividis_r']
color_map_grey = ['grayC', 'grayC_r']

def returnColors(n, cmap):
    """Return a specified number of colors selected with even distances
    from a specified color map.

    Parameter:
    ------
    n (int): 
        Number of colors to return.
    cmap (str):
        Name of colormap to select from.
    Return:
    ------
    List of list containing four float values: for r, g, b and a
    """
    try:
        c = 'cm.'+cmap
        cmap = eval(c)
    except:
        cmap = plt.get_cmap(cmap)
    if n > 256:
        pass
    else:
        return (cmap(np.linspace(0, 1, n)))