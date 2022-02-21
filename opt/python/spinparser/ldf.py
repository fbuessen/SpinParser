import re
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy as np

class Arrow3D(FancyArrowPatch):
    """
    Draw 3d arrows. 
    Based on https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
    with fix for Matplotlib v3.5.0 from https://github.com/matplotlib/matplotlib/issues/21688
    """
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        return np.min(zs)

def plot(filename, interactions="all"):
    """
    Plot the contents of an ldf file. 

    Parameters
    ----------
    filename : string
        Input file path. 
    
    Attributes
    ----------
    interactions : {"all" | "none" | int}
        Specify the interactions to include in the plot. 
        If an integer is provided, plot all interactions that involve the lattice site with that identifier. 
    """
    #read ldf file
    sites = []
    bonds = []
    interacs = []
    with open(filename, "r") as file:
        for line in file:
            if re.search("(?<=<)site", line):
                id = int(re.search("(?<=\W)id=\"(\d+)\"", line)[1])
                x = float(re.search("(?<=\W)x=\"(-?\d+\.\d*)\"", line)[1])
                y = float(re.search("(?<=\W)y=\"(-?\d+\.\d*)\"", line)[1])
                z = float(re.search("(?<=\W)z=\"(-?\d+\.\d*)\"", line)[1])
                parametrized = True if re.search("(?<=\W)parametrized=\"(true|false)\"", line)[1] == "true" else False
                sites.append({"id":id, "x":x, "y":y, "z":z, "parametrized":parametrized})
            elif re.search("(?<=<)bond", line):
                fromId = int(re.search("(?<=\W)from=\"(\d+)\"", line)[1])
                toId = int(re.search("(?<=\W)to=\"(\d+)\"", line)[1])
                bonds.append({"from":fromId, "to":toId})
            elif re.search("(?<=<)interaction", line):
                fromId = int(re.search("(?<=\W)from=\"(\d+)\"", line)[1])
                toId = int(re.search("(?<=\W)to=\"(\d+)\"", line)[1])
                value = eval(re.search("(?<=\W)value=\"([\[\d\.\]\,\-]+)\"", line)[1])
                interacs.append({"from":fromId, "to":toId, "value":value})

    #prepare plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.axis('off')

    #plot lattice bonds
    for bond in bonds:
        s1 = next(i for i in sites if i["id"]==bond["from"])
        s2 = next(i for i in sites if i["id"]==bond["to"])
        ax.plot3D([s1["x"],s2["x"]], [s1["y"],s2["y"]], [s1["z"],s2["z"]], 'gray')

    #plot lattice sites
    ax.scatter([i["x"] for i in sites], [i["y"] for i in sites], [i["z"] for i in sites], c=[("#ff7f0e" if i["parametrized"]==True else "#1f77b4") for i in sites], alpha=1.0)

    #plot interactions
    for interaction in interacs:
        s1 = next(i for i in sites if i["id"]==interaction["from"])
        s2 = next(i for i in sites if i["id"]==interaction["to"])

        if interactions == "all" or interaction["from"] == interactions or interaction["to"] == interactions:
            a = Arrow3D([s1["x"], s2["x"]], [s1["y"], s2["y"]], [s1["z"], s2["z"]], mutation_scale=10, lw=3, arrowstyle="-|>", color="r")
            ax.add_artist(a)

            label = "%6.2f %6.2f %6.2f\n%6.2f %6.2f %6.2f\n%6.2f %6.2f %6.2f" % (interaction["value"][0][0], interaction["value"][0][1], interaction["value"][0][2], interaction["value"][1][0], interaction["value"][1][1], interaction["value"][1][2], interaction["value"][2][0], interaction["value"][2][1], interaction["value"][2][2])
            ax.text(0.5*(s1["x"]+s2["x"]), 0.5*(s1["y"]+s2["y"]), 0.5*(s1["z"]+s2["z"]), label, fontsize="xx-small", horizontalalignment='center', verticalalignment='center',zorder=float('inf'))

    return plt