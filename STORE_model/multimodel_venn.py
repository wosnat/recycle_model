#!/usr/bin/env python
# coding: utf-8


"""Paths and patches"""

from matplotlib.patches import PathPatch
from matplotlib.path import Path
from numpy import asarray, concatenate, ones
import seaborn as sns
from shapely.geometry import Point
from shapely.affinity import scale, rotate
import shapely.geometry as sg
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
import matplotlib.cm as cm
from matplotlib.colors import Normalize

######################################################
# translate shapely polygon to plt patch
######################################################

def PolygonPath(polygon):
    """Constructs a compound matplotlib path from a Shapely or GeoJSON-like
    geometric object"""
    #this = Polygon(polygon)
    this = polygon
    assert this.geom_type == 'Polygon'
    def coding(ob):
        # The codes will be all "LINETO" commands, except for "MOVETO"s at the
        # beginning of each subpath
        n = len(getattr(ob, 'coords', None) or ob)
        vals = ones(n, dtype=Path.code_type) * Path.LINETO
        vals[0] = Path.MOVETO
        return vals
    vertices = concatenate(
                    [asarray(this.exterior.coords)[:, :2]] 
                    + [asarray(r.coords)[:, :2] for r in this.interiors])
    codes = concatenate(
                [coding(this.exterior)] 
                + [coding(r) for r in this.interiors])
    return Path(vertices, codes)


def PolygonPatch(polygon, **kwargs):
    """Constructs a matplotlib patch from a geometric object
    
    The `polygon` may be a Shapely or GeoJSON-like object with or without holes.
    The `kwargs` are those supported by the matplotlib.patches.Polygon class
    constructor. Returns an instance of matplotlib.patches.PathPatch.

    Example (using Shapely Point and a matplotlib axes):

      >>> b = Point(0, 0).buffer(1.0)
      >>> patch = PolygonPatch(b, fc='blue', ec='blue', alpha=0.5)
      >>> axis.add_patch(patch)

    """
    
    return PathPatch(PolygonPath(polygon), **kwargs)

def addPolygonPatch(ax, polygon, text='', **kwargs):
    #print(polygon.geom_type)
    if polygon.geom_type == 'Polygon':
        ax.add_patch(PolygonPatch(polygon, **kwargs))
        xy = polygon.centroid.coords.xy
        ax.text( xy[0][0], xy[1][0], text, ha='center', va='center', fontsize=12)
    elif polygon.geom_type == 'MultiPolygon':
        #return
        for geom in polygon.geoms:
            xy = geom.centroid.coords.xy
            ax.text( xy[0][0], xy[1][0], text, ha='center', va='center', fontsize=12)
            #print(geom)
            ax.add_patch(PolygonPatch(geom, **kwargs))
            break
    


# given sliver specification (0/1 for each model), return shapely polygon

def get_sliver(circles_list, sliver_membership):
    sliver_membership = list (sliver_membership)
    # get a circle to start with
    for first_circle, member in zip (circles_list, sliver_membership):
        if member:
            break
    final_sliver = first_circle
    for circle, member in zip (circles_list, sliver_membership):
        #if (circle != first_circle) and final_sliver.intersects(circle):
        if member:
            final_sliver = final_sliver.intersection(circle)
        else:
            final_sliver = final_sliver.difference(circle)
    #print(sliver_membership, final_sliver)
    return final_sliver

    
        


# create shapely ellipse
def create_ellipse(x, y, w, h, inclination):
    p = Point(x,y)
    c = p.buffer(1)
    ellipse = scale(c, w, h, origin='centroid')
    ellipse = rotate(ellipse, inclination)
    return ellipse






def multimodel_venn(model_names,  cmap, data_dict, ax=None, vmin=None, vmax=None, annfmt = '.2g', cbar_label='', alpha=1, fontsize=18, ec='k', **kwargs):
    
    if ax is None:
        ax= plt.gca()
        
    ax.set_ylim(bottom=-0.2, top=1.5)
    ax.set_xlim(left=-1, right=2)
    ax.set_aspect('equal')
    ax.set_axis_off()

    radius = 1.1
    # create the circles with shapely
    circles_list = [
        create_ellipse(0.400, 0.400, 0.72, 0.3, 150.0),
        create_ellipse(0.540, 0.600, 0.72, 0.3, 140.0),
        create_ellipse(0.504, 0.600, 0.72, 0.3, 50.0),
        create_ellipse(0.644, 0.400, 0.72, 0.3, 40.0),
    ]
    # get a set of slivers
    slivers = {i:get_sliver(circles_list, i) for i in product([False,True], repeat=4) if i != (0,0,0,0)}
               
    # figure out colors/text per sliver
    #vmin = np.amin()
    if vmin is None:
        vmin = data_dict.min()
    if vmax is None:
        vmax = data_dict.max()
        
    cmap_norm = Normalize(vmin=vmin, vmax=vmax)
    

    for sliver, polygon in slivers.items():
        #print (sliver, c)
        #print (morder, sliver, data_dict[sliver])
        annotation=''
        c='white'
        if sliver in data_dict:
            c = cmap(cmap_norm(data_dict[sliver]))
            annotation = ("{:" + annfmt + "}").format(data_dict[sliver])
        addPolygonPatch(ax, polygon, fc=c, ec=ec,  alpha=alpha, text=annotation, **kwargs)
    
    ax.text(0.1, 1.25, model_names[1], fontsize=fontsize, ha="center", va='top')
    ax.text(-0.3, 0.88, model_names[0], fontsize=fontsize, ha="center", va="bottom")
    ax.text(1.28, 0.92, model_names[3], fontsize=fontsize, ha="center", va="bottom")
    ax.text(1.0, 1.28, model_names[2], fontsize=fontsize, ha="center", va="top")
        
    # for sliver in circles_list:
    #     #print (sliver, c)
    #     #ax.add_patch(PolygonPatch(sliver, fc=c, ec='k', alpha=0.3))
    #     # TODO
    #     addPolygonPatch(ax, sliver, ec='k', fill=None)

    plt.colorbar(cm.ScalarMappable(norm=cmap_norm, cmap=cmap), ax=ax, label=cbar_label)
    plt.title(cbar_label)
