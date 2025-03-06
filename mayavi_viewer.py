#!/usr/bin/env python

import os, sys
# In case this is Bruno's machine where this variable exists, add to path
if 'PY_MODULES' in os.environ.keys():
    sys.path.insert(1, os.environ['PY_MODULES'])

# Custom module
import mod_objects as m_obj

import numpy as  np
from mayavi import mlab
import mayavi.tools

# Global parameters
see_through_options=[ x.upper() for x in ['front','back','none']]

def make_colors(values):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap

    cmap_name="RdBu"
    average=np.average(values)
    values_rel_to_average=values-average
    max_val, min_val=[max(values_rel_to_average), min(values_rel_to_average) ]

    abs_val=max( abs(max_val), abs(min_val))
    values_shifted=values_rel_to_average+abs_val


    cmap = cm.get_cmap(cmap_name, 2*abs_val)

    colors=cmap(values_shifted)

    return cm.get_cmap(cmap_name, 12)

def get_scale(values):
    average=np.average(values)
    values_rel_to_average=values-average
    max_val, min_val=[max(values_rel_to_average), min(values_rel_to_average) ]

    abs_val=max( abs(max_val), abs(min_val))

    return average-abs_val, average+abs_val

# Plot molecules
def plot_mol(geometry_object=None, coords=None, atom_types=None, dont_show_indices=False ):
    atom_color_di={ 'O':(1,0,0), 'H':(0.9,0.9,0.9) , 'N':(0,0,1), 'C':(0,0,0), 'Cl':(0,1,0), 'Br':(0.6,0.3,0)}
    atom_scale_di={ 'O':1, 'H':0.75 , 'N':1, 'C':1, 'Cl':1.5, 'Br':2.0}
    if isinstance(geometry_object, m_obj.geometry):
        coords=geometry_object.coordinates
        atom_types=geometry_object.atom_types
        adj_mat=geometry_object.adjacency_matrix
    else:
        quit()
    
    # Plot the atoms
    for i,z in enumerate( zip(coords,atom_types) ):
        coor, tag= z
        at_col=atom_color_di[tag]    
        at_scl=atom_scale_di[tag]
        atom = mlab.points3d( *list(coor), scale_factor=at_scl, resolution=50, color=at_col, scale_mode='none')
        if not dont_show_indices:
            label = mlab.text3d(*list(coor) , f"{i}", color=at_col)

    # Plot the connections
    if type(adj_mat)==type(None):
        plain_xyz= np.array(coords).T
        mlab.plot3d(*plain_xyz,
                    tube_radius=0.2, color=(0.5,0.5,0.5), tube_sides=20)
    else:
        for i,con in enumerate(adj_mat):
            for j,target in enumerate(con):
                plot=False
                tube_radius=0.2
                if i==j or j>i: # avoid plotting to itself or doubling
                    pass
                elif con[j] in [0,1,2,3]:
                    if con[j]==0:
                        pass
                    elif con[j] in [1,2,3]:
                        plot=True
                        if con[j]==2:
                            tube_radius=tube_radius*2
                        elif con[j]==2:
                            tube_radius=tube_radius*3
                else:
                    print(con[j])
                    raise Exception

                if plot:
                    origin=coords[i]
                    target=coords[j]
                    plain_xyz=np.array([ origin, target]).T
                    mlab.plot3d(*plain_xyz,
                                tube_radius=tube_radius, color=(0.5,0.5,0.5), tube_sides=20)
        # Have a look into this to fuse it together
        # https://stackoverflow.com/questions/54144002/drawing-disconnected-lines-in-mayavi-calling-mlab-plot3d-once
    return atom 
    
def plot_surface(xyz=None, esp=None, triangles=None, opacity=1, see_through=None, plot_colorbar=True):
    if not isinstance(see_through,type(None)):
        if not see_through.upper() in see_through_options:
            raise Exception(f"Provided option for see through is not recognized: {see_through}, options are \
                    {see_through_options}")
    def manage_color(esp, transparent=False, opacity=1):
        cmap='RdBu'
        vmin, vmax=get_scale(esp)
        from functools import partial
        my_mesh=partial( mlab.triangular_mesh, opacity=opacity, colormap=cmap, transparent=transparent, vmax=vmax, vmin=vmin)
        return my_mesh

    def make_color_bar(surf):
        lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()
        surf.module_manager.scalar_lut_manager.lut.table = lut
        title='Energy in kj/mol\n(centered around average)'
        cb=mayavi.mlab.colorbar(object=surf, title=title, orientation='vertical', nb_labels=None, nb_colors=None,
                label_fmt='%.1f')
        cb.scalar_bar.unconstrained_font_size = True
        cb.label_text_property.font_size=24
        return surf

    def set_see_through(surf, see_through):
        if isinstance(see_through, type(None)):
            pass
        elif see_through.lower() in ['back']:
            surf.actor.property.frontface_culling = True
        elif see_through.lower() in ['front']:
            surf.actor.property.backface_culling = True
        elif see_through.lower() in ['none']:
            pass
        else:
            raise Exception(f"Received key for see_through option \'{see_through}\' but no case implemented")
        return surf 

    if isinstance(xyz,list):
        xyz=np.array(xyz)
    
    # Set the colors for the mesh
    make_mesh=manage_color(esp, transparent=False, opacity=opacity)
    # plot the surface
    surf=make_mesh( *xyz.T, triangles, representation='surface', scalars=esp )

    # colorbar (for heatmap numeric values)
    if plot_colorbar:
        surf=make_color_bar(surf)
    # set see through
    surf = set_see_through(surf, see_through)

    return surf

#molecule=mlab.figure(figure='molecule')
def plot_esp_surface(esp=None, esp_ref=None, geom=None, opacity=1, see_through=None,
        dont_show_indices=False, title=None):
    """ Provide an esp map (as .map file) and geometry (as .mom file).
    This function will generate a mayavi plot!"""

    def check_coordiantes(array_0, array_1):
        comp_thres=1.e-8
        assert array_0.shape==array_1.shape
        assert max( (array_1-array_0).reshape(-1) )< comp_thres

    def read_map(the_list):
        esp_map=None
        for item in the_list:
            the_map=m_obj.grid(item)
            if isinstance(esp_map, type(None)):
                esp_map=the_map
            else:
                check_coordiantes(esp_map.grid_points, the_map.grid_points)
                esp_map.values+=the_map.values
        return esp_map
    
    def make_figure(the_title=None):
        # Plot settings
        #   Colors
        background_color=(1,1,1)
        foreground_color=(0,0,0)
        #   Mayavi figure
        figure = mlab.figure(the_title, bgcolor=background_color, fgcolor=foreground_color, size=(500, 500))
        mlab.clf()
        return figure

    
    # Input checks
    if isinstance(esp, type(None)):
        raise Exception(f"Provided esp")
    elif isinstance(esp, list):
        esp_map=read_map(esp)
    else:
        raise Exception(f"Unkown format for \'esp\' argument: {type(esp_map)}")
    
    #if isinstance(geom, type(None)):
    #    raise Exception(f"Provide geometry")

    if isinstance(esp_ref, type(None)):
        esp_values=esp_map.values
        pass
    else:
        if isinstance(esp, type(None)):
            raise Exception(f"Provided reference map but not actual map!")
        if isinstance(esp, list):
            esp_ref_map=read_map(esp_ref)
            check_coordiantes(esp_map.grid_points, esp_ref_map.grid_points)
        else:
            raise Exception(f"Unkown type of \'esp_ref\': {type(esp_ref)}")    
        esp_values=esp_map.values-esp_ref_map.values

    figure=make_figure(title)

    # Get the map data
    xyz,esp, triangles=[esp_map.grid_points, esp_values, esp_map.triangles]
    # Get the atom data
    if not isinstance(geom, type(None)):
        multipoles=m_obj.multipoles(geom)
        geom = multipoles.return_geom(connectivity=True)
        geom.connec()
        molecule=plot_mol(geometry_object=geom, dont_show_indices=dont_show_indices)

    # Plot
    surface=plot_surface(xyz=xyz, esp=esp, triangles=triangles, opacity=opacity, see_through=see_through)
    return [surface] 




if __name__ == '__main__':

    # Setup parser
    indent=4*' '
    description=f"This routine takes a molecule and a map and plots these object using the mayavi library"
    epilog=f"Usage:\n{indent} {__file__} --map_ref <ESP-WFN*.map> --map <[ESP-DMP*.map]> --geom <*.mom> --opacity 0.75 --see_through back --title <my_title>"
    import argparse as ap; par=ap.ArgumentParser(description=description, epilog=epilog, formatter_class=ap.RawDescriptionHelpFormatter)
    # Define Arguments
    adar=par.add_argument
    adar(
        '--map', nargs='*', help='ESP MAP: esp values to be plotted. Expects .map extension'
    )
    adar(
        '--map_ref', nargs='*', help='ESP REF MAP: esp reference values -> plot difference between ESP MAP and ESP REF MAP. Expects .map extension. (coordinates between ESP MAP and ESP MAP DEF neeed to agree!)'
    )
    adar(
        '--geom', help='file with molecular geometry (for molecule visualization)'
    )
    adar(
        '--opacity', help='Opacity of surface', type=float, default=1.
    )
    adar(
        '--see_through', help='Allows one-side see through', type=str, choices=see_through_options,
            default='back',
    )
    adar(
        '-no_indices', help='Do_not show indices', action='store_true', default=False
    )
    adar(
        '--title', '-t', help=f"Title of the plot (at the moment only displayed as window title, useful for distinguish plots)"
    )
    # PARSE ARGUMENTS
    args=par.parse_args()
    map_fi=args.map
    map_ref_fi=args.map_ref
    dont_show_indices=args.no_indices
    mom_fi=args.geom
    opacity=args.opacity
    see_through=args.see_through
    title=args.title

    plot_esp_surface(esp=map_fi, esp_ref=map_ref_fi, geom=mom_fi, opacity=opacity, see_through=see_through,
            dont_show_indices=dont_show_indices, title=title)
    mlab.show()
    # mlab keyboard actions https://docs.enthought.com/mayavi/mayavi/auto/mlab_figure.html

