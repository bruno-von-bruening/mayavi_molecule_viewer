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

# Plot molecules
def plot_mol(geometry_object=None, coords=None, atom_types=None ):
    atom_color_di={ 'O':(1,0,0), 'H':(0.9,0.9,0.9) , 'N':(0,0,1), 'C':(0,0,0)}
    atom_scale_di={ 'O':1, 'H':0.75 , 'N':1, 'C':1}
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
    if not type(see_through)==type(None):
        if not see_through.upper() in see_through_options:
            raise Exception(f"Provided option for see through is not recognized: {see_through}, options are \
                    {see_through_options}")
    xyz=[ list(x) for x in xyz ]
    esp=[ x for x in esp ]
    xyz*=2
    esp*=2
    xyz=np.array(xyz)
    esp=np.array(esp)
    #mlab.triangular_mesh(0.99*xyz.T[0], 0.99*xyz.T[1], 0.99*xyz.T[2], triangles[:], opacity=0.75, representation='surface',
    #        color=(1,1,1),)
    surf=mlab.triangular_mesh(xyz.T[0], xyz.T[1], xyz.T[2], triangles[:], opacity=opacity, representation='surface', scalars=esp,
            colormap='RdBu', transparent=False)
    lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()
    for i,entry in enumerate(lut):
        minimum=min( entry[:3] )
    #    lut[i, -1] = 255-minimum*0.9
    surf.module_manager.scalar_lut_manager.lut.table = lut

    if type(see_through)==type(None):
        pass
    elif see_through.lower() in ['back']:
        surf.actor.property.frontface_culling = True
    elif see_through.lower() in ['front']:
        surf.actor.property.backface_culling = True
    elif see_through.lower() in ['none']:
        pass
    else:
        raise Exception(f"Received key for see_through option \'{see_through}\' but no case implemented")
    #surf=mayavi.tools.pipeline.triangular_mesh_source(xyz.T[0], xyz.T[1], xyz.T[2], triangles[:],scalars=esp,
    #name='surface')
    if plot_colorbar:
        mayavi.mlab.colorbar(object=surf, title=None, orientation='vertical', nb_labels=None, nb_colors=None,
                label_fmt='%.1f')
    return surf
    #mlab.triangular_mesh(xyz.T[0], xyz.T[1], xyz.T[2], triangles[:], opacity=0.95, representation='wireframe', scalars=esp)
    #mlab.triangular_mesh(xyz.T[0], xyz.T[1], xyz.T[2], triangles[0::2])

#molecule=mlab.figure(figure='molecule')
def plot_esp_surface(esp=None, geom=None, opacity=1, see_through=None):
    """ Provide an esp map (as .map file) and geometry (as .mom file).
    This function will generate a mayavi plot!"""
    # Input checks
    if type(esp)==None:
        raise Exception(f"Provided esp")
    elif type(geom)==None:
        raise Exception(f"Provide geometry")

    # Plot settings
    #   Colors
    background_color=(1,1,1)
    foreground_color=(0,0,0)
    #   Mayavi figure
    figure = mlab.figure(1, bgcolor=background_color, fgcolor=foreground_color, size=(350, 350))
    mlab.clf()

    # Get the map data
    esp_map=m_obj.grid(esp)
    xyz,esp, triangles=[esp_map.grid_points, esp_map.values, esp_map.triangles]


    # Get the atom data
    multipoles=m_obj.multipoles(geom)
    geom = multipoles.return_geom(connectivity=True)
    geom.connec()

    # Plot
    surface=plot_surface(xyz=xyz, esp=esp, triangles=triangles, opacity=opacity, see_through=see_through)
    molecule=plot_mol(geometry_object=geom)#,atom_coords, atom_types)
    #mlab.show()
    return [surface] 




if __name__ == '__main__':
    import argparse as ap; par=ap.ArgumentParser()
    par.add_argument('-map', help='.map file')
    par.add_argument('-geom', help='file with molecular geometry (for molecule visualization)')
    par.add_argument('-opacity', help='Opacity of surface', type=float, default=1.)
    par.add_argument('-see_through', help='Allows one-side see through', type=str, choices=see_through_options,
            default=None)
    args=par.parse_args()
    map_fi=args.map
    mom_fi=args.geom
    opacity=args.opacity
    see_through=args.see_through

    plot_esp_surface(esp=map_fi, geom=mom_fi, opacity=opacity, see_through=see_through)
    mlab.show()
