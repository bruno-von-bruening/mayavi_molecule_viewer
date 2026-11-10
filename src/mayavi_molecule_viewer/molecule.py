from . import *

from .decoration import get_atom_color, get_atom_scale, get_bond_thickness
from .config import mayavi_config

# Plot molecules
@validate_call
def plot_mol(geometry_object:m_obj.geometry, coords=None, atom_types=None, 
    config: mayavi_config=mayavi_config()
):
    
    geometry_object.units.LENGTH='BOHR'
    coords=geometry_object.coordinates
    atom_types=geometry_object.atom_types
    adj_mat=geometry_object.adjacency_matrix
    
    # Plot the atoms
    for i,z in enumerate( zip(coords,atom_types) ):
        coor, tag= z
        at_col=get_atom_color(tag)
        at_scl=get_atom_scale(tag)
        atom = mlab.points3d( *list(coor), scale_factor=at_scl, resolution=50, color=at_col, scale_mode='none')

        # Indices
        if config.molecule.show_indices:
            label = mlab.text3d(*list(coor) , f"{i}", color=at_col)

    # Plot the connections
    if type(adj_mat)!=type(None) and config.molecule.show_bonds:
    #    plain_xyz= np.array(coords).T
    #    mlab.plot3d(*plain_xyz,
    #                tube_radius=0.2, color=(0.5,0.5,0.5), tube_sides=20)
    #else:
        for i,con in enumerate(adj_mat):
            for j,bond_order in enumerate(con):

                if (i==j or j>i): # avoid plotting to itself or doubling
                    continue

                if bond_order==0:
                    continue

                # make the vector of the bond
                origin=coords[i]
                target=coords[j]
                plain_xyz=np.array([ origin, target]).T

                # plot the tube
                mlab.plot3d(*plain_xyz,
                            tube_radius=config.molecule.bond_thickness(bond_order), 
                            color=config.molecule.bond_color,
                            tube_sides=config.molecule.tube_sides)


        # Have a look into this to fuse it together
        # https://stackoverflow.com/questions/54144002/drawing-disconnected-lines-in-mayavi-calling-mlab-plot3d-once
    return atom 