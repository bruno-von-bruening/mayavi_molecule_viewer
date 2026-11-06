from . import *
from .config import mayavi_config
import json


def manage_color(esp, config):

    def get_min_max(esp, config):
        vmin, vmax=config.scalars.get_scale(esp)
        provided_scale=config.scalars.range
        if provided_scale is None:
            return vmin, vmax
        else:
            provided_min, provided_max=provided_scale
            if not provided_min is None:
                vmin=provided_min
            if not provided_max is None:
                vmax=provided_max
            return vmin, vmax
        
    from functools import partial
    vmin, vmax = get_min_max(esp, config)
    my_mesh=partial( 
                    mlab.triangular_mesh, 
                    opacity=config.surface.opacity, 
                    colormap=config.surface.color_map, 
                    #transparent=transparent, 
                    vmax=vmax, 
                    vmin=vmin
    )
    return my_mesh

@validate_call
def make_color_bar(surf, config: mayavi_config=mayavi_config()):

    
    lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()
    surf.module_manager.scalar_lut_manager.lut.table = lut

    cb=mlab.colorbar(
        object=surf, 
        title           =config.scalars.colorbar_title, 
        orientation     =config.scalars.orientation, 
        nb_labels       =None, 
        nb_colors       =None,
        label_fmt       =config.scalars.label_format,
    )
    cb.scalar_bar.unconstrained_font_size = True
    cb.label_text_property.orientation=1#'reverse'
    cb.label_text_property.font_size=24
    #if not config.scalars.labels is None:
    #    cb.set_ticklabels(config.scalars.labels)
    return surf

def set_see_through(surf, config):
    see_through=config.surface.see_through
    if see_through.lower() in ['back']:
        surf.actor.property.frontface_culling = True
    elif see_through.lower() in ['front']:
        surf.actor.property.backface_culling = True
    elif see_through.lower() in ['none']:
        pass
    else:
        raise Exception(f"Received key for see_through option \'{see_through}\' but no case implemented")
    return surf 
@validate_call
def plot_surface(
        xyz: List[Tuple[float,float,float]], esp: List[float], triangles: List[Tuple[int,int,int]], 
        config: mayavi_config=mayavi_config(),
):
    
    xyz=np.array(xyz)

    # Set the colors for the mesh
    make_mesh=manage_color(esp, config)

    # plot the surface
    surf=make_mesh( *xyz.T, triangles, representation='surface', scalars=esp )

    # colorbar (for heatmap numeric values)
    if config.scalars.show_colorbar:
        surf=make_color_bar(surf, config)

    if config.scalars.drop_colorbar:
        lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()
        with open('colorbar.json', 'w') as wr:
            json.dump(lut.tolist(), wr)
    
    # set see through
    surf = set_see_through(surf, config)

    return surf