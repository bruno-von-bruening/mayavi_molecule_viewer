from . import *

from .decoration import atom_color_di, atom_scale_di, make_colors
from .config import mayavi_config



import mayavi.tools

# Global parameters
see_through_options=[ x.upper() for x in ['front','back','none']]

@validate_call
def check_coordiantes(array_0:np.ndarray, array_1:np.ndarray):
    comp_thres=1.e-8
    assert array_0.shape==array_1.shape
    assert max( (array_1-array_0).reshape(-1) )< comp_thres

def read_esp(map_input):
    # Input checks
    def read_map(the_list):
        esp_map=None
        for item in the_list:
            the_map=m_obj.grid(item)
            if isinstance(esp_map, type(None)):
                esp_map=the_map
            else:
                check_coordiantes(np.array(esp_map.grid_points), np.array(the_map.grid_points))
                esp_map.values+= np.array(esp_map.values)+np.array(the_map.values)
        return esp_map
    if isinstance(map_input, type(None)):
        raise Exception(f"Provided invalid input for esp_map {map_input}")
    elif isinstance(map_input, str):
        map_esp=read_map([map_input])
    elif isinstance(map_input, list):
        map_esp=read_map(map_input)
    else:
        raise Exception(f"Unkown format for \'map_input\' argument: {type(map_input)}")
    
    return map_esp
def read_geom(geom_input):
    if geom_input.endswith('.xyz'):
        geom=m_obj.geometry(geom_input)
    elif geom_input.endswith('.mom'):
        multipoles=m_obj.multipoles(geom_input)
        geom = multipoles.return_geom(connectivity=True)
    elif geom_input.lower().endswith('.fchk'):
        from qcp_objects.objects.properties import geometry
        geom=geometry(geom_input)
    else:
        raise NotImplementedError(f"Cannot process file extension of {geom_input}")
    try:
        geom.connec()
    except Exception as ex:
        print(f"Warning could not auto generated connectivity: {str(ex)}")
    return geom


@validate_call
def plot_esp_surface(
        esp:List[str]|str|None=None, 
        esp_ref=None, 
        geom=None, 
        config:mayavi_config=mayavi_config(),
        #opacity=1, see_through=None,
        #dont_show_indices=False, title=None,
):
    """ Provide an esp map (as .map file) and geometry (as .mom file).
    This function will generate a mayavi plot!"""



    # read the esp map
    esp_map=read_esp(esp)
    if not esp_ref is None:
        esp_ref_map=read_esp(esp_ref)
        check_coordiantes(np.array(esp_map.grid_points), np.array(esp_ref_map.grid_points))
        esp_values=np.array(esp_map.values)-np.array(esp_ref_map.values)
    else:
        esp_values=esp_map.values


    figure=config.figure.make_figure()

    if not isinstance(geom, type(None)):
        from .molecule import plot_mol
        geom=read_geom(geom)
        molecule=plot_mol(geometry_object=geom, config=config)

    # Plot
    from .surface import plot_surface
    xyz,esp, triangles=[esp_map.grid_points, esp_values, esp_map.triangles]
    surface=plot_surface(
            xyz=xyz, esp=esp, triangles=triangles, 
            config=config,
    )
    #mlab.show()

    config.view.set_view()
    return [surface] 
