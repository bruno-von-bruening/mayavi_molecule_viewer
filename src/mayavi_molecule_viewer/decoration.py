from . import *

atom_color_di=dict(
    #https://www.figma.com/colors/ : look for RGB value
    H   = (0.9,0.9,0.9) , 
    #
    Li  = (.7,.11,.11),
    Be  = (.44, .26, .08),
    B   = (1,.75,0.8),
    C   = (0,0,0), 
    N   = (0,0,1), 
    O   = (1,0,0), 
    F   =(0.0,1,0.25), 
    #
    Si  = (.28, 51, .71),
    P   = (0.53,0.2,.38),
    S   = (0.9,.82,.04),
    Cl  = (0,1,0), 
    #
    Br  = (0.6,0.3,0),
    default=(0,1,1), # I dont know about cyan to be used
)
atom_scale_di=dict(
    H   = 0.75 , 
    #
    Li  = 1,
    Be  = 1,
    B   = 1,
    C   = 1, 
    N   = 1, 
    O   = 1, 
    F   = 1,
    #
    Si  = 1.5,
    P   = 1.5,
    S   = 1.5, 
    Cl  = 1.5, 
    #
    Br  = 2.0, 
    default=1.,
)
@validate_call
def get_atom_scale(atom_type:str) -> float:
    if not atom_type in atom_scale_di.keys():
        logging.warning(f"Could not find atom scale for atom \'{atom_type}\'")
        return atom_scale_di['default']
    else:
        return atom_scale_di[atom_type]
@validate_call
def get_atom_color(atom_type:str) -> Tuple[float,float, float]:
    if not atom_type in atom_color_di.keys():
        logging.warning(f"Could not find atom color for atom \'{atom_type}\'")
        return atom_color_di['default']
    else:
        return atom_color_di[atom_type]

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

allowed_values=[1,1.5,2,3] # 1.5 is for aromatic
bond_thickness_mapper={
    1:1.,
    1.5: 1.33,
    2  : 1.66,
    3  : 2.5,
}
    #dict([  (v,v) for v in allowed_values])

@validate_call
def get_bond_thickness(bond_order: int|float) -> float:

    if not bond_order in bond_thickness_mapper.keys():
        logging.warning(f"Bond order {bond_order} is unkown")

    tube_factor=bond_thickness_mapper[bond_order]
    return tube_factor