from . import *
color=Tuple[float,float,float]

class mayavi_view(BaseModel):
    azimuthal       : float     = Field(1,ge=0., le=360.)
    elevation       : float     = Field(1, ge=0, le=180.)
    distance        : float     = Field(50,ge=0)
    focalpoint      : Tuple[float,float,float]=(0,0,0)
    roll            : float|None     = None
    reset_roll      : bool      = False
    figure          : None      = None
    parallel_projection: bool=True # parallel projection (as viewer would have infinite distance from object)
    parallel_scale  : float|None=None 
    light_intensity : List[float]|float = Field(1.,ge=0., le=1.)

class mavayi_scalars(BaseModel):
    show_colorbar   : bool=True
    drop_colorbar   : bool=False
    colorbar_title  : str='Energy in kj/mol\n(centered around average)'
    orientation     : Literal['vertical','horizontal']='vertical'
    decimal_places  : int=Field(1,ge=0)
    range           : Tuple[float|None,float|None]|None=None
    labels          : List[float|int|str]|None=None
    
class mayavi_figure(BaseModel):
    title               : str|None=None
    background_color    : color = (1,1,1)
    foreground_color    : color = (0,0,0)
    size                : Tuple[float,float]=(500,500)
    output_size         : Tuple[float,float]=(500,500)

class mayavi_molecule(BaseModel):
    show_bonds: bool=True
    show_indices: bool=True
    base_bond_thickness: float=Field(0.2, gt=0)
    tube_sides: int=20
    bond_color: color=(0.5,0.5,.5)

class mayavi_surface(BaseModel):
    color_map       : str='RdBu'
    opacity         : float=Field(.75,ge=0,le=1)
    see_through     : Literal['FRONT','BACK','NONE']='BACK'


class mayavi_view(mayavi_view):
    def set_view(self):
        fig=mlab.gcf()
        fig.scene.parallel_projection = self.parallel_projection
        if not self.parallel_scale is None:
            fig.scene.camera.parallel_scale = self.parallel_scale
        def set_lights():
            """ Set light according to provided light intensity key (could also set elevation (angle to z-axis) and azimuthal (lcoation on xy plane))
            I got this from following forum post: https://github.com/enthought/mayavi/issues/333
            Apparently there are 3 lights by standard
            """
            camera_lights = fig.scene.light_manager.lights
            if isinstance(self.light_intensity, float):
                for cam_light in camera_lights:
                    cam_light.intensity = self.light_intensity
            elif isinstance(self.light_intensity, (list, tuple)):
                assert len(camera_lights)==len(self.light_intensity), f"Provided an tuple of lkist for light intensity corresponding length unequal number of lights"
                f"({len(self.light_intensity)}!={len(camera_lights)})"
                for cam_light, intensi in zip(camera_lights, self.light_intensity):
                    cam_light.intensity = intensi 
            #camera_light0.elevation = elevation_light0
            #camera_light0.azimuth = azimuth_light0
        set_lights()
        mlab.view(
            azimuth=self.azimuthal, 
            elevation=self.elevation, 
            distance=self.distance, 
            focalpoint=self.focalpoint,
            roll=self.roll, 
            reset_roll=self.reset_roll, 
            figure=self.figure
        )
class mayavi_figure(mayavi_figure):
    def make_figure(self):
        figure = mlab.figure(
            self.title, 
            bgcolor=self.background_color, 
            fgcolor=self.foreground_color, 
            size=self.size,
        )
        mlab.clf()
        
        return figure
class mayavi_molecule(mayavi_molecule):
    def bond_thickness(self, bond_order):
        from .decoration import get_bond_thickness
        return get_bond_thickness(bond_order)*self.base_bond_thickness
class mayavi_surface(mayavi_surface):
    @field_validator('see_through', mode='before')
    @classmethod
    def validate_command(cls, v: str|None) -> str:
        if v is None:
            v='none'
        return v.upper()
class mavayi_scalars(mavayi_scalars):
    @property
    def label_format(self):
        return "%.{}f".format(self.decimal_places)

    @validate_call
    def get_scale(self, values: List[float]):
        average=np.average(values)
        values_rel_to_average=values-average
        max_val, min_val=[max(values_rel_to_average), min(values_rel_to_average) ]

        abs_val=max( abs(max_val), abs(min_val))

        return average-abs_val, average+abs_val

class mayavi_config(BaseModel):
    view        : mayavi_view=mayavi_view()
    scalars     : mavayi_scalars=mavayi_scalars()
    figure      : mayavi_figure=mayavi_figure()
    molecule    : mayavi_molecule=mayavi_molecule()
    surface     : mayavi_surface=mayavi_surface()


    def __init__(self, *args, **kwargs):
        #   # Add empty directory for unpresent subclasses
        #   keys=['view','scalars','layout','molecule','figure']
        #   for key in keys:
        #       if not key in kwargs.keys():
        #           kwargs.update({key:{}})
        #   # init
        super().__init__(*args, **kwargs)
