#!/usr/bin/env python
from mayavi_molecule_viewer.main import plot_esp_surface, see_through_options
from mayavi import mlab
import mayavi
import os
import yaml
from mayavi_molecule_viewer.config import mayavi_config
from pydantic import Field,BeforeValidator,validate_call; val_call=validate_call(config=dict(arbitrary_types_allowed=True))
from typing import Annotated

def check_file(input):
    if not isinstance(input,str):
        raise Exception(f"not a string")
    else:
        return os.path.isfile(input)
pdtc_file=Annotated[str, BeforeValidator(check_file)]
def main(
    # options
    config: mayavi_config|str=mayavi_config(),
    print_config:bool=True,
    # data
    map_file:pdtc_file|None=None,
    geometry_file:pdtc_file|None=None,
    reference_map_file:pdtc_file|None=None,
    dump_file:pdtc_file|None=None,
    opacity: float|None=Field(default=None,ge=0,le=1),
    see_through: bool|None=None,
    no_indices: bool|None=None,
    title:str|None=None,
):

    def set_config(config):
        if config is  None:
            config=mayavi_config()
        else:
            assert os.path.isfile(config)
            with open(config,'r') as rd:
                config=yaml.safe_load(rd)
                config=mayavi_config(**config)
        return config
    def drop_default_config(config):
        config_file='default_config.yaml'
        with open(config_file,'w') as wr:
            yaml.safe_dump(config.model_dump(), wr)
        print(f"Printed default config file to {config_file}")

    # Either use default config or load config from provided file
    config=set_config(config)
    
    # Either print default config or plot
    if print_config:
        drop_default_config(config)
    else:
        from functools import partial
        the_func=partial(plot_esp_surface, esp=map_file, esp_ref=reference_map_file, geom=geometry_file, config=config)
            #opacity=opacity, see_through=see_through,
            #    dont_show_indices=dont_show_indices, title=title)

        if dump_file is not None:
            # do not open mayavi window
            mlab.options.offscreen = True
            the_func()
            mlab.savefig(dump_file, size=config.figure.output_size)
        else:
            the_func()
            mlab.show()

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
        '--no_indices', help='Do_not show indices', action='store_true', default=False
    )
    adar(
        '--title', '-t', help=f"Title of the plot (at the moment only displayed as window title, useful for distinguish plots)"
    )
    adar(
        '--config', '-c', type=str, help=f"Name of config_file"
    )
    adar(
        '--print_config','--pc', action='store_true', help=f"Prints example of config file with default values"
    )
    adar(
        '--dump', type=str, help=f"File name under which to dump file"
    )
    # PARSE ARGUMENTS
    mapper=dict(
        geom='geometry_file',
        dump='dump_file',
        map_ref='reference_map_file',
        map='map_file'
    )
    args_di=dict([ 
        (mapper[k],v)  if k in mapper.keys() else (k,v)
        for k,v in 
        vars(par.parse_args()).items()
    ])
    #args=vars(args)
    #map_fi=args.map
    #map_ref_fi=args.map_ref
    #dont_show_indices=args.no_indices
    #mom_fi=args.geom
    #opacity=args.opacity
    #see_through=args.see_through
    #title=args.title
    #print_config=args.print_config
    #config=args.config
    #dump_path=args.dump


    main(
        **args_di
    )

    # mlab keyboard actions https://docs.enthought.com/mayavi/mayavi/auto/mlab_figure.html

