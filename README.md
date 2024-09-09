# mayavi_molecule_viewer

## Installation
The python scripts require a conda environment which can be installed from a `yaml` file:
```
conda env create --file=install/mayavi_conda-env_essential-export.yml
```
This install has been tested on Sep-09-2024.
Note that a full export of the developer's environment can be found under `install/mayavi_conda-env_full-export.yml`

## Running the Code
The code can be found under `mayavi_viewer.py` and references an object that can be found under `mod_objects.py`.
The code can be invoked either by launching the supplied jupyter notebook `interactive_jupyter_plot.ipynb` or by direct exectuion.
```
./mayavi_viewer.py -geom moment_files/conf_1_ISA-GRID.mom -map map_files/conf_1_ISA-GRID_dip_iso-p001.map -opacity=1
```

Note that example files have been provided in the `moment_files/` and `map_files/` directories (a peptide and the water
molecule) respectively. Both file types are derived from the CamCASP suite developed by Alston Misquitta.
[https://wiki.ph.qmul.ac.uk/ccmmp/AJMPublic/camcasp]

The user may customize the processing of arbitrary input files by modification of the classes in the provided
`mod_objects.py` file.