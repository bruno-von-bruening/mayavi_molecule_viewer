# mayavi_molecule_viewer

## Installation
First create a conda environment with the provided file:
```
conda env create --file=install/mayavi_mol-view_conda_env.yaml
```
at the moment the package versions are quite tightly defined and if one want to install newer version that may require
some trial-and-error with these settings.

Then install mayavi_molecule_viewer as a local package:
```bash
conda activate mayavi
python -m pip install --upgrade build
python -m build
python -m pip install .
```

## Running the Code
The code can be found under `mayavi_viewer.py` and references an object that can be found under `mod_objects.py`.
The code can be invoked either by launching the supplied jupyter notebook `interactive_jupyter_plot.ipynb` or by direct exectuion (currently not possible since it links to module not included in this repository).
```
./mayavi_viewer.py -geom moment_files/conf_1_ISA-GRID.mom -map map_files/conf_1_ISA-GRID_dip_iso-p001.map 
```
Latter will open the mayavi GUI. This example should show an inversion in one way see through (wron orientation).
The effect can be fixed by opening the scene settings in the GUI (mayavi logo on left-upper-hand) and setting the check
under 'Mayavi Scene 1 -> TriangularMeshSource -> PolyDataNormals -> Auto orient normals'. (One may activate this
settings when launching the GUI from the python script but I have not figure out yet how...)

Note that example files have been provided in the `moment_files/` and `map_files/` directories (a peptide and the water
molecule) respectively. Both file types are derived from the CamCASP suite developed by Alston Misquitta.
[https://wiki.ph.qmul.ac.uk/ccmmp/AJMPublic/camcasp]

The user may customize the processing of arbitrary input files by modification of the classes in the provided
`mod_objects.py` file.

## Changing parameters
*(Command line options are currently disabled but could overwrite the config)*

The configuration of the plot can be controlled by a config file. An example with the defaults can be obtained by
running:
```bash
   bin/mayavi_viewer.py --print_config
```
Then the file can be edited:
```bash
   mv <auto generated config> my_config.yaml 
   vim my_config.yaml
```
And the mayavi plot be rerun by providing the file:
```bash
   bin/mayavi_viewer.py -geom moment_files/conf_1_ISA-GRID.mom -map map_files/conf_1_ISA-GRID_dip_iso-p001.map --config
   my_config.yaml
```

