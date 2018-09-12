#PyBasin: basin evolution, heat flow and thermochronology


## Downloading PyBasin

* click the download link on the right for a zip file of the source code
* or clone the repository




## Installing & running PyBasin

* Make sure you have installed the necessary python modules (see list below). The easiest is to use a python distribution that already includes all the necessary modules, like Anaconda (https://www.continuum.io/downloads), pythonxy (https://code.google.com/p/pythonxy/) or canopy (https://www.enthought.com/products/canopy/).
* Download and unzip the PyBasin source code
* Or clone the repository to keep your version up to date with the PyBasin version on the bitbucket. See "Cloning a Git repository" on this webpage for more info on cloning a repository: https://confluence.atlassian.com/bitbucket/clone-a-repository-223217891.html  
* Navigate to the PyBasin directory and execute PyBasin by executing the following command from the command line:

````sh
python pybasin.py
````	

* The model will now run with the default input dataset from the Roer Valley Graben. Check the model output in the directory ``model_output/example_dataset``
* Optionally: If you want to model apatite fission track data and you want this to run relatively fast, you will have to compile the fission track annealing module first. The fission track annealing module was also written in Fortran instead of Python to reduce computation time. Not ethat the pure Python version is still available, but is relatively slow. Compile the Fortran file by navigating to the subdirectory `pybasin/lib``, opening a terminal and running the following command:

````sh
f2py -c calculate_reduced_AFT_lengths.f90 -m calculate_reduced_AFT_lengths
````	

* If you want to use the Fortran version of the fission track annealing library you may need to install a fortran compiler first. For linux operating systems this should be included in your distribution, for Mac OSX you can find installers for Gfortran here: https://gcc.gnu.org/wiki/GFortranBinariesMacOS. I have not tested compiling Fortran code on windows systems. Note that you do not need to compile anything for modeling vitrinite reflectance and apatite (U-Th)/He data

## Dependencies

PyBasin requires the following Python packages to be installed:

- Python2.x: http://www.python.org/
- Numpy: http://www.numpy.org/
- Pandas: https://pandas.pydata.org/
- Scipy: https://www.scipy.org/
- Matplotlib: https://matplotlib.org/


## Running single or multiple models:

The command for running pybasin models is:

````sh
python pybasin.py input_directory -w well1,well2,well3
````

``input_directory`` is the directory that contains all input files. For the example dataset this should be ``example_dataset``. With the command line option -w you can specify which wells to run. This can either be a single well or a list of well separated by a comma. If you do not specify an input directory, PyBasin will use the default input directory defined in the file ``default_folder.txt``. If you do not specify which well to run at the command line PyBasin will look for a list of wells in the file ``model_scenarios.py`` in your input directory.

	

## Model input 

* Adjust the default location of the model input files by changing the file ``default_input_directory.txt``. The directory name specified in this file should be a subdirectory of the directory ``pybasin/input_data``
* Adjust the text files in the directory to adjust the well stratigraphy, lithology and thermal properties of the rocks that you are trying to model
* Adjust the model parameters in the file ``pybasin_params.py``. Use a text editor or preferably a dedicated Python editor like Spyder, PyCharm (https://www.jetbrains.com/pycharm/) or Enthought Canopy (https://www.enthought.com/products/canopy/) to adjust the model parameters.
* If you want to run multiple parameter sets and/or wells in one go, adjust the ``model_scenarios.py`` file and change the parameters you wish to change. For instance, changing from ``exhumation_magnitudes = [None]`` to ``exhumation_magnitudes = [500.0, 1000.0, 2000.0]`` will tell PyBasin to execute three a single model runs in one go with exhumation magnitudes of 500, 1000 and 2000 m, respectively. For automatically executing model runs for multiple wells change the ``wells`` list in the same file to include several wells.

## Model output:

* The model generates a series of .csv files that show the stratigraphy, burial depths and formation thicknesses over time. These files are stored in ``model_output/model_directory/burial_history_csv_files/``, where ``model_directory`` is the name of the input and output directory that you have specified in the file ``default_input_folder.txt``
* If ``make_model_data_fig = True`` in the file PyBasin_params.py, the model script will generates a single figure for each model run showing burial depth, temperature and modeled vitrinite relfectance, apatite fission track age and apatite (U-Th)/He data. THis figure is saved in the directory ``model_output/model_directory/model_data_fig/``
* If ``save_model_run_data = True`` in the file PyBasin_params.py, the model script will store all model simulation data for each model run in a datafile that uses the python pickle module. This file can be read later using a additional python script (still in the works...). The location where these files should be saved is specified in the ``PyBasin_params.py`` in line 19: ``datafile_output_dir = '../../heavy_data/PyBasin_MB'``
* You can also make a figure later of a model run by running the script ``make_figure.py``. You can either specify a directory that contains output datafiles or an output datafile directly using:

````sh
python make_figure output_directory_or_file
````

Output files can be recognized by the file type .pck. If you specify a directory the script will list all .pck files in this directory and ask you which one to use for making a figure.




## Reference

Please cite the following paper if you publish work that uses PyBasin:

Luijendijk, E., R.T. Van Balen, M. Ter Voorde, P.A.M. Andriessen.
Reconstructing the Late Cretaceous inversion of the Roer Valley Graben
(southern Netherlands) using a new model that integrates burial and
provenance history with fission track thermochronology.
Journal of Geophysical Research 116 (B6)

The paper is freely available and can be found here: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JB008071

A bibtex file with the citation is included in the PyBasin directory 



## Authors
* **Elco Luijendijk**, <elco.luijendijk@geo.uni-goettingen.de>

## License
This project is licensed under the GNU general public license (GPL v3). See the [LICENSE.txt](LICENSE.txt) file for details.
