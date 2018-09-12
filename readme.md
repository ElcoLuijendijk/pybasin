#PyBasin: basin evolution, heat flow and thermochronology

:Authors: Elco Luijendijk, Georg-August-Universitaet Goettigen

## 

## Downloading pybasin

* click the download link on the right for a zip file of the source code
* or clone the repository

## Dependencies

PyBasin requires the following packages to be installed:

- Python2.x_	Contributes the Py to PyBasin		 
- Numpy_	Fast arrays
- Pandas_	Reading model input and output and data processing
- Scipy_	Statistical and model calibration algorithms
- Matplotlib_	External library for making model result figures

.. _Python2.x: http://www.python.org/
.. _Numpy: http://www.scipy.org/NumPy
.. _Pandas: http://pandas.pydata.org
.. _Scipy: http://www.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/


## Installing & running pybasin

* make sure you have installed the necessary python modules. the easiest is to use a python distribution that already includes all the necessary modules, like Anaconda (https://www.continuum.io/downloads), pythonxy (https://code.google.com/p/pythonxy/) or canopy (https://www.enthought.com/products/canopy/).
* Download and unzip the pybasin source code
* or clone the repository to keep your version up to date with the pybasin version on the bitbucket. See "Cloning a Git repository" on this webpage for more info on cloning a repository: https://confluence.atlassian.com/bitbucket/clone-a-repository-223217891.html  
* Navigate to the pybasin folder and execute pybasin by executing the following command from the command line:

	>>> python pybasin.py
	

* The model will now run with the default input dataset from the Swiss Molasse Basin. Check the model output in the folder ``model_output/MB``
* Optionally:
	* If you want to model apatite fission track data, you will have to compile the fission track annealing module first. The fission track annealing module was written in fortran instead of python to reduce computation time. Compile the fortran file by navigating to the subfolder `pybasin/lib``, opening a terminal and running the following command:

    f2py -c calculate_reduced_AFT_lengths.f90 -m calculate_reduced_AFT_lengths

	* you may need to install a fortran compiler first. For linux operating systems this should be included in your distribution, for Mac OSX you can find installers for Gfortran here: https://gcc.gnu.org/wiki/GFortranBinariesMacOS. I have not tested compiling fortran code on windows systems.
	* Note that you do not need to compile anything for modeling vitrinite reflectance and apatite (U-Th)/He data

## Model input 

* adjust the default location of the model input files by changing the file ``default_input_folder.txt``. The folder name specified in this file should be a subfolder of the folder ``pybasin/input_data``
* adjust the text files in the folder to adjust the well stratigraphy, lithology and thermal properties of the rocks that you are trying to model
* adjust the model parameters in the file ``pybasin_params.py``. Use a text editor or preferably a dedicated Python editor like Spyder, PyCharm (https://www.jetbrains.com/pycharm/) or Enthought Canopy (https://www.enthought.com/products/canopy/) to adjust the model parameters.
* if you want to run multiple parameter sets and/or wells in one go, adjust the ``model_scenarios.py`` file and change the parameters you wish to change. For instance, changing from ``exhumation_magnitudes = [None]`` to ``exhumation_magnitudes = [500.0, 1000.0, 2000.0]`` will tell pybasin to execute three a single model runs in one go with exhumation magnitudes of 500, 1000 and 2000 m, respectively. For automatically executing model runs for multiple wells change the ``wells`` list in the same file to include several wells.

## Model output:

* The model generates a series of .csv files that show the stratigraphy, burial depths and formation thicknesses over time. These files are stored in ``model_output/model_folder/burial_history_csv_files/``, where ``model_folder`` is the name of the input and output folder that you have specified in the file ``default_input_folder.txt``
* If ``make_model_data_fig = True`` in the file pybasin_params.py, the model script will generates a single figure for each model run showing burial depth, temperature and modeled vitrinite relfectance, apatite fission track age and apatite (U-Th)/He data. THis figure is saved in the folder ``model_output/model_folder/model_data_fig/``
* If ``save_model_run_data = True`` in the file pybasin_params.py, the model script will store all model simulation data for each model run in a datafile that uses the python pickle module. This file can be read later using a additional python script (still in the works...). The lcoation where these files should be saved is specified in the ``pybasin_params.py`` in line 19: ``datafile_output_dir = '../../heavy_data/pybasin_MB'``


## Reference

Please cite the following paper if you publish work that uses pybasin:

Luijendijk, E., R.T. Van Balen, M. Ter Voorde, P.A.M. Andriessen.
Reconstructing the Late Cretaceous inversion of the Roer Valley Graben
(southern Netherlands) using a new model that integrates burial and
provenance history with fission track thermochronology.
Journal of Geophysical Research 116 (B6)

a bibtex file with the citation is included in the pybasin folder 



## Authors
* **Elco Luijendijk**, <elco.luijendijk@geo.uni-goettingen.de>

## License
This project is licensed under the GNU general public license (GPL v3). See the [LICENSE.txt](LICENSE.txt) file for details.
