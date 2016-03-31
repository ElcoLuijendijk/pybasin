========================================================
pybasin: basin evolution, heat flow and thermochronology
========================================================

:Authors: Elco Luijendijk, Georg-August-Universitaet Goettigen
:version: 0.2, march 2016
 

Downloading pybasin
-------------------

* click the download link on the right for a zip file of the source code
* or clone the repository

Dependencies
------------

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


Installing & running pybasin
----------------------------
* make sure you have installed the necessary python modules. the easiest is to use a python distribution that already includes all the necessary modules, like Anaconda (), pythonxy (https://code.google.com/p/pythonxy/) or canopy (https://www.enthought.com/products/canopy/).
* unzip the pybasin source code
* navigate to the pybasin folder and execute pybasin by executing the following command from the command line:

	>>> python pybasin.py
	

* the model will now run with the default input dataset from the Swiss Molasse Basin. Check the model output in the folder ``model_output/MB``
* the first time you run pybasin, you will have to compile the fission track annealing module. The fission track annealing module was written in fortran instead of python to reduce computation time. Compile the fortran file by navigating to the subfolder `pybasin/lib``, opening a terminal and running the following command:

    f2py -c calculate_reduced_AFT_lengths.f90 -m calculate_reduced_AFT_lengths

* you may need to install a fortran compiler first. For linux operating systems this should be included in your distribution, for Mac OSX you can find installers for Gfortran here: https://gcc.gnu.org/wiki/GFortranBinariesMacOS


Model input & output
--------------------

input:
~~~~~~

* adjust the default location of the model input files by changing the file ``default_input_folder.txt``. The folder name specified in this file should be a subfolder of the folder ``pybasin/input_data``
* adjust the text files in the folder
* adjust the model parameters in the file ``pybasin_params.py``. Use a text editor or preferably a dedicated Python editor like Spyder, PyCharm (https://www.jetbrains.com/pycharm/) or Enthought Canopy (https://www.enthought.com/products/canopy/) to adjust the model parameters.
* if you want to run multiple parameter sets and/or wells in one go, adjust the ``model_scenarios.py`` file and change the parameters you wish to change from ``exhumation_magnitudes = [None]`` to ``exhumation_magnitudes = [500.0, 1000.0, 2000.0]`` to run the model with exhumation magnitudes of 500, 1000 and 2000 m, respectively. For running multiple wells change the ``wells`` list in the same file.

output:
~~~~~~~

* The model generates a series of .csv files that shows the stratigraphy, burial depths and formation thicknesses over time. These files are stored in ``model_output/model_folder/burial_history_csv_files/``, where ``model_folder`` is the name of the input and output folder that you have specified in the file ``default_input_folder.txt``
* If ``make_model_data_fig = True`` in the file pybasin_params.py, the model script will generates a single figure for each model run showing burial depth, temperature and modeled vitrinite relfectacne and apatite fission track age data. THis figure is saved in the folder ``model_output/model_folder/model_data_fig/``
* If ``save_model_run_data = True`` in the file pybasin_params.py, the model script will store all model simulation data for each model run in a datafile that uses the python pickle module. THis file can be read later using a additional python script (still in the works...). The lcoation where these files should be saved is specified in the ``pybasin_params.py`` in line 19: ``datafile_output_dir = '../../heavy_data/pybasin_MB'``


Reference
---------

Please cite the following paper if you publish work that uses pybasin:

Luijendijk, E., R.T. Van Balen, M. Ter Voorde, P.A.M. Andriessen.
Reconstructing the Late Cretaceous inversion of the Roer Valley Graben
(southern Netherlands) using a new model that integrates burial and
provenance history with fission track thermochronology.
Journal of Geophysical Research 116 (B6)

a bibtex file with the citation is included in the pybasin folder 



