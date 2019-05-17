# PyBasin: numerical model for basin evolution, heat flow and thermochronology


## Downloading PyBasin

* click the download link on the right for a zip file of the source code
* or clone the repository


## Installing & running PyBasin

* Make sure you have installed the necessary python modules (see list below). The easiest is to use a python distribution that already includes all the necessary modules, like Anaconda (https://www.continuum.io/downloads), pythonxy (https://code.google.com/p/pythonxy/) or canopy (https://www.enthought.com/products/canopy/).
* Navigate to the PyBasin directory and execute PyBasin by executing the following command from the command line:

````sh
python pybasin.py
`````
	

* The model will now run with the default input dataset from the Roer Valley Graben. Check the model output in the directory ``model_output/example_dataset_1``
* Optional: PyBasin includes a version of the fission track module that is written in Fortran instead of Python. The Fortran version runs much faster than the relatively slow Python version. However, the Fortran module needs to be compiled first to be able to use it. Compile the Fortran file by navigating to the subdirectory ``lib``, opening a terminal and running the following command:

````sh
f2py -c calculate_reduced_AFT_lengths.f90 -m calculate_reduced_AFT_lengths
`````
	

* To use the Fortran version of the fission track annealing library you may need to install a fortran compiler first. For linux operating systems this should be included in your distribution, for Mac OSX you can find installers for Gfortran here: https://gcc.gnu.org/wiki/GFortranBinariesMacOS. I have not tested compiling Fortran code on windows systems. Note that you do not need to compile anything to model vitrinite reflectance and apatite (U-Th)/He data.


## Dependencies

PyBasin requires the following Python packages:

- Numpy: http://www.numpy.org/
- Pandas: https://pandas.pydata.org/
- Scipy: https://www.scipy.org/
- Matplotlib: https://matplotlib.org/


## Running single or multiple models:

The command for running pybasin models is:

````sh
python pybasin.py input_directory -w well1,well2,well3
````

``input_directory`` is the directory that contains all input files. For one of the example datasets this should be ``example_dataset_1`` or ``example_dataset_2``. With the optional command line option -w you can specify which wells to run. This can either be a single well or a list of wells separated by a comma. If you do not specify an input directory, PyBasin will use the default input directory defined in the file ``default_folder.txt``. If you do not specify which well to run at the command line PyBasin will look for a list of wells in the file ``pybasin_params.py`` in your input directory.

	

## Model input 

PyBasin reads input parameters from a single parameter file, ``pybasin_params.py``, along with stratigraphy, lithology and thermal properties of the subsurface from a series of .csv files. 


### Example datasets

PyBasin contains two example datasets, one dataset from the Roer Valley Graben in the southern Netherlands that can be found in the directory ``example_dataset_1`` (Luijendijk et al., 2011, JGR 116(B6)) and a second dataset from the Molasse Basin in Switzerland that can be found in the directory ``example_dataset_2`` (von Hagke et al. 2012, Tectonics 31(5)). You can run these example models using:

````sh
python pybasin.py example_dataset_1
````
and
````sh
python pybasin.py example_dataset_2
````

The model runs will automatically generate figures of the modelled burial and thermal history and vitrinite reflectance, apatite fission track and/or apatite (U-Th)/He data, which can be found in the directory ``model_output/example_dataset_x``.


### Model parameter file

The main model parameters are located in a file called ``pybasin_params.py``. See the directory ``example_dataset_1`` or ``example_dataset_2`` for two examples. This is a Python file, use a text editor or dedicated Python editor like Spyder or PyCharm to adjust the model parameters. The file contains two classes. The class ``ModelParameters`` contains all parameters needed for a single model run. The class ``ParameterRanges`` is optional and can be used to set up multiple model runs. See the section on running multiple models below for more information.


### Input data files

* ``well_stratigraphy.csv``. Contains the depth and name of stratigraphic units for one or more wells or surface outcrops.
* ``stratigraphy_info.csv``. Contains the age and lithology for each stratigraphic unit. Optionally one can include a series of provenance age histories here, with one or more values for provenance_age start and provenance_age_end, that denote the age (Ma) at which the sediments that make up each unit were at 120 ˚C (=roughly 4 km depth at normal geothermal gradients) and the surface, respectively.
* ``lithology_properties.csv``. Contains data on the density, porosity at the surface, compressibility, thermal conductivity, heat capacity and heat production values for each lithological unit. This also should include one row with the properties of pore water.

**optional input files**
* ``vitrinite_reflectance.csv``. Contains columns for depth, vitrinite reflectance and the 1 sigma uncertainty of vitrinite reflectance for samples from one or more wells or surface outcrops.
* Apatite fission track (AFT) data:
	* ``aft_samples.csv``: Contains data on sample names, depth, age, length for one or more AFT samples. Each sample takes one row.
	* ``aft_data.csv``: Contains data on the sample name, single grain AFT age, plus/minus one standard error for each single grain age, and Dpar values. This file contains one row for each single grain age.
* Apatite (U-Th)/He (AHe) data:
	* ``ahe_samples.csv``: Sample names and depths
	* ``ahe_data.csv``: sample name, uncorrected age (Ma) and ±1 standard error, corrected age and ±1 standard error, grain radius (μm), U, Th and Sm content (ppm). This file should contain one row per single grain age.


## Running multiple models

Optionally you can start automated runs to test a range of parameter combinations. This is useful for automated sensitivity or uncertainty analysis. 

The model input file contains a class called ``ParameterRanges``. Any parameter value that is included in this class will be used as input for a single model run. All results will be stored and written to a comma separated (.csv) file names ``model_output/model_params_and_results_x_runs.csv``. 

You can include any model parameter in the automated runs. Simply copy a parameter from the ``ModelParameters`` class into the ``ParameterRanges`` class, add _s to the parameter name and add square brackets around the parameter value. For instance to test multiple values of exhumed thickness, add `exhumed_thicknesses_s = [[500.0], [1000.0]]` to test the effect of exhuming 500 and 1000 m on VR, AFT or AHe thermochronometers. See the file ``pybasin_params.py`` in ``model_input/example_dataset_1`` for an example.

There are two options for running multiple model runs. The default is a sensitivity run. In each model run a single parameter will be changed, while all other parameters are kept constant at the default value specified in the ``ModelParameters`` class. Alternatively you can test all parameter combinations by changing `parameter_combinations = False` to `parameter_combinations = True`. Note that this will generate a lot of model runs, testing ten parameter values for two parameters each will generate 10*10 = 100 model runs, for three parameters this increase to a 1000 model runs, etc...



### Using multiple processors

PyBasin includes an option to distribute model runs over multiple processors. To enable parallel processing, change parameter ``parallel_model_runs = False`` to ``parallel_model_runs = True`` in the ``ParameterRanges`` class in the ``pybasin_params.py`` file. You can specify how many processes you want to use simultaneously using the parameter ``max_number_of_processes``. PyBasin will generate one process for each individual model run. Using multiple processes can signficantly speed up multiple model runs. Note that by default if multiple processes is enabled the extensive screen output that is generated during model runs is redirected to a log file that is saved in a subdirectory of the model output directory. 


## Model output

Each model run will generate a number of output files that are saved to the directory ``model_output/model_directory/``, where ``model_directory`` is a directory that can be specified in the ``pybasin_params.py`` file.

### Figures

If ``make_model_data_fig = True`` in the file PyBasin_params.py, the model script will generate a single figure for each model run showing burial depth, temperature and modeled vitrinite reflectance, apatite fission track age and apatite (U-Th)/He data. 


### Output data files:

The results of each model run are stored in a file named ``model_results_date_well_name_ms0-x_final.csv``. The file contains a copy of all input parameters for each model run, along with model statistics on the goodness of fit (GOF) of the modelled and measured temperature, vitrinite reflectance, AFT or AHe data and values for the modelled temperatures. More detailed model output can be found in a series of .csv files that record the stratigraphy, burial depths, formation thicknesses and modelled temperatures over time. These files are stored in ``model_output/model_directory/thermal_history_csv_files``.


### Binary data file

If ``save_model_run_data = True`` in the file PyBasin_params.py, the model script will store all model simulation data for each model run in a datafile that uses the python pickle module. This file can be read later using a additional python script (still in the works...). The location where these files should be saved is specified in the ``PyBasin_params.py`` in line 19: ``datafile_output_dir = '../../heavy_data/PyBasin_MB'``. You can also make a figure later of a model run by running the script ``make_figure.py``. You can either specify a directory that contains output datafiles or an output datafile directly using:

````sh
python make_figure output_directory_or_file
````

Output files can be recognised by the file type .pck. If you specify a directory the script will list all .pck files in this directory and ask you which one to use for making a figure.



## Reference

Please cite the following paper if you publish work that uses PyBasin:

Luijendijk, E., R.T. Van Balen, M. Ter Voorde, P.A.M. Andriessen. 2011. Reconstructing the Late Cretaceous inversion of the Roer Valley Graben (southern Netherlands) using a new model that integrates burial and provenance history with fission track thermochronology. Journal of Geophysical Research 116 (B6). DOI:10.1029/2010JB008071 

The paper is freely available and can be found here: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JB008071

A bibtex file with the citation is included in the PyBasin directory 



## Authors
* **Elco Luijendijk**, <elco.luijendijk-at-geo.uni-goettingen.de>

## License
This project is licensed under the GNU lesser general public license (LGPL v3). See the [LICENSE.txt](LICENSE.txt) file for details.
