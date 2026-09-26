# PyBasin repository review

Review of the PyBasin repository for improvements and user friendliness, 26 September 2026. Items are grouped by priority. Each item names the location in the code, the effect for the user and the suggested fix.

## High priority

### 1. Inconsistent path handling

The input folder is joined to the script directory (`pybasin.py:1057`), so a relative path pointing outside the repository silently resolves inside the repository and fails; only absolute paths work for input data kept elsewhere. At the same time `output_dir` is used exactly as given (`pybasin.py:1109`), so it is relative to the current working directory, and `runtime_<well>.txt` is written straight into the current working directory (`pybasin.py:1290`), which is why files such as `runtime_E40.txt` end up in the repository root.

Fix: resolve the input path relative to the current working directory first and fall back to the script directory, resolve `output_dir` relative to the script directory when it is a relative path, and write the runtime file into `output_dir`.

### 2. A misspelled parameter in pybasin_params.py fails silently

There are around 30 `getattr(pybasin_params, ...)` calls with default values (17 in `lib/pybasin_lib.py`, 6 each in `pybasin.py` and `lib/model_input_io.py`). Writing `simulate_fluid_flwo = True` or `lithosphere_base_temp = 1330` therefore runs the default model without any warning, and the user gets a plausible but wrong result.

Fix: a startup check that compares the attribute names in the user's `ModelParameters` class against the set of names that PyBasin actually reads, and that logs a warning for every name it does not recognise.

### 3. Well names are validated too late

`select_well_strat` is called inside the per-well loop (`pybasin.py:1226`), so a typo in the second of two wells only surfaces after the first well has run to completion, which can take many minutes. The error message itself (`lib/model_input_io.py:433`) runs two sentences together without a separator and does not list the well names that are available. In addition, `-w "NDW-01, NDW-02"` is split on commas without stripping whitespace (`pybasin.py:1153`), so the second name becomes `" NDW-02"` and is not found.

Fix: validate every requested well directly after `read_model_input_data`, list the available wells in the error message, and strip whitespace from each name given on the command line.

### 4. The manual has fallen behind the code

`manual/PyBasin_manual.md:77` still states that PyBasin contains two example datasets, and the `simulate_lithosphere` group of parameters used by example dataset 5 (`lithosphere_base_temperature`, the heat production parameters and the lithosphere grid parameters) is not documented. The front matter date is still 22 May 2019. The parameter reference sections for `simulate_fluid_flow` and `compaction_method` are current, so it is mainly the getting started and example dataset sections plus the new lithosphere parameters that need an update.

Fix: update those sections and regenerate the pdf.

## Medium priority

### 5. Stale counts and pointers

`readme.md:53` says that there are four example datasets and then lists five. The error message in `pybasin.py:1080` refers the user to "example_dataset_1 to input_data/example_dataset_4". The header of `pybasin.py` still gives version 0.1 and an old Goettingen email address, while `CITATION.cff` and the readme give the current one.

### 6. There is no pip installation route

There is no `pyproject.toml` or `setup.py`, so the code can only be run from the source directory and cannot be imported as a package from a script or notebook kept elsewhere. `environment.yml` also pins exact versions, which is good for reproducibility but will eventually fail to solve.

Fix: a minimal `pyproject.toml` with a console script entry point, plus a `requirements.txt` with lower version bounds for users who prefer pip over conda.

### 7. Two bare except clauses hide real problems

`lib/AFTannealingLib.py:43` catches every exception, so a Fortran module that was compiled against a different NumPy version is reported as "not found" and the user silently falls back to the slow Python implementation. `pybasin.py:1787` swallows any failure while splitting output columns.

Fix: catch `ImportError` specifically and report the underlying exception at info level instead of debug level.

### 8. The Fortran fallback message is printed twice at startup

The same message is emitted from `lib/AFTannealingLib.py:49` and from `lib/helium_diffusion_models.py:17`. `AFTannealingLib` already contains a `_warn_once` helper that the import time code path does not use.

Fix: emit the message from a single place.

## Low priority

### 9. Command line gaps

`-w` has no `--wells` long form (`pybasin.py:1040`), there is no option to override the output directory, and there is no way to list the bundled example datasets. An output directory option in particular would let a user run variations without editing the parameter file.

### 10. default_input_folder.txt can be replaced by an argparse default

The file, the file read and the unused `scenario_name` variable (`pybasin.py:1057` to `pybasin.py:1064`) can all be replaced by a default value in the argument parser.

### 11. Clutter in the repository root

`10.1029%2F2010JB008071.bib` in the root is byte identical to `references/10.1029%2F2010JB008071.bib`, and the empty `__init__.py` in the root no longer serves a purpose now that all imports go through the `lib` package.

### 12. The parameter files are long and largely duplicated

The five `pybasin_params.py` files are 250 to 360 lines each and mostly identical, so a new option has to be copied into all of them and a user cannot easily tell which lines matter for their own case. A single module that lists every parameter with its default value would give one source of truth, and pairs naturally with the unknown parameter check in item 2. Note that the comments in the existing parameter files are useful documentation for users who edit them, so they are worth keeping in place rather than replacing them with a bare set of overrides.
