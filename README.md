# PyFDAP

PyFDAP is a Python-based analysis framework for Fluorescence Decay After Photoconversion (FDAP) experiments. 

<img style="float: right;" width=300px src="pyfdap/logo/pyfdap_icon.png">

FDAP is a microscopy-based technique for mea-suring protein half-lives. In FDAP experiments, a protein of interest
is tagged with a photoconvertible fluorescent protein and expressed in vivo. The fluorescent
fusion protein is then photoconverted, and the decrease in fluorescence intensity over time is
monitored. The resulting intensity data is fitted with a decay function, and half-lives can be
calculated from the fits.

Both intracellular and extracellular protein half-lives can be determined using FDAP. A static
intracellular signal (e.g. Alexa488-dextran) can be used to create an intracellular mask, such
that only intracellular pixels are considered when calculating intracellular intensity. The mask
can be inverted to calculate extracellular intensities

To learn more about FDAP, please refer to 

[Rogers KW, Bläßle A, Schier AF, Müller P (2015). Measuring protein stability in living zebrafish embryos using Fluorescence Decay After Photoconversion (FDAP). J Vis Exp, 95:52266. ](http://www.jove.com/video/52266/measuring-protein-stability-living-zebrafish-embryos-using)

## Main Features

1. A comprehensive data format for handling, sorting, and anno-tating large FDAP datasets 
2. The capability to separate fluorescence intensities in FDAP datasets into intra- and extracellular compartments based on counter-labeling
3. Model fitting via established fitting algorithms 
4. A user-friendly environment that allows researchers from a non-computational background to easily evaluate FDAP datasets.

## Documentation

PyFDAPs documentation is available online on the [PyFDAP website](http://people.tuebingen.mpg.de/mueller-lab/#content/home.html ) or as a PDF [here](http://people.tuebingen.mpg.de/mueller-lab/downloads/manual.pdf).

## Obtaining PyFDAP

You can download the PyFDAP source code via

	git clone https://github.com/alexblaessle/PyFDAP.git
	
or download PyFDAP as  binary executables from [here](http://people.tuebingen.mpg.de/mueller-lab/#content/downloads.html).

## Required packages

If you run PyFDAP from source, you will need to have the following packages installed

- numpy
- scipy
- matplotlib
- PyQt4
- scikit-image

If you want to make use of PyFDAPs video output, you will need

- mencoder (ffmpeg).

If you run PyFDAP from an executable, all necessary packages will come bundled inside the executable (except mencoder).

## Installing PyFDAP

Detailed explanation on howto install all necessary dependencies can be found in the [installation section of the documentation](http://people.tuebingen.mpg.de/mueller-lab/#content/installation.html).

## Running PyFDAP

PyFDAP can be run in three different ways:

1. As a Python package via

```python
import pyfdap
pyfdap.pyfdap_app.main()
```

2. By executing the downloaded files directly

``bash
cd path/to/pyfdap_root_folder
python pyfdap/pyfdap_app.py
``

3. Or by downloading an executable binary and launching PyFDAP by double-clicking.


## Citing PyFDAP

If you use PyFDAP for your research, please cite the original publication:

[Bläßle A, Müller P. PyFDAP: automated analysis of fluorescence decay after photoconversion (FDAP) experiments. *Bioinformatics*. 2015 Mar 15;31(6):972-4. doi: 10.1093/bioinformatics/btu735. Epub 2014 Nov 6.](http://bioinformatics.oxfordjournals.org/content/early/2014/12/01/bioinformatics.btu735.abstract)

## Updates Version 1.1

### Major Changes:
- Changed the way average fits are computed:
  -> Data is now pooled in time bins
  -> Average fit is newly computed and is a fit to average data, instead previously an average of the fits.
- PyFDAP automatically imports useful_fcts as uf into terminal
- PyFDAP is now a python packages (has __init__.py and setup.py)
- Added constrained Nelder-Mead algorithm.

### Minor Changes:
- Changed the way some data is getting backuped before saving a molecule or an embryo to increase performance and find a workaround to memory overloads.

### Major Additions:
- Added mutliple statistical tests:
	* Shirapo
	* t-test
	* Welch's t-test
	* Wilcoxon
	* Mann-Whitney-U
- Added option to load a list of timeseries from a csv and created a list of embryos out of it
- Added option to load a list of timeseries from a csv and created a list of bkgds out of it
- Added Corresponding CSV Input dialogs
- Added option for different threshholding algorithms (fixed,otsu,adaptive)
- Added threshhold selector for better threshhold selection
- Added option to fill values outside of embryo with nan to become better histogram and improve threshholding performance
- Added function that checks if objects (molecule,embryo,fits,...) are of newest version and adds latest properties to object if not yet present
- Added copy/paste embryo feature. Now can also copy embryo objects between molecule objects.
- Added copy/paste background feature. Now can also copy background objects between molecule objects.
- Added shortcut to copy and paste quickly any objects.
- Added plot extrapolated fit function. Plots extrapolated fits from (0,tend) and highlights c0, y0 and tau.

#Minor Additions:
- Added plot all ext/int/slice data. Plots all data + background datasets in one single plot.
- Now fit object keeps track of number of iterations and function calls during fitting.
- fit object now has print_results().
- Added I_pre as a option for x0_y0.
- Added auto completion for pre folder dialog.
- Added functions get_max_int_pxs, label_pxs to create range indicator masks .
- Added range indicator to dataset_dialog.
- Added two new export functions: Export selected object to csv and export selected fits to csv.
- Added the option to also ignore frames from background datasets
- Added option to launch PyFDAP GUI and not redirect terminal stout/sterr to PyFDAP terminal
- Added option to launch PyFDAP GUI and ignore warning messages for FutureWarning and DeprecatedWarning

#Bugfixes: 
- Fixed Bug where analyze all datasets resulted in error if dataset had no background dataset.
- Fixed some unnecessary print outs.
- Fixed save dialog to show .pk files and also autocomplete filename.
- Fixed a bug in copy_fit where fit.fit_number was getting assigned in the wrong way.   
