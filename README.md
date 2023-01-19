## SLSN-bol


This repository contains scripts and results of work devoted to processing of photometric observations of superluminous supernovae presented in the Open Supernova Catalog (OSC, https://sne.space/). The snad https://github.com/snad-space/snad package was used for convenient work with data from OSC. Multicolor light curves was approximated by vector Gaussian processes implemented in https://gp.snad.space/. Then bolometric light curves in the black body approximation were calculated using obtained continuous light curves and the least squares method.

## Sample

`SLSN.csv` - table downloaded from the main page of the OSC. It includes objects for which there are spectral data. Initial SLSNe sample is made in `sampling.py`, result is written to `sample_SLSN.csv`. Secondary (according to the results of approximation of multicolor light curves from the primary sample) selection of SLSNe takes place in `second_cut.py`. Final list of selected objects: `second_cut.csv`.

## Approximation

File `approx.py` performs an approximation of all objects included in our sample of superluminous supernovae objects. The results of this script ( `.csv` tables with time and corresponding magnitudes, and figures of multicolor light curves) are in `approx_magn/` directory.

## Bolometric light curves

File `bb.py` defines all formulas needed to build a blackbody model. Black body temperature and radius are determined by the least squares method from the approximated multicolor light curves. File `bolometric.py` initializes methods described in `bb.py` for objects from the sample and writes the result in `bol_output/`. `bol_output/data/` contains tables with following columns: mjd - modified julian date, T - temperature, R - radius (meters), fun - squared errors sum, success - did least squares method converge. `bol_output/figures/` contains figures of temperature, radius and bolometric light curve over time.

## Methods comparison

`PTF12dam+superbol` directory is devoted to comparison our bolometric light curve for superluminous supernova PTF12dam, pseudobolometric and obtained by superbol package (https://github.com/mnicholl/superbol). `PTF12dam+superbol/bb_compare_superbol.py` draws figures of luminosity depending on time obtained by different models. `PTF12dam+superbol/PTF12dam figures/` contains resulting figures.

In order to demonstrate that the approximation of multicolor light curves by common GPs, which are applied to each filter independently, is not applicable in our problem, we build a figure (`gp1d_compare.py`) for one of the SLSNe, which shows results of approximations obtained by one-dimensional GPs and vector GPs.
