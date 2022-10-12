## SLSN-bol


This repository contains scripts and results of work devoted to processing of photometric observations of superluminous supernovae presented in the Open Supernova Catalog. Multicolor light curves was approximated by vector Gaussian processes implemented in https://gp.snad.space/. Then bolometric light curves in the black body approximation were calculated using obtained continuous light curves and the least squares method.

## Approximation

File approx.py performs an approximation of all objects included in our sample of superluminous supernovae objects. The results of this script (.CSV tables with time and corresponding magnitudes, and figures of multicolor light curves) are in approx_magn/ directory.
