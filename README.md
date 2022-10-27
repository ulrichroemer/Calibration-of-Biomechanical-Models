# Calibration-of-Biomechanical-Models

The repository contains MATLAB code related to the article:

Römer, Ulrich, Jintian Liu, and Markus Böl. "Surrogate‐Based Bayesian Calibration of Biomechanical Models with Isotropic Material Behavior." International Journal for Numerical Methods in Biomedical Engineering (2022): e3575.

You will need: 

- MATLAB (tests run with version R2020a) 
- UQLab (tests run with version 1.4.0)
- An AIES sampler (UQLab has a build in version)
  the current version includes a routine provided here (https://github.com/grinsted/gwmcmc)

There are three examples: 
	
- anlytical_uq.m can be used to produce the analytical example of Section 3.1
- gel_"xyz".m relate to the example of Section 3.2; gel_uq.m is the main file
- oocyte_"xyz".m relate to the example of Section 3.3; oocyte_uq.m is the main file

Further information: 

- uqlab_"xyz".m are functions calling UQLab routines (Bayesian update, surrogate modeling,...)

- The gel example uses FEM input-output data (folder data) simulated beforehand 
	- The original data is provided in exp_WPG.xlsx
	- This data is preprocessed and loaded in the main file gel_uq.m

- The oocyte example uses FEM input-output data (folder data) simulated beforehand 
	- The original data is provided in daten_fig5_indent_raw
	- This data is preprocessed and loaded in the main file oocyte_uq.m





