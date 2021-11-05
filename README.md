
# Earth System Model Aerosol-Cloud Diagnostics Package

This Earth System Model (ESM) aerosol-cloud diagnostics package (ESMAC Diags) is currently used to evaluate aerosols, clouds and aerosol-cloud interactions simulated by the Department of Energy’s (DOE) Energy Exascale Earth System Model (E3SM). The first version (v1.0) focuses on comparing simulated aerosol properties with in-situ aircraft, ship and surface measurements. Various types of diagnostics and metrics are performed for aerosol number, size distribution, chemical composition, and CCN concentration to assess how well E3SM represents observed aerosol properties across spatial scales. Metrics for various meteorological and aerosol precursor quantities from the same field campaigns are also included. Version 2 is under development focusing on aerosol-cloud interactions.

More information can be found in README_ESMAC_Diags_v1.0.pdf

# Package dependencies
This code is dependent on the following python packages:

os
sys
glob
time
numpy
scipy
matplotlib
netCDF4


# Test run
To verify the package, enter scripts/ directory and run scripts_testcase.csh. Then check the directory in testcase/figures/. There should be three figures generated:

flighttrack_ACEENA_20170630a.png

flightheight_ACEENA_20170630a.png

AerosolComposition_ACEENA_20170630a.png

Directory testcase/figures_verify/ contains what the three figures should look like. If the three figures in testcase/figures/ are consistent with figures_verify/, the testcase is successfully run.


