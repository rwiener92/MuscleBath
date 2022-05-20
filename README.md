# MuscleBath

This repository contains

1. v2p6: a folder with the MATLAB package for reading LabView TDMS files.

2. A version of the MTC_Analysis_Function
    - This was designed for MTC and FMD tissues measured on the Costa Lab WPI Muscle Bath. It reads in the raw LabView TDMS file (v2p6) which should contain two columns (motor position and force). It will output stress-strain values at discrete motor position intervals. *Look inside the function for more details*

3. A version of sigep_fitting for "modulus calculations"
    - This will interpolate the stress-strain output from the MTC_Analysis_Function and fit data for extracting different metrics.

4. Prototyping folder: prototyping code, scripts, and ideas.

