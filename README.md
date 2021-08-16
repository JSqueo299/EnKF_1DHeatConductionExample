# EnKF-1D Transient Heat Conduction Example
Esemble Kalman filter data assimilation example involving a 1D transient heat conduction in an aluminum rod.

------------------
### Contributor
Name: Joseph N. Squeo

Contact: joseph.squeo@uconn.edu

Institution: Univeristy of Connectiut (UConn)

Country: USA

------------------
### References
**[1] The 1D heat conduction example and explanation**
S. Gillijns, O. B. Mendoza, J. Chandrasekar, B. L. De Moor, D. S. Bernstein, and A. Ridley,
What Is the Ensemble Kalman Filter and How Well Does it Work?, (2006).

**[2] Deeper analysis and explanation of the ensemble Kalman filter, used to develop future codes with CFD simulations**
J. F. Labahn, H. Wu, B. Coriton, J. H. Frank, and M. Ihme, Data Assimilation Using High-
Speed Measurements and LES to Examine Local Extinction Events in Turbulent Flames,
Proceedings of the Combustion Institute (2017).

**[3] My Thesis - Explanation of the 1D heat conduction example and further examples**
Squeo, J.N. [2021], Data-based soot modeling in buoyancy-driven diffusion flames, Dissertation, The University of Connecticut.

Available here: [Squeo Thesis](https://collections.ctdigitalarchive.org/islandora/object/20002:UConnTheses)
Type "Squeo" into the search bar.

**[4] For a more complete understanding of Kalman filter and the mathematics/Bayesian statistics behind the algorithms**
Geir Evensen. Data assimilation: The ensemble Kalmanf filter. Springer.


------------------
### Explanation of the Code
Data assimilation combines a model, observations and the corresponding uncertainties of a dynamical system to provide anu pdate to the system’s state.
These sequential time-stepping algorithms update the model state to reflect newly received state observations, or measurements, using the previous model forecast. Many different data assimilation techniques exist, serving as a useful tool for state and parameterestimation.

The ensemble Kalman code and example implemented is based on the work of Gillijns et al. [1]. 

Please see section 4.3 from reference [3] for the explanation of the code and example provided here.


------------------
### Navigation of GitHub Repo
This repo contains a number of MATLAB files required to run the 1D heat conduction example. However, this code requires open source OpenFOAM CFD software to run on a Linux or MacOS system.

The main code is titled **OpenFOAM_1DheatConductionExample.m**, where various functions are called to implement the EnKF, navigate repositories and set up the OpenFOAM case. The EnKF code itself is contained in the MATLAB file titled **ensemblekfilter.m**, although this is only used when a finite difference method in MATLAB is used as the model to evolve the state forward in time. If OpenFOAM is used to evolve the state forward in time, the EnKF is contained in the file **OF_ensemblekfilter2.m**. 

**ensemblekfilter_plots.m** is a file to produce plots of the results.

**finiteDiff_1DheatConductionExample.m** is an analytical solution to the 1D heat conduction in an aluminum rod toy problem.

**tempDifference_finiteDiff_OpenFoam.m** is a post-processing script to verify the analytical finite difference model solution of temperature calculated in MATLAB against OpenFOAM's finite volume method. This process is detailed in section 4.3.1 of [3].

Please note, this is a relatively outdated and poorly written version of the EnKF code used simply for testing and understanding the EnKF algorithm from a fundamental level with a simply 1D toy problem. Improved versions of the code were developed for coupling with flame simulations in OpenFOAM on a parallel HPC cluster. However, these codes are much more difficult to understand with more functions and complications due to MATLAB's coupling with OpenFOAM CFD software on an HPC cluster. More information can be found in [3] and [4].

