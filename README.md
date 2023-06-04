# cryo-evaporation-foam
OpenFOAM and Python code to simulate the non-isobaric evaporation of cryogenic liquids in storage tanks. This code is part of the supplementary material for the research article "CFD modelling of the non-isobaric evaporation of cryogenic liquids in storage tanks" written by Felipe Huerta and Velisa Vesovic.
   
Cite the code:  [![DOI](https://zenodo.org/badge/649334155.svg)](https://zenodo.org/badge/latestdoi/649334155)

## Features
* Simulates the non-isobaric evaporation of pure cryogenic liquids in storage tanks.
* The model is valid for vertically orientated cylinders.
* The heat ingress through the bottom, walls and top are input parameters.
* An evaporative fraction that partitions the vapour wall heat ingress in wall boiling and vapour heating.
* The code is written in Python and uses the OpenFOAM finite volume library.
The code is accompanied by a set of example cases that demonstrate its use.

## Prerequisites

* Python >= 3.8.
* OpenFOAM >= v2006. For installation, follow the instructions [in the official OpenFOAM website](https://develop.openfoam.com/Development/openfoam/-/blob/master/doc/Build.md). Installation from source is recommended.
* Paraview >= 5.11. For installation, follow the instrutions [in the official Paraview website](https://www.paraview.org/)
* Gmsh

## Installation
1. Install Python, OpenFOAM and paraview.
2. Clone the cryo-evaporation-foam repository.
3. Run a shell in the source directory and navigate to the solver folder

   `cd buoyantBoussinesqPvapFoam`

4. Compile the buoyantBoussinesqPvapFoam solver
```shell
wclean
wmake
```

## Usage

To illustrate the use of the solver an example to run the case with the shortest simulation time is presented.

1. From the root directory, navigate to Seo/LF_30 directory

   `cd Seo/LF_30`

2. Run the shell script to generate the mesh from the cyl_structured.geo file

   `./update_mesh`

3. Run the shell script to run the simulation. Add "&" to separate the process from the shell.

   `./buoy_run &`

4. Create an empty file to monitor the simulation with paraview

   `touch tank.foam`

5. Open ParaView to visualise the results as they are saved every 60 seconds of simulation time

   `paraview tank.foam`

6. You may need to rotate the view in ParaView  to ensure that the z-axis  is pointing upwards and the x-axis increases to the right. 

7. The following figures show the figures displaying the velocity magnitude (m/s) and temperature (K) profiles after 60 s of evaporation.

![U_Seo_LF30_120s](https://github.com/felipehuerta17/cryo-evaporation-foam/assets/33637198/daf67cf7-2639-4ed8-9da7-54c9783c15ab)

![TL_Seo_LF30_120s](https://github.com/felipehuerta17/cryo-evaporation-foam/assets/33637198/c65ecda9-0c01-4f23-841b-4c5bbd6bcce0)

8. The pressure build-up data is generated the directory postProcessing/p_average/0/volFieldValue.dat. The data has two columns and is separated by tabs. You can visualise the data as the simulation progresses from the shell with your favourite shell text editor. For instance, with Vim:

   `vim postProcessing/p_average/0/volFieldValue.dat`

9. The vapour height and vapour temperature are saved for each written time-step directory in the file vap_mesh.csv and vap_temp.csv. For example, to examine the vapour length and temperature values, input the following commands:

   `vim 120/vap_mesh.csv`
   
   `vim 120/vap_temp.csv`


## References for the experimental scenarios:

1. M. Seo and S. Jeong, "Analysis of self-pressurization phenomenon of cryogenic fluid storage tank with thermal diffusion model," Cryogenics, vol. 50, no. 9, pp. 549-555, Sep 2010, doi: https://doi.org/10.1016/j.cryogenics.2010.02.021. 

2. M. Kang, J. Kim, H. You, and D. Chang, "Experimental investigation of thermal stratification in cryogenic tanks," Exp. Therm. Fluid Sci., vol. 96, pp. 371-382, 2018/09/01/ 2018, doi: https://doi.org/10.1016/j.expthermflusci.2017.12.017. 
