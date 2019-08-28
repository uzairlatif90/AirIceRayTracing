# AnalyticRayTracing from Air to Ice
C++ code which uses analytic ray tracing for tracing rays from any point in the air to any point in the ice

- MultiRayAirIceRefraction.C : This script takes in as arguments the antenna depth and the height of the ice layer. Then it loops over values of Tx Height and the ray launch angle and prints a file called "TableValues.txt" which contains the following columns:

  - the entry number, the Tx height, ice layer height, Tx height above the ice layer height, total horizontal distance traveled on the surface above ice, total horizontal distance traveled in ice, the initial ray launch angle at Tx, incident angle on ice and received angle in ice at the antenna inside this file

  - The step size for Tx height loop is set at 20 m and it starts at the maximum available height from the refractive index data to a step above the ice surface.

  - The step size for the ray launch angle loop goes from 91 deg to 179 deg in steps of 1 deg. Here 0 deg is vertically upwards

- SingleRayAirIceRefraction.C : This script takes in as arguments the antenna depth, ice layer height, initial launch angle of the ray and height of the Tx. It traces the ray for this particular configuration and then prints out the total horizontal distance that was travelled by the ray above the ice and inside the ice. It also makes a text file called "RayPathinAirnIce.txt" which contains the x (distance), y (height) values (in m) of the ray path as it traverses through the atmosphere and ice.

- SingleRayAirIceRefraction_wROOTGr.C : It is the same as the above script except that the only difference is that it makes plots of the ray path in ROOT and it can be run with just ROOT by doing: root -l SingleRayAirIceRefraction_wROOTGr.C

- Atmosphere.dat: This file has been generated with GDAS tool that comes with CORSIKA. The command that was used to generate this file was: ./gdastool -t 1533600000 -o Atmosphere.dat -c -89.9588 -109.794 -m -5 -v -g

  - These are the coordinates for the ARA2 station at the South Pole.

## Prerequisites
You will need to have a functioning installation of [GSL](https://www.gnu.org/software/gsl/) ([2.4](https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz) is verified to work).
- You will need to set the enviromnent variable `GSLDIR` to your local installation of GSL.
- You will also need to have `GSLDIR` in your `LD_LIBRARY_PATH`.
- For Mac users: you can locate your GSL installation via `gsl-config --prefix`, i.e. `export GSLDIR=$(gsl-config --prefix)`
- If you have Ubuntu or Linux you can skip all of the above and just get it from the repository by doing: "sudo apt install libgsl-dev"

## Install instructions

### SingleRayAirIceRefraction.C as standalone package
To run you just have to do:
- Make it: `make SingleRayAirIceRefraction`
- Run it: `./SingleRayAirIceRefraction 200 170 20000 3000`
- In this case the example arguments are: Antenna Depth is set at 200 m, The Ray Launch Angle is set at 170 deg, Tx Height is set at 20000 m, Ice Layer Height is set as 3000 m
- The main is at the bottom of the code, which you can modify to your liking.

### MultiRayAirIceRefraction.C as standalone package
To run you just have to do:
- Make it: `make MultiRayAirIceRefraction`
- Run it: `./MultiRayAirIceRefraction 200 3000`
- In this case the example arguments are: Antenna Depth is set at 200 m, Ice Layer Height is set as 3000 m
- The main is at the bottom of the code, which you can modify to your liking.
