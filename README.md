# AnalyticRayTracing from Air to Ice
C++ code that uses analytic ray tracing for tracing rays from any point in air to any location in the ice

- RayTracingFunctions.cc , RayTracingFunctions.h: These two files contain functions for a namespace "RayTracingFunctions" that I made for all of these scripts. The functions use the raytracing algorithms to trace rays from any point in ice to any point in air.

- MultiRayAirIceRefraction.cc , MultiRayAirIceRefraction.h: These two files contain functions for a namespace "MultiRayAirIceRefraction" that I made for CoREAS. The functions use the raytracing algorithms to trace rays from any point above the ice surface in air to any point below the ice surface. There are two ways to do that:
  1) is to use a function that generates a table first and then interpolates values through raytracing. To make the tables, the step size for Tx height loop is set at 20 m, and it starts at the maximum available height from the refractive index data to a step above the ice surface. And the step size for the ray launch angle loop goes from 92 degrees to 180 degrees in steps of 0.5 deg. Here 0 deg is vertically upwards. The table is stored c++ vectors and contains:
     - the entry number, the Tx height, total horizontal distance traveled, total horizontal distance traveled on the surface above the ice, total horizontal distance traveled in ice, total propagation time, total propagation time in air, total propagation time in ice, the initial ray launch angle at Tx, incident angle on ice, received angle in ice at the antenna inside this file, Reflection Coefficient S, Reflection Coefficient P

  2) is to use raytracing function directly that contains the minimization

- RunMultiRayCode.C: This script takes in as arguments the antenna depth, the height of the ice layer, Air Tx Height and Horizontal Distance btw Air Tx and Antenna Rx in ice. Then you can use 1) or 2) to find a possible solution (ray path) between the given two points.

- SingleRayAirIceRefraction.C : This script takes in as arguments the antenna depth, ice layer height, initial launch angle of the ray, and height of the Tx. It traces the ray for this particular configuration and then prints out the total horizontal distance traveled by the ray above the ice and inside the ice. It also makes a text file called "RayPathinAirnIce.txt" which contains the x (distance), y (height) values (in m) of the ray path as it traverses through the atmosphere and ice.

- SingleRayAirIceRefraction_wROOTGr.C : It is the same as the above script except that the only difference is that it makes plots of the ray path in ROOT.

- Air2IceRayTracing.C : This script does the same thing as SingleRayAirIceRefraction.C, but the only difference is that it calculates the launch angle itself given a point in air and a point in ice. So it will calculate the whole path for you if you tell it the coordinates of those two points.

  - It takes in as arguments the Tx height in air, horizontal distance btw Tx in air and Rx in ice, the height of the ice layer and depth of the antenna or Rx in ice. It traces the ray for this particular configuration and then prints out the total horizontal distance traveled by the ray above the ice and inside the ice. It also makes a text file called "RayPathinAirnIce.txt" which contains the x (distance), y (height) values (in m) of the ray path as it traverses through the atmosphere and ice.

- Air2IceRayTracing_wROOTplot.C : It is the same as the above script except that the only difference is that it makes plots of the ray path in ROOT.

- AirRayTracing.C : This script calculates the ray path between any two possible points in air

  - It takes in as arguments the Tx height in air, Rx height in air, horizontal distance btw Tx and Rx in air and height of the ice layer. It traces the ray for this particular configuration and then prints out the total horizontal distance that was traveled by the ray. It also makes a text file called "RayPathinAir.txt" which contains the x (distance), y (height) values (in m) of the ray path as it traverses through the atmosphere.

- AirRayTracing_wROOTplot.C : It is the same as the above script except that the only difference is that it makes plots of the ray path in ROOT.

- MakeMultiRayPlot.C : This plots a Multi-Ray diagram of rays launching from a given transmitter (Tx) position. It takes in the given Tx depth and angle step size to launch rays and launches rays from 0 to 90 degrees and traces the rays from ice to air. Here the angle is measured w.r.t vertical where 0 deg is straight up. This macro depends on IceRayTracing.cc and RayTracingFunctions.cc.

- Atmosphere.dat: This file has been generated with GDAS tool that comes with CORSIKA. The command that was used to generate this file was: ./gdastool -t 1533600000 -o Atmosphere.dat -c -89.9588 -109.794 -m -5 -v -g

  - These are the coordinates for the ARA2 station at the South Pole.

## Prerequisites
You will need to have a functioning installation of [GSL](https://www.gnu.org/software/gsl/) ([2.4](https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz) is verified to work).
- You will need to set the environment variable `GSLDIR` to your local installation of GSL.
- You will also need to have `GSLDIR` in your `LD_LIBRARY_PATH`.
- For Mac users: you can locate your GSL installation via `gsl-config --prefix`, i.e. `export GSLDIR=$(gsl-config --prefix)`
- If you have Ubuntu or Linux you can skip all of the above and get it from the repository by doing: "sudo apt install libgsl-dev"

## Install instructions

### SingleRayAirIceRefraction.C as standalone package
To run, you have to do:
- Make it: `make SingleRayAirIceRefraction`
- Run it: `./SingleRayAirIceRefraction 200 170 20000 3000`
- In this case, the example arguments are: Antenna Depth is set at 200 m, The Ray Launch Angle is set at 170 deg, Tx Height is set at 20000 m, Ice Layer Height is set as 3000 m
- The main is at the bottom of the code, which you can modify to your liking.

### MultiRayAirIceRefraction namespace functions from .cc and .h files
To use the functions you just have to call the example RunMultiRayCode.C that calls the two main important functions:
- Run it: 'root -l RunMultiRayCode.C'

### SingleRayAirIceRefraction_wROOTGr.C as standalone package
To run you just have to do:
- `root -l 'SingleRayAirIceRefraction_wROOTGr.C('200','170','20000','3000')'`
- In this case the example arguments are: Antenna Depth is set at 200 m, The Ray Launch Angle is set at 170 deg, Tx Height is set at 20000 m, Ice Layer Height is set as 3000 m

### Air2IceRayTracing.C as standalone package
To run, you just have to do:
- Make it: `make Air2IceRayTracing`
- Run it: `./Air2IceRayTracing 5000 1000 3000 200`
- In this case, the example arguments are: Tx Height in air is set at 5000 m, the horizontal distance btw Tx in air and Rx in ice is set at 1000 m, Ice Layer Height is set at 3000 m, Antenna Depth is set at 200 m
- The main is at the bottom of the code, which you can modify to your liking.

### Air2IceRayTracing_wROOTplot.C as standalone package
To run, you have to do:
- `root -l 'Air2IceRayTracing_wROOTplot.C('5000','1000','3000','200')'`
- In this case, the example arguments are: Tx Height in air is set at 5000 m, the horizontal distance btw Tx in air and Rx in ice is set at 1000 m, Ice Layer Height is set at 3000 m, Antenna Depth is set at 200 m

### AirRayTracing.C as standalone package
To run, you have to do:
- Make it: `make AirRayTracing`
- Run it: `./AirRayTracing 5000 3100 1000 3000`
- In this case, the example arguments are: Here 5000 m is Tx Height in air in m, 3100 m is Rx Height in air, 1000 m is the horizontal distance btw Tx in air and Rx in air in m and 3000 m is Ice Layer Height
- The main is at the bottom of the code, which you can modify to your liking.

### AirRayTracing_wROOTplot.C as standalone package
To run, you have to do:
- `root -l 'AirRayTracing_wROOTplot.C('5000','3100','1000','3000')'`
- In this case, the example arguments are: Here 5000 m is Tx Height in air in m, 3100 m is Rx Height in air, 1000 m is the horizontal distance btw Tx in air and Rx in air in m and 3000 m is Ice Layer Height

### MakeMultiRayPlot.C as standalone package
To run, you have to do:
- `root -l 'MakeMultiRayPlot.C('-180','0.25')'`
- Here the two arguments are the Tx depth (i.e., Tx_z=-180 m) and the launch angle step size (i.e., 0.25 deg)
- Please note that MakeMultiRayPlot depends on IceRayTracing.cc and RayTracingFunctions.cc, so you will have to include all three files in the same directory.