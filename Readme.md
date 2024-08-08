## Safety and Efficiency in Robotics: the Control Barrier Functions Approach

Source code of the simulations described in the Part I of a tutorial paper on the use of Control Barrier Functions (CBFs) in the design of safety controllers for dynamic and robotic systems.

The code implements in C++ the setup of a Quadratic Programming (QP) problem embedding the safety contraints for the system under control, formalized as CBFs. The solution of the QP problem corresponds to the control value to be applied to the system to preserve its safety, while minimizing at the same time the difference between the applied control and a nominal, potentially unsafe, control input. In practice, the CBF-based optimal controller acts a safety filter for the nominal controller, whatever is the behavior of the latter.

The code included in this capsule relies on OSQP as an optimization solver, but provides a wrapping layer that can be easily extended to support other QP solvers.

## Requirements
The code is designed to be compiled on almost any platform (Linux, Windows, Mac) with either GCC or Microsoft Visual C++, provided that the following tools and libraries are installed:
* CMake
* Eigen 3 (only for the PUMA 560 kinematics included by the manipulator example)
* Python 3 and Matplotlib (only for plotting with the included version of matplotlib-cpp, adapted from https://github.com/lava/matplotlib-cpp)
* OSQP (installed from sources https://osqp.org/docs/get_started/sources.html with the option ENABLE_MKL_PARDISO=OFF)

## Examples included
* Single integrator with state space bounds
* Double integrator with state space bounds
* 2-DOF single integrator avoding an obstacle
* 6-DOF manipulator (PUMA 560) avoiding an obstacle

## Instructions
This was clone from so that i could run it locally and i was running this on my mac 
so build the image like this remove the --platform if you are on linux
```
   docker build --platform linux/amd64 -t cbf_tutorial_image_2 .
```

to run the actual script 

```
    docker run --rm -v $(pwd)/results:/results cbf_tutorial_image_2 /bin/bash -c "/code/run_script"
```
The "Reproducbile Run" command executes the "run" file in the code directory. The example executed by default is the single integrator with state space bounds. Edit the capsule and modify the "run" file by commenting/uncommenting its last lines to execute other examples.

Each example generates some plots in PNG format and a CSV file with simulation results to be download for further analysis with other tools (e.g. Matlab). Details about plotted and logged data can be found in the source file of each example. 
