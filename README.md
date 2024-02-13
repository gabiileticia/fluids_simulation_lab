# fluids_simulation_lab
Project for implementing SPH using different pressure solvers and surface tension. Additionally surface reconstruction using Marching Cubes algorithm was also implemented.

## Use
### Install/Compile
Need CMake and modern GlibC

First use CMake to create makefiles with

```
cmake CMakeList.txt
```

When done recommend changing compile mode in CmakeCache.txt to Release.

Then compile executables with ```make```

### What now?

Execute main simulator by invoking ```./app/app_learnSPH <Simulation Number>```. 
The simulation number is simply the number that is mapped to a simulation setup function in the switch-case at the beginning of main.cpp.
If you want to create your own simulation scenarios simply create a new function in simulation_setup.cpp and .h and then add it in the switch-case of main.cpp.

### fast clearing of output directory
```yes | rm res/wcsph/*.vtk ```