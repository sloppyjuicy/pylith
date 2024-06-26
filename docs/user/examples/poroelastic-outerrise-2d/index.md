# 2D Outer-Rise Hydration Using Poroelasticity

The files are in the directory `examples/poroelastic-outerrise-2d`.
The files and directories for this set of examples includes:

:`README.md`: README file containing a brief description of the various examples.
:`*.cfg`: PyLith parameter files.
:`generate_gmsh.py`: Python script to generate mesh using Gmsh.
:`*.msh`: Gmsh mesh files generated by Gmsh.
:`*.spatialdb`: Spatial database filesFiles associated with the spatial databases.
:`output`: Directory containing simulation output. It is created automatically when running the simulations.

## Overview

This example examines how bending stresses and permeability structure influence the hydration state of a slab of oceanic lithosphere in the outer-rise of a subduction zone.
The poroelastic slab is deformed based on bathymetry along a trench-orthogonal transect of the Central-American subduction system.
Enhanced permeability is introduced in the outer-rise to simulate the effect of normal faulting.

:Step 1: Underformed lithosphere with a uniform fluid-pressure applied to the top-boundary
:Step 2: Lithosphere is deformed to the current bathymetry over 0.3 Myr. 
:Step 3: Lithosphere is deformed to the current bathymetry over 0.3 Myr with enhanced permeability from outer-rise faults. 

:::{figure-md} fig:example:poroelastic:outerrise:2d:overview
<img src="figs/geometry.*" alt="Diagram of geometry for poroelastic-outerrise-2d." width="75%"/>

Diagram of geometry used for poroelastic-outerrise-2d.
The domain extends from 0 km to +150 km in the x direction and from -30 km to 0 km in the y direction.
We refer to the domain boundaries using the names shown in the diagram.
:::

## Example Workflow

:::{toctree}
meshing-gmsh.md
common-information.md
step01-no-faults-no-flexure.md
step02-no-faults-flexure.md
step03-faults-flexure.md
exercises.md
:::
