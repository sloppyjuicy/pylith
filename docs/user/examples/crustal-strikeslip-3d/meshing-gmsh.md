# Gmsh Mesh

## Geometry

We use a the OpenCASCADE geometry engine to create a rectangular prism for the domain and construct vertical fault surfaces from the surface traces.

## Meshing using Python Script

We use the Python script `generate_gmsh.py` to create the geometry and generate the mesh.
The script makes use of the `gmsh_utils.GenerateMesh` class (discussed in {ref}`sec-user-meshing-gmsh-utils`), which provides the command line arguments and boilerplate methods.
In our `generate_gmsh.py` Python script, we create a class `App` that implements the functionality missing in `gmsh_utils.GenerateMesh`.
We must implement the `create_geometry()`, `mark()`, and `generate_mesh()` methods that are abstract in the `GenerateMesh` base class.

:::{figure-md} fig:example:crustal:strikeslip:3d:gmsh:geometry
<img src="figs/gmsh-geometry.*" alt="Geometry created in Gmsh for generating the mesh." width="500px"/>

Geometry created in Gmsh for generating the finite-element mesh.
:::
tro
:::{important}
In most cases, we store geometric entities in variables, so that we can refer to them later.
However, when we subdivide the domain using the fault surfaces, we do not know the ids of the newly created entities.
In these cases we use the GUI to identify the curves and surfaces and save the ids in variables.
:::

:::{important}
In order for Gmsh to create cells with faces that align with the fault surfaces, we embed the fault surfaces in the domain.
We use the Gmsh OpenCASCADE geometry function `fragment()` to embed the fault surfaces.
:::

We use the Gmsh MeshSize options to define a discretization size that grows slowly at a geometric rate with distance from the fault.
See [6.3.1 Specifying mesh element sizes in the Gmsh documentation](https://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes) for more information.

```{code-block} console
---
caption: Run the `generate_gmsh.py` Python script to generate the mesh.
---
# Generate a mesh with tetrahedral cells and save it to `mesh_tet.msh` (default filename).
$ ./generate_gmsh.py --write

# Save as above but start the Gmsh graphical interface after saving the mesh.
$ ./generate_gmsh.py --write --gui

# Create only the geometry and start the Gmsh graphical interface.
$ ./generate_gmsh.py --geometry --gui

# Show available command line arguments.
$ ./generate_gmsh.py --help
```

By default the Python script will generate a finite-element mesh with tetrahedral cells and save it to the file `mesh_tet.msh`.
You can view the mesh using Gmsh either by using the `--gui` command line argument when you generate the mesh or running Gmsh from the command line and opening the file.

```{code-block} console
---
caption: View the Gmsh mesh file `mesh_tet.msh` using Gmsh.
---
gmsh -open mesh_tet.msh
```

:::{figure-md} fig:example:crustal:strikeslip:3d:gmsh:mesh
<img src="figs/gmsh-tet.*" alt="Finite-element mesh with tetrahedral cells generated by Gmsh." width="500px"/>

Finite-element mesh with tetrahedral cells generated by Gmsh.
:::