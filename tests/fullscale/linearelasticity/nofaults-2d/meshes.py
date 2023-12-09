# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.testing.FullTestApp import MeshEntity

class TriGmsh(object):
    """Mesh information for tri mesh generated by Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=184, ncorners=3, nvertices=81),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "elastic_xneg": MeshEntity(ncells=64, ncorners=3, nvertices=45),
        "elastic_xpos": MeshEntity(ncells=64, ncorners=3, nvertices=45),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }


class QuadGmsh(object):
    """Mesh information for quad mesh generated by Gmsh.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=99, ncorners=4, nvertices=120),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        "elastic_xneg": MeshEntity(ncells=54, ncorners=4, nvertices=70),
        "elastic_xpos": MeshEntity(ncells=45, ncorners=4, nvertices=60),

        "bc_xneg": MeshEntity(ncells=9, ncorners=2, nvertices=10),
        "bc_xpos": MeshEntity(ncells=9, ncorners=2, nvertices=10),
        "bc_yneg": MeshEntity(ncells=11, ncorners=2, nvertices=12),
        "bc_ypos": MeshEntity(ncells=11, ncorners=2, nvertices=12),
    }


class TriCubit(object):
    """Mesh information for tri mesh generated by Cubit.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=124, ncorners=3, nvertices=79),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        # Materials
        "elastic_xneg": MeshEntity(ncells=60, ncorners=3, nvertices=43),
        "elastic_xpos": MeshEntity(ncells=64, ncorners=3, nvertices=45),

        # Boundaries
        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }


class QuadCubit(object):
    """Mesh information for quad mesh generated by Cubit.
    """
    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=4, nvertices=81),
        "points": MeshEntity(ncells=3, ncorners=1, nvertices=3),

        "elastic_xneg": MeshEntity(ncells=32, ncorners=4, nvertices=45),
        "elastic_xpos": MeshEntity(ncells=32, ncorners=4, nvertices=45),

        "bc_xneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_xpos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }



# End of file
