#!/usr/bin/env bash

# create boundary vtu's
~/w/ogs/r/bin/constructMeshesFromGeometry -g geometry.gml -m mesh_rev1.vtu -s 1e-8
rename bar_ mesh_rev1_ bar_*.vtu

# create METIS mesh
~/w/ogs/r/bin/partmesh -s -i mesh_rev1.vtu

# partition the main mesh and the boundary meshes
~/w/ogs/r/bin/partmesh -n 2 -m -i mesh_rev1.vtu -- mesh_rev1_*vtu
