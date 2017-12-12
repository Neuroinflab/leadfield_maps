gmsh -3 -optimize_netgen setup_refined.geo
dolfin-convert setup_refined.msh mesh_setup.xml

