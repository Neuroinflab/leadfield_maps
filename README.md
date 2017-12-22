LEADFIELD ATLAS
===============

# Installation
## Docker

sudo apt-get install docker docker.io
[Reference](https://fenics.readthedocs.io/projects/containers/en/latest/)

### First go
docker run -ti -v $(pwd):/home/fenics/shared --name fenics-container quay.io/fenicsproject/stable:2016.2.0

### Subsequent access
docker start fenics-container
docker exec -ti -u fenics fenics-container /bin/bash -l

### To stop
docker stop fenics-container

## List containers
docker ps -a

### Remove a container
docker rm -f CONTAINER

## Build from source
Required to access the barycenters of the meshes.
Critical step - heavily unstable and unreliable.
Only used to obtain the bary points at this stage.

See ./build_script.sh for full details on building

This can be also avoided by moving def compute_barycenters to
elsewhere [TODO - not crucial]


# Mesh generation
This step requires gmsh and dolfin. - do not have to
work together Can be skipped if the .xml files are already provided

## ITKSnap to stl
to export the rough surface of the brain

## Meshlab
stl smoothing and surface fixing

## blender
boolean subtracting of the ground electrode.

## GMSH
version 3.0+ (installed via sudo apt-get install) see
/mesh/mesh_it.sh

Use setup_refined.geo to obtain setup_refined.msh dolfin then converts
it setup_refined.msh into mesh_setup*.xml

# Probe points

This requires vtk 6.3.0 version installed correctly
access to numpy and nifti file This ALSO requires the right fenics
version installed along side. Hence the need to mesh.

See ./probe_points.py

Goes through the meshes, fetches the barycenters, goes though the
imgage files and fills the necessary fields corresponding to those
barycenters and dumps these files inside /points folder

Can skip this step if the points are already provided.

# Conductivity

This requires the output from the probepoinst to be in place
Can be avoided if simgas are already provided.

See ./conductivity.py

Populates ./sigmas/*/*.xml.gz Which correspond to the sigmas at the
mesh barycentes populated via conductivity_c.py file C-code

# Forward method

Currently works with point sources sent in as a list of points
corresponding to the nifti world coordinates. Specifying the electrode
positions is not implemented. Currently dumps as hdf5 files for the
phi vector of the potential space.

Critical dependency is with HDF5 - works smooth with docker.

Works only with hypre_amg preconditioner via PETSCL with hypre
libraries - exists in docker image luckily.

