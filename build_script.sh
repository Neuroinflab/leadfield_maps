wget http://repo.continuum.io/archive/Anaconda2-4.3.0-Linux-x86_64.sh
bash Anaconda2-4.3.0-Linux-x86_64.sh
conda create --name fenics2016_2 --file fenics_2016_2.txt
source activate fenics2016_2


sudo apt-get install emacs git libglu1-mesa libxi-dev libxmu-dev libglu1-mesa-dev libboost-dev libboost-filesystem-dev libboost-all-dev cmake libeigen3-dev swig libparmetis-dev scotch ptscotch fcc

mkdir fenics
cd fenics
git clone https://bitbucket.org/fenics-project/ffc.git
git clone https://bitbucket.org/fenics-project/fiat.git
git clone https://bitbucket.org/fenics-project/dijitso.git
git clone https://bitbucket.org/fenics-project/dolfin.git
git clone https://bitbucket.org/fenics-project/instant.git
git clone https://bitbucket.org/fenics-project/ufl.git
git clone https://bitbucket.org/fenics-project/fenics.git
git clone https://bitbucket.org/fenics-project/mshr.git

git clone -b maint https://bitbucket.org/petsc/petsc petsc
git clone https://bitbucket.org/slepc/slepc


cd petsc
git checkout tags/v3.7
# ./configure --download-f2cblaslapack
./configure --download-superlu_dist --download-mumps --download-hypre --download-f2cblaslapack --download-scalapack --download-parmetis --download-metis
make all test
export PETSC_DIR=$HOME/fenics/petsc
export PETSC_ARCH=arch-linux2-c-debug
cd ..

cd slepc
git checkout tags/v3.7
./configure
make SLEPC_DIR=$PWD
make test SLEPC_DIR=$PWD
export SLEPC_DIR=$HOME/fenics/slepc
cd ..


cd instant
git checkout tags/instant-2016.2.0
pip install .
cd ..

cd ufl
git checkout tags/ufl-2016.2.0
pip install .
cd ..

cd dijitso
git checkout tags/dijitso-2016.2.0
pip install .
cd ..

cd fiat
git checkout tags/fiat-2016.2.0
pip install .
cd ..

cd ffc
git checkout tags/ffc-2016.2.0
pip install .
cd ..

cd dolfin
git checkout tags/dolfin-2016.2.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/anaconda2/envs/fenics2016_2 ..
sed -i 's/lib64/lib\/x86_64-linux-gnu/g' dolfin/CMakeFiles/dolfin.dir/build.make
# search for string /usr/lib64/ and replace with /usr/lib/x86_64-linux-gnu/
make install
cd ..
