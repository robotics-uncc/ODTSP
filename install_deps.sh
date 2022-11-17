# Orbiting Dubins Traveling Salesman Problem
# A. Wolek, J. McMahon, Sep-2018
#
# LICENSE:
#
# The source code is in the public domain and not licensed or under
# copyright. The information and software may be used freely by the public.
# As required by 17 U.S.C. 403, third parties producing copyrighted works
# consisting predominantly of the material produced by U.S. government
# agencies must provide notice with such work(s) identifying the U.S.
# Government material incorporated and stating that such material is not
# subject to copyright protection.
#
# Derived works shall not identify themselves in a manner that implies an
# endorsement by or an affiliation with the Naval Research Laboratory.
#
# RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
# SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
# RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
# OF RECIPIENT IN THE USE OF THE SOFTWARE.
#
# (UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

#!/bin/sh
echo "***********************************************************"
echo "Starting dependency installer for ODTSP Library"
echo "  Last Updated: Sep-2018, A. Wolek"
echo "***********************************************************"

echo "==========================================================="
echo " Updating apt-get..."
echo "==========================================================="
sudo apt-get update
./clean -all
./clean -ext

echo "==========================================================="
echo " Installing g++, cmake, make..."
echo "==========================================================="
sudo apt-get install g++ cmake build-essential -y

echo "================================================================="
echo "Installing math libraries (LAPACK,BLAS,ARPACK,EIGEN,ATLAS)..."
echo "================================================================="
sudo apt-get install libopenblas-dev liblapack-dev libblas-dev libarpack++2-dev -y
sudo apt-get install libeigen3-dev
sudo apt-get install libatlas-base-dev

echo "================================================================="
echo "Installing math libraries (python)..."
echo "================================================================="
sudo apt-get install python2.7 python-pip -y
sudo apt-get install python3-pip -y
sudo apt-get install fftw3 fftw3-dev pkg-config -y

echo "================================================================="
echo "Installing google-glog + gflags..."
echo "================================================================="
sudo apt-get install libgoogle-glog-dev

echo "================================================================="
echo "Installing SuiteSparse and CXSparse..."
echo "================================================================="
sudo apt-get install libsuitesparse-dev

echo "================================================================="
echo "Installing Armadillo..."
echo "================================================================="
cd external 
tar xf armadillo-8.500.0.tar.xz
cd armadillo-8.500.0
cmake .
make 
sudo make install	
cd ..
cd ..

echo "================================================================="
echo "Installing Ceres..."
echo "================================================================="
cd external 
tar -xvzf ceres-solver.tar.gz
cd ceres-solver
mkdir ceres-bin
cd ceres-bin
cmake ..
make -j3
#make test
sudo make install
cd ..
cd ..

echo "================================================================="
echo "Installing Cliquer..."
echo "================================================================="
tar -xvzf cliquer-1.21.tar.gz
cd cliquer-1.21
make
cd ..
cd ..

echo "================================================================="
echo "Installing GLKH..."
echo "================================================================="
cd external
tar -xvf GLKH-1.0.tar
cd GLKH-1.0/
make
cd ..
cd ..

echo "================================================================="
echo "Installing GLKH..."
echo "================================================================="
cd external
tar -xvf GLKH-1.0.tar
cd GLKH-1.0/
make
cd ..
cd ..


echo "================================================================="
echo "Adding environemnt variables, symbolic links..."
echo "================================================================="
cwd=$(pwd)
echo "# PATH update, autogenerated from libODTSP install_deps.sh" >> ~/.bashrc
echo "export PATH=$PATH:${cwd}/external/GLKH-1.0:${cwd}/external/GLKH-1.0/LKH-2.0.7" >> ~/.bashrc
. ~/.bashrc
ln -s ./external/GLKH-1.0/LKH LKH


echo "================================================================="
echo "Running build script..."
echo "================================================================="
./clean.sh -all
./build.sh

echo "==========================================================="
echo " Install script complete."
echo "==========================================================="



