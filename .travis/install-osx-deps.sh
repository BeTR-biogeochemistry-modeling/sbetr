# Install required software
brew update
brew install git
brew upgrade cmake
brew tap homebrew/science
brew unlink gcc
brew install gcc netcdf --with-fortran --without-cxx

# Make sure the weird gfortran library links are in place.
ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

cat >> ~/compilers.sh <<EOF
CC=gcc
CXX=g++
FC=gfortran
export CC CXX FC
EOF
