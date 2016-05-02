sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
sudo apt-get install -y cmake git gcc-5 gfortran-5 g++-5 liblapack-dev netcdf-bin libnetcdf-dev

cat >> ~/compilers.sh <<EOF
CC=gcc-5
CXX=g++-5
FC=gfortran-5
export CC CXX FC
EOF
