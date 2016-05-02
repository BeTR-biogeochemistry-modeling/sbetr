sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
sudo apt-get install -y cmake gcc-5 liblapack-dev gfortran-5 git netcdf-bin libnetcdf-dev

ls -al /usr/bin/gfor*
