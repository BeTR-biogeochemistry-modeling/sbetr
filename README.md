# BeTR

BeTR is a standalone reactive transport libary designed to be
integrated ino land surface models such as CLM and ALM.

Jinyun Tang, jinyuntang@lbl.gov

## Building

There is a cmake based build system that has been tested on Mac OS X
with gfortran installed into the users path. To build:

    cd ${SBETR_ROOT}
    make config
    make debug=1 all

This will do an out of source build in:

    ${SBETR_ROOT}/build/Xyz/src

where Xyz is the build configuration.

The standalone `sbetr` executable is in:

    ${SBETR_(ROOT}/bulid/Xyz/src/driver/sbetr

## Testing

Integration of pFUnit based unit tests and a systems level regression
test driver are in progress.

## Running

sbetr takes exactly one command line arguement, the name of the input
namelist file. The namelist input file specifies runtime configuration
and paths to other input data files. NOTE: paths are relative to the
directory where sbetr is executed!


    cd ${SBETR_ROOT}/example_input
    ../build/Darwin-x86_64-static-double-cc-Debug/src/driver/sbetr mock.namelist


---------------------------------------------

README (jinyuntang@lbl.gov)
direcotries:
code/src contains model code
code/bld contains makefile and configuration scripts
scripts contains scripts to create run directories
example_input contains example input data

the example is set with mock run that transport five tracers: N2, O2, AR, CO2 and CH4
before compile the code, make sure the set up in Macros.compiler is appropriate,
then check the new_case script in directory scripts/ specifies the compiler to be used.
running ./new_case will then generate the correct Makefile for your case.
before building the code, make sure the data directory in CLMForcType.F90 is setup
correctly, right now it is hardwired, but it eventually will be replaced with
input from ascii file or namelist file.
To configure a new bgc implementation, follow the example in BGCReactionsMockRunType
and add the new configuration name to BGCReactionsFactoryMod, and set the correct
bgc_method in betr_initializeMod (because it is not read in from namelist in this offline
code, one has to change it manually).
After the above procedure, run new_case.build in script to build the code, the executalbe
sbetr.exe will be in a local directory build/, then copy it to whereever you want to
run the code. 
