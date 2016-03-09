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

BeTR testing inclueds [pFUnit](http://pfunit.sourceforge.net/) based
unit tests and a systems level regression test driver.

pFUnit tests are created automatical during the build. Run the unit
tests with

    make test

or by calling ctest in the build directory.

Regression tests are based on calling the standalone sbetr executable
and checking the results are within a specified epsilon of a baseline.


### Creating new tests

*Document procedure for new pFUnit tests*

## Running

sbetr takes exactly one command line arguement, the name of the input
namelist file. The namelist input file specifies runtime configuration
and paths to other input data files. NOTE: paths are relative to the
directory where sbetr is executed!


    cd ${SBETR_ROOT}/example_input
    ../build/Darwin-x86_64-static-double-cc-Debug/src/driver/sbetr mock.namelist

The example is set with mock run that transport five tracers: N2, O2,
AR, CO2 and CH4



## Development

Key direcotries:
 * 3rd-party - select 3rd-party sources that betr depends on.
 * src - contains model code
 * cmake - contains utilities for the configuration and build system.
 * regression_tests - regressiont test input and baselines.


To configure a new bgc implementation, follow the example in
BGCReactionsMockRunType and add the new configuration name to
BGCReactionsFactoryMod.

