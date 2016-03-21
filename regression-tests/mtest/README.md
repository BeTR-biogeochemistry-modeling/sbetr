# BeTR regression test "meta-tests"

The test system needs to be the most bomb proof piece of code we
have. The 'mtest' directory contains the tests for the test system
itself.

## Running tests

If you make changes to the testing infracture, you *MUST* run the meta-test suite:

    # from the regression-tests directory, not mtest!
    make mtest
    python -m unittest discover --buffer --start-directory mtest/unit/


## Test Design

The testing infrastructure is designed to be heavily tested. Key principles are:

* read data into memory as soon as possible and isolate calles to the
  file system. Pass in-memory data objects to data processing routines
  to avoid working with files that are hard to test.

* no global variables

* no left justified code. all code must be in a subroutine so it can
  be isolated and tested, including main.

There are two types of tests for the test system, unit tests and
system tests. all are launched through a unit test framework.

The system tests are primarily expected failure tests. We seed the
test system with bad data and ensure that it catches it. system tests
should be run in 'check-only' mode when ever possible, or be extremely
short model runs.

## Code Coverage

A key aspect of being confident in the meta-test system is code
coverage. If you make changes to the test infrastructure, you must run
the code coverage tool and ensure that you are maintaining or
increasing code coverage before your changes will be accepted. Code
coverage reports are obtained by running:

    coverage
