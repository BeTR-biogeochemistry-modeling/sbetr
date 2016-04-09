SBETR = ../build/Darwin-x86_64-static-double-cc-Debug/src/driver/sbetr

RTEST = ./rtest_betr.py
RTEST_ARGS = --backtrace
RTEST_ARGS += --executable $(SBETR)

COVERAGERC = --rcfile mtest/.coveragerc
COVERAGE_RUN_FLAGS = $(COVERAGERC) -m unittest discover

test rtest : FORCE
	$(RTEST) $(RTEST_ARGS)

# NOTE(bja, 201603) rcheck - runs the regression test manager with
# --check-only, to recheck existing results. It does NOT run
# sbetr. Useful for development of the regression test manager itself.
rcheck : FORCE
	$(RTEST) $(RTEST_ARGS) --check-only

# NOTE(bja, 201603) run the meta-tests for the test manager.
mtest : FORCE
	python -m unittest discover --buffer --start-directory mtest/unit/

test-coverage : FORCE
	coverage erase
	COVERAGE_FILE=.coverage.unit coverage run $(COVERAGE_RUN_FLAGS) mtest/unit
	COVERAGE_FILE=.coverage.xfail coverage run $(COVERAGE_RUN_FLAGS) mtest/xfail
	coverage combine $(COVERAGERC)
	coverage html $(COVERAGERC)
	coverage report $(COVERAGERC)
	@echo 'For detailed analysis, see the html report at:'
	@echo '    coverage-html/index.html'

clean : FORCE
	-find . -name '*~' -print0 | xargs -0 rm


# try to remove all generated files. Note, these globs should be as
# specific as possible to avoid removing unintended files!
clobber : clean
	-find . -name '*.cdl.nc' -print0 | xargs -0 rm
	-find . -name '*.output.nc' -print0 | xargs -0 rm
	-find . -name '*.bak' -print0 | xargs -0 rm
	-find . -name '*.stdout' -print0 | xargs -0 rm
	-find . -name '*.regression' -print0 | xargs -0 rm


FORCE : 
