SBETR = ../build/Darwin-x86_64-static-double-cc-Debug/src/driver/sbetr

RTEST = ./rtest_betr.py
RTEST_ARGS = --backtrace
RTEST_ARGS += --executable $(SBETR)

rtest : FORCE
	$(RTEST) $(RTEST_ARGS)

# NOTE(bja, 201603) rcheck - runs the regression test manager with
# --check-only, to recheck existing results. It does NOT run
# sbetr. Useful for development of the regression test manager itself.
rcheck : FORCE
	$(RTEST) $(RTEST_ARGS) --check-only

# NOTE(bja, 201603) run the meta-tests for the test manager.
mtest : FORCE
	python -m unittest discover --buffer --start-directory mtest/unit/

clean : FORCE
	-rm *~

FORCE : 
