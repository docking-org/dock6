#  Makefile for the DOCK_GA (conf_gen_ga) unit test suite.
#
#  Self-contained ON PURPOSE: it does NOT include ../../install/config.h, so the
#  test suite builds and runs even in a checkout that has not been ./configure'd.
#  This is what lets `make test` work in the serial configuration on any machine
#  with a C++ compiler, adding no dependency to the normal DOCK build.
#
#  The parent Makefile's `test` target delegates here and forwards the configured
#  compiler (CXX) so a lab build uses its own toolchain. Standalone usage:
#
#      cd src/dock/tests && make run     # build + run the suite
#      make build                        # just build ./run_tests
#      make clean
#
#  As Phase 2/3 extract low-dependency units (ga_filters.*, ga_selection.*, ...),
#  add each new tests/test_<unit>.cpp to TEST_SRCS below and, when it needs the
#  extracted object, add that object to UNIT_OBJS (built by the parent build).

CXX      ?= c++
# Kept minimal and serial: no RDKit/MPI defines, no optimization surprises.
CXXFLAGS ?= -std=c++11 -O0 -g -Wall -Wno-vla-cxx-extension

# doctest entry point + every test translation unit.
TEST_SRCS = test_main.cpp test_smoke.cpp test_ga_descriptors.cpp test_ga_filters.cpp test_ga_selection.cpp test_ga_mutation.cpp test_ga_naming.cpp

# Extracted GA units (from ../) that the tests exercise. These are self-contained
# translation units decoupled from DOCKMol/globals during Phase 2, so they compile
# straight into the test binary with no full DOCK build. Add each new extracted unit
# here as it lands.
UNIT_SRCS = ../ga_descriptors.cpp ../ga_filters.cpp ../ga_selection.cpp ../ga_mutation.cpp ../ga_naming.cpp

BIN = run_tests

.PHONY: run build clean check-compile
run: build
	./$(BIN)

build: $(BIN)

# Fast syntax-only regression gate for conf_gen_ga.cpp (see compile_check.sh).
# Runnable even where a full dock6 build is impossible (no Fortran/RDKit toolchain).
check-compile:
	CXX="$(CXX)" ./compile_check.sh

$(BIN): $(TEST_SRCS) $(UNIT_SRCS) doctest.h
	$(CXX) $(CXXFLAGS) -I. -I.. -o $(BIN) $(TEST_SRCS) $(UNIT_SRCS)

clean:
	/bin/rm -f $(BIN)
	/bin/rm -rf $(BIN).dSYM
