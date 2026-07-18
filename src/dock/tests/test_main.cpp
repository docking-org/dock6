// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// test_main.cpp
//
// doctest entry point for the DOCK_GA (conf_gen_ga) unit test suite.
//
// This is the ONLY translation unit that defines DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN,
// so the doctest runtime is compiled exactly once. All other tests/*.cpp files
// just `#include "doctest.h"` and register TEST_CASEs.
//
// This suite is an OPTIONAL, SEPARATE build target (`make test`) — it is not part
// of the normal `make` / `make install` path and adds no dependency to the DOCK build.
// See tests/README.md.
//
// This software is copyrighted, 2004-2026, and follows the DOCK6 license terms
// (see the repository LICENSE file). doctest.h is vendored under the MIT license
// (Copyright (c) 2016-2023 Viktor Kirilov).
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
