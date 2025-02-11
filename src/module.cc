// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include <pybind11/pybind11.h>
// TODO: Include the header files of classes that will be exported to Python.
#include "ChiralPairPotential.h"

using namespace hoomd::hpmc::detail;

// TODO: Set the name of the python module to match ${COMPONENT_NAME} (set in
// CMakeLists.txt), prefixed with an underscore.
PYBIND11_MODULE(_hpmc_chiral, m)
    {
    // TODO: Call export_Class(m) for each C++ class to be exported to Python.
    export_ChiralPairPotential(m);

#ifdef ENABLE_HIP
    // TODO: Call export_ClassGPU(m) for each GPU enabled C++ class to be exported
    // to Python.
#endif
    }
