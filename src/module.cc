// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include <pybind11/pybind11.h>
#include "ChiralPairPotential.h"

using namespace hoomd::hpmc::detail;

PYBIND11_MODULE(_hpmc_chiral, m)
    {
    export_ChiralPairPotential(m);

#ifdef ENABLE_HIP
    // TODO: Call export_ClassGPU(m) for each GPU enabled C++ class to be exported
    // to Python.
#endif
    }
