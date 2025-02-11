// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#include <hoomd/VectorMath.h>
#include <hoomd/hpmc/PairPotential.h>
#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace hpmc
    {

/** Chiral pair potential for use with HPMC simulations.
 */
class ChiralPairPotential : public PairPotential
    {
    public:
    ChiralPairPotential(std::shared_ptr<SystemDefinition> sysdef);
    virtual ~ChiralPairPotential() { }

    virtual LongReal energy(const LongReal r_squared,
                            const vec3<LongReal>& r_ij,
                            const unsigned int type_i,
                            const quat<LongReal>& q_i,
                            const LongReal charge_i,
                            const unsigned int type_j,
                            const quat<LongReal>& q_j,
                            const LongReal charge_j) const;

    /// Compute the non-additive cuttoff radius
    LongReal computeRCutNonAdditive(unsigned int type_i, unsigned int type_j) const;

    /// Set type-pair-dependent parameters to the potential.
    void setParamsPython(pybind11::tuple particle_types, pybind11::dict params);

    /// Get type-pair-dependent parameters.
    pybind11::dict getParamsPython(pybind11::tuple particle_types);

    void setMode(const std::string& mode_str)
        {
        if (mode_str == "none")
            {
            m_mode = none;
            }
        else if (mode_str == "cubic")
            {
            m_mode = cubic;
            }
        else
            {
            throw std::domain_error("Invalid mode " + mode_str);
            }
        }

    std::string getMode()
        {
        std::string result = "none";

        if (m_mode == cubic)
            {
            result = "cubic";
            }

        return result;
        }

    protected:
    /// Shifting modes that can be applied to the energy
    enum SymmetryMode
        {
        none = 0,
        cubic
        };

    protected:
    /// per-type-pair parameters
    struct ParamType
        {
        ParamType() { }

        /// Construct a parameter set from a dictionary.
        ParamType(pybind11::dict v);

        /// Convert a parameter set to a dictionary.
        pybind11::dict asDict();
        
        LongReal m_epsilon;
        LongReal m_alpha;
        LongReal m_theta;
        LongReal m_r_cut;
        };

    /// Parameters per type pair.
    std::vector<ParamType> m_params;

    SymmetryMode m_mode = none;
    };

namespace detail
    {
//! Export the ExampleUpdater class to python
void export_ChiralPairPotential(pybind11::module& m);

    } // end namespace detail

    } // end namespace hpmc
    } // end namespace hoomd
