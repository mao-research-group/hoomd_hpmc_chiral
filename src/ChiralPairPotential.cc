// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "ChiralPairPotential.h"

namespace hoomd
    {
namespace hpmc
    {

ChiralPairPotential::ChiralPairPotential(std::shared_ptr<SystemDefinition> sysdef)
    : PairPotential(sysdef), m_params(m_type_param_index.getNumElements())
    {
    }

std::vector<quat<LongReal>> CubicSymmetries(){
    std::vector<quat<LongReal>> squat_list;
    squat_list.resize(24);  
    squat_list[0].s = 0;
    squat_list[0].v = vec3<LongReal>(1, 0, 0);
    
    return squat_list;
}

LongReal ChiralPairPotential::energy(const LongReal r_squared,
                                      const vec3<LongReal>& r_ij,
                                      const unsigned int type_i,
                                      const quat<LongReal>& q_i,
                                      const LongReal charge_i,
                                      const unsigned int type_j,
                                      const quat<LongReal>& q_j,
                                      const LongReal charge_j) const
    {
    unsigned int param_index = m_type_param_index(type_i, type_j);
    const auto& param = m_params[param_index];

    quat qj_inv = conj(q_j);
    vec3<LongReal> rhat_ij = r_ij / fast::sqrt(r_squared);
    quat qaa_inv = quat<LongReal>::fromAxisAngle(rhat_ij, -param.m_theta);
    quat q_rel = q_i*qj_inv*qaa_inv;

    LongReal ori_factor;
    LongReal current_highest = 0; // minimum value of q_s^2

    std::vector<quat<LongReal>> squat_list;
    if (m_mode == cubic){
        squat_list = CubicSymmetries();
    } else {
        squat_list.resize(1);  // identity quat by default
    }

    size_t symm_size = squat_list.size();
    for (size_t i = 0; i < symm_size; i++){
        quat sym_q_rel = squat_list[i] * q_rel;
        ori_factor = sym_q_rel.s * sym_q_rel.s;
        if (ori_factor>current_highest){
            current_highest = ori_factor;
        }
    }
    // 4(1-cos(ang)) = 4(1-(2cos^2(ang/2)-1)) = 8(1-q_s^2)
    current_highest = LongReal(8.0) * (LongReal(1.0) - current_highest);
    LongReal shifted_factor = param.m_alpha - current_highest; // Range moved from 0:max to alpha:alpha-max
    LongReal energy = param.m_epsilon * shifted_factor; // epsilon should be negative usually

    return energy;
    }

LongReal ChiralPairPotential::computeRCutNonAdditive(unsigned int type_i,
                                                      unsigned int type_j) const
    {
    unsigned int param_index = m_type_param_index(type_i, type_j);
    return m_params[param_index].m_r_cut;
    }

void ChiralPairPotential::setParamsPython(pybind11::tuple particle_types, pybind11::dict params)
    {
    auto pdata = m_sysdef->getParticleData();
    auto type_i = pdata->getTypeByName(particle_types[0].cast<std::string>());
    auto type_j = pdata->getTypeByName(particle_types[1].cast<std::string>());
    unsigned int param_index_1 = m_type_param_index(type_i, type_j);
    m_params[param_index_1] = ParamType(params);
    unsigned int param_index_2 = m_type_param_index(type_j, type_i);
    m_params[param_index_2] = ParamType(params);

    notifyRCutChanged();
    }

pybind11::dict ChiralPairPotential::getParamsPython(pybind11::tuple particle_types)
    {
    auto pdata = m_sysdef->getParticleData();
    auto type_i = pdata->getTypeByName(particle_types[0].cast<std::string>());
    auto type_j = pdata->getTypeByName(particle_types[1].cast<std::string>());
    unsigned int param_index = m_type_param_index(type_i, type_j);
    return m_params[param_index].asDict();
    }

ChiralPairPotential::ParamType::ParamType(pybind11::dict params)
    {
    // TODO: unpack per-type-pair quanties from the Python dictionary to the ParamType struct.
    m_epsilon = params["epsilon"].cast<LongReal>();
    m_alpha = params["alpha"].cast<LongReal>();
    m_theta = params["theta"].cast<LongReal>();
    m_r_cut = params["r_cut"].cast<LongReal>();
    }

pybind11::dict ChiralPairPotential::ParamType::asDict()
    {
    pybind11::dict pydict;
    // TODO; pack per-type-pair quantities from the ParamType struct to the Python dictionary.
    pydict["epsilon"] = m_epsilon;
    pydict["alpha"] = m_alpha;
    pydict["theta"] = m_theta;
    pydict["r_cut"] = m_r_cut;
    return pydict;
    }

namespace detail
    {
void export_ChiralPairPotential(pybind11::module& m)
    {
    pybind11::class_<ChiralPairPotential, PairPotential, std::shared_ptr<ChiralPairPotential>>(
        m,
        "ChiralPairPotential")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>>())
        .def("setParams", &ChiralPairPotential::setParamsPython)
        .def("getParams", &ChiralPairPotential::getParamsPython)
        .def_property("mode",
            &ChiralPairPotential::getMode,
            &ChiralPairPotential::setMode);
    }
    } // end namespace detail
    } // end namespace hpmc
    } // end namespace hoomd
