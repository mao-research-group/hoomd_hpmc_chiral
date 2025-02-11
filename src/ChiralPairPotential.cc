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

std::vector<rotmat3<LongReal>> CubicSymmetries(){
    std::vector<rotmat3<LongReal>> smat_list;
    smat_list.resize(24);  // all initialized as identity matrices
    
    std::vector<int> signs = {1, -1};
    int idx = 0;
    for (int row0=0; row0<3; row0++){
        for (auto sign0 : signs){
            for (auto sign1 : signs){
                for (auto sign2 : signs){
                    int row1 = 0;
                    int row2 = 0;
                    if ((sign0*sign1*sign2)>0){
                        row1 = (row0+1)%3;
                        row2 = (row1+1)%3;
                    }
                    else {
                        row1 = (row0+2)%3;
                        row2 = (row1+2)%3;
                    }
                    smat_list[idx].row0[0] = 0;
                    smat_list[idx].row1[1] = 0;
                    smat_list[idx].row2[2] = 0;
                    smat_list[idx].row0[row0] = sign0;
                    smat_list[idx].row1[row1] = sign1;
                    smat_list[idx].row2[row2] = sign2;

                    idx = idx+1;
                }
            }
        }
    }
    
    return smat_list;
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

    rotmat3 rmati(q_i);
    rotmat3 rmatj(q_j);
    rotmat3 rmatj_inv = transpose(rmatj);

    LongReal ori_factor;
    LongReal current_lowest = 8; // maximum value of ori before shift, ~5.9 for cubic and 8 for no symmetries
    vec3<LongReal> rhat_ij = r_ij / fast::sqrt(r_squared);

    std::vector<rotmat3<LongReal>> smat_list;
    if (m_mode == cubic){
        smat_list = CubicSymmetries();
    } else {
        smat_list.resize(1);  // identity matrix by default
    }
    size_t symm_size = smat_list.size();
    for (size_t i = 0; i < symm_size; i++){

        rotmat3 rmat_times = rmati * smat_list[i];
        rmat_times = rmat_times * rmatj_inv;

        rotmat3 rmat_bet = rmati.fromAxisAngle(rhat_ij, param.m_theta);   

        vec3<LongReal> row0 = rmat_times.row0 - rmat_bet.row0;
        vec3<LongReal> row1 = rmat_times.row1 - rmat_bet.row1;
        vec3<LongReal> row2 = rmat_times.row2 - rmat_bet.row2;
        ori_factor = dot(row0, row0) + dot(row1, row1) + dot(row2, row2);
        if (ori_factor<current_lowest){
            current_lowest = ori_factor;
        }
    }
    ori_factor = param.m_alpha - current_lowest; // Range moved from 0:max to alpha:alpha-max
    LongReal energy = param.m_epsilon * ori_factor; // epsilon should be negative usually

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
