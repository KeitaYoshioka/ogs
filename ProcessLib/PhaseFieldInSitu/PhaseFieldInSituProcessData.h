/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include <Eigen/Eigen>

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace PhaseFieldInSitu
{
template <int DisplacementDim>
struct PhaseFieldInSituProcessData
{
    PhaseFieldInSituProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<
                          DisplacementDim>>>&& solid_materials_,
        ParameterLib::Parameter<double> const& residual_stiffness_,
        ParameterLib::Parameter<double> const& crack_resistance_,
        ParameterLib::Parameter<double> const& crack_length_scale_,
        ParameterLib::Parameter<double> const& solid_density_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_,
        bool const propagating_crack_, int const split_method_,
        bool const crack_pressure_, double const pf_irrv_, int const at_param_)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          solid_density(solid_density_),
          specific_body_force(specific_body_force_),
          propagating_crack(propagating_crack_),
          crack_pressure(crack_pressure_),
          pf_irrv(pf_irrv_),
          at_param(at_param_),
          split_method(split_method_)
    {
    }

    PhaseFieldInSituProcessData(PhaseFieldInSituProcessData&& other) = default;

    //! Copies are forbidden.
    PhaseFieldInSituProcessData(PhaseFieldInSituProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldInSituProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldInSituProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;

    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& residual_stiffness;
    ParameterLib::Parameter<double> const& crack_resistance;
    ParameterLib::Parameter<double> const& crack_length_scale;
    ParameterLib::Parameter<double> const& solid_density;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0.0;
    double t = 0.0;
    double const unity_pressure = 1.0;
    double pressure = 0.0;
    double pressure_old = 0.0;
    double pressure_error = 0.0;
    double injected_volume = 0.0;
    double crack_volume0 = 0.0;
    double crack_volume1 = 0.0;
    double elastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
    bool propagating_crack = false;
    bool crack_pressure = false;
    double pf_irrv = 0.05;
    int at_param = 2;
    int split_method = 0;
};

}  // namespace PhaseFieldInSitu
}  // namespace ProcessLib
