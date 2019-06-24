/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include <memory>
#include <utility>

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

namespace PhaseFieldExternal
{
template <int DisplacementDim>
struct PhaseFieldExternalProcessData
{
    PhaseFieldExternalProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<
                          DisplacementDim>>>&& solid_materials_,
        ParameterLib::Parameter<double> const& residual_stiffness_,
        ParameterLib::Parameter<double> const& crack_resistance_,
        ParameterLib::Parameter<double> const& crack_length_scale_,
        ParameterLib::Parameter<double> const& solid_density_,
        ParameterLib::Parameter<double> const&
            linear_thermal_expansion_coefficient_,
        ParameterLib::Parameter<double> const& pressure_ext_,
        ParameterLib::Parameter<double> const& temperature_ext_,
        ParameterLib::Parameter<double> const& biot_coefficient_,
        double const reference_temperature_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_,
        int const split_method_,double const reg_param_, double const pf_irrv_, int const at_param_)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          solid_density(solid_density_),
          linear_thermal_expansion_coefficient(
              linear_thermal_expansion_coefficient_),
          pressure_ext(pressure_ext_),
          temperature_ext(temperature_ext_),
          biot_coefficient(biot_coefficient_),
          reference_temperature(reference_temperature_),
          specific_body_force(specific_body_force_),
          split_method(split_method_),
          reg_param(reg_param_),
          pf_irrv(pf_irrv_),
          at_param(at_param_)
    {
    }

    PhaseFieldExternalProcessData(PhaseFieldExternalProcessData&& other) =
        default;

    //! Copies are forbidden.
    PhaseFieldExternalProcessData(PhaseFieldExternalProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldExternalProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldExternalProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;

    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& residual_stiffness;
    ParameterLib::Parameter<double> const& crack_resistance;
    ParameterLib::Parameter<double> const& crack_length_scale;
    ParameterLib::Parameter<double> const& solid_density;
    ParameterLib::Parameter<double> const& linear_thermal_expansion_coefficient;
    ParameterLib::Parameter<double> const& pressure_ext;
    ParameterLib::Parameter<double> const& temperature_ext;
    ParameterLib::Parameter<double> const& biot_coefficient;
    double const reference_temperature;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt;
    double t;
    int split_method = 0;
    double reg_param = 0.01;
    double pf_irrv = 0.05;
    int at_param = 2;
};

}  // namespace PhaseFieldExternal
}  // namespace ProcessLib
