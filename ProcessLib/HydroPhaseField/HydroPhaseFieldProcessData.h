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

namespace HydroPhaseField
{
template <int DisplacementDim>
struct HydroPhaseFieldProcessData
{
    HydroPhaseFieldProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& residual_stiffness_,
        Parameter<double> const& crack_resistance_,
        Parameter<double> const& crack_length_scale_,
        Parameter<double> const& kinetic_coefficient_,
        Parameter<double> const& solid_density_,
        Parameter<double>& history_field_,
        Parameter<double> const& flow_rate_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_)
        : material{std::move(material_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          kinetic_coefficient(kinetic_coefficient_),
          solid_density(solid_density_),
          history_field(history_field_),
          specific_body_force(specific_body_force_),
          flow_rate(flow_rate_)
    {
    }

    HydroPhaseFieldProcessData(HydroPhaseFieldProcessData&& other)
        : material{std::move(other.material)},
          residual_stiffness(other.residual_stiffness),
          crack_resistance(other.crack_resistance),
          crack_length_scale(other.crack_length_scale),
          kinetic_coefficient(other.kinetic_coefficient),
          solid_density(other.solid_density),
          history_field(other.history_field),
          specific_body_force(other.specific_body_force),
          flow_rate(other.flow_rate),
          dt(other.dt),
          t(other.t)
    {
    }
    s
    //! Copies are forbidden.
    HydroPhaseFieldProcessData(HydroPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroPhaseFieldProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    Parameter<double> const& residual_stiffness;
    Parameter<double> const& crack_resistance;
    Parameter<double> const& crack_length_scale;
    Parameter<double> const& kinetic_coefficient;
    Parameter<double> const& solid_density;
    Parameter<double>& history_field;
    Parameter<double>& flow_rate;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0.0;
    double t = 0.0;
    double injected_volume = 0.0;
    double unity_pressure = 1.0;
};

}  // namespace HydroPhaseField
}  // namespace ProcessLib
