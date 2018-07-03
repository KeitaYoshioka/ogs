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

namespace HydroMechanicalPhaseField
{
template <int DisplacementDim>
struct HydroMechanicalPhaseFieldProcessData
{
    HydroMechanicalPhaseFieldProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& residual_stiffness_,
        Parameter<double> const& crack_resistance_,
        Parameter<double> const& crack_length_scale_,
        Parameter<double> const& solid_density_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_,
        double pf_irrv_, int at_param_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& fluid_viscosity_,
        Parameter<double> const& fluid_density_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& biot_modulus_,
        Parameter<double> const& porosity_)
        : material{std::move(material_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          solid_density(solid_density_),
          specific_body_force(specific_body_force_),
          pf_irrv(pf_irrv_),
          at_param(at_param_),
          intrinsic_permeability(intrinsic_permeability_),
          fluid_viscosity(fluid_viscosity_),
          fluid_density(fluid_density_),
          biot_coefficient(biot_coefficient_),
          biot_modulus(biot_modulus_),
          porosity(porosity_)
    {
    }

    HydroMechanicalPhaseFieldProcessData(
        HydroMechanicalPhaseFieldProcessData&& other)
        : material{std::move(other.material)},
          residual_stiffness(other.residual_stiffness),
          crack_resistance(other.crack_resistance),
          crack_length_scale(other.crack_length_scale),
          solid_density(other.solid_density),
          specific_body_force(other.specific_body_force),
          pf_irrv(other.pf_irrv),
          at_param(other.at_param),
          intrinsic_permeability(other.intrinsic_permeability),
          fluid_viscosity(other.fluid_viscosity),
          fluid_density(other.fluid_density),
          biot_coefficient(other.biot_coefficient),
          biot_modulus(other.biot_modulus),
          porosity(other.porosity),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    HydroMechanicalPhaseFieldProcessData(
        HydroMechanicalPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicalPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicalPhaseFieldProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    Parameter<double> const& residual_stiffness;
    Parameter<double> const& crack_resistance;
    Parameter<double> const& crack_length_scale;
    Parameter<double> const& solid_density;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& fluid_viscosity;
    Parameter<double> const& fluid_density;
    Parameter<double> const& biot_coefficient;
    Parameter<double> const& biot_modulus;
    Parameter<double> const& porosity;
    double poroelastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
    double pf_irrv = 0.05;
    int at_param = 2;
    double dt;
    double t;
};

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib
