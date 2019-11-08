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

#include "MaterialLib/Fluid/FluidType/FluidType.h"
#include "ParameterLib/Parameter.h"

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
namespace HydroMechanicalPhaseField
{
template <int DisplacementDim>
struct HydroMechanicalPhaseFieldProcessData
{
    HydroMechanicalPhaseFieldProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<
                          DisplacementDim>>>&& solid_materials_,
        ParameterLib::Parameter<double> const& residual_stiffness_,
        ParameterLib::Parameter<double> const& crack_resistance_,
        ParameterLib::Parameter<double> const& crack_length_scale_,
        ParameterLib::Parameter<double> const& solid_density_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_,
        int split_method_, double const reg_param_, double const pf_irrv_,
        double const li_disc_, double const cum_grad_d_CutOff_,
        int const at_param_,
        ParameterLib::Parameter<double> const& intrinsic_permeability_,
        ParameterLib::Parameter<double> const& fluid_viscosity_,
        ParameterLib::Parameter<double> const& fluid_density_,
        ParameterLib::Parameter<double> const& biot_coefficient_,
        ParameterLib::Parameter<double> const& biot_modulus_,
        ParameterLib::Parameter<double> const& drained_modulus_,
        ParameterLib::Parameter<double> const& porosity_,
        FluidType::Fluid_Type const fluid_type_,
        double const fluid_compressibility_,
        double const specific_gas_constant_,
        double const reference_temperature_,
        Eigen::Vector3d const& source_location_, double const source_)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          solid_density(solid_density_),
          specific_body_force(specific_body_force_),
          split_method(split_method_),
          reg_param(reg_param_),
          pf_irrv(pf_irrv_),
          li_disc(li_disc_),
          cum_grad_d_CutOff(cum_grad_d_CutOff_),
          at_param(at_param_),
          intrinsic_permeability(intrinsic_permeability_),
          fluid_viscosity(fluid_viscosity_),
          fluid_density(fluid_density_),
          biot_coefficient(biot_coefficient_),
          biot_modulus(biot_modulus_),
          drained_modulus(drained_modulus_),
          porosity(porosity_),
          fluid_type(fluid_type_),
          fluid_compressibility(fluid_compressibility_),
          specific_gas_constant(specific_gas_constant_),
          reference_temperature(reference_temperature_),
          source_location(source_location_),
          source(source_)
    {
    }

    HydroMechanicalPhaseFieldProcessData(
        HydroMechanicalPhaseFieldProcessData&& other) = default;

    //! Copies are forbidden.
    HydroMechanicalPhaseFieldProcessData(
        HydroMechanicalPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicalPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicalPhaseFieldProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;

    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& residual_stiffness;
    ParameterLib::Parameter<double> const& crack_resistance;
    ParameterLib::Parameter<double> const& crack_length_scale;
    ParameterLib::Parameter<double> const& solid_density;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    int split_method = 0;
    double reg_param = 0.01;
    double const pf_irrv = 0.05;
    double const li_disc = 60;
    double cum_grad_d_CutOff = 0.5;
    int const at_param = 2;
    ParameterLib::Parameter<double> const& intrinsic_permeability;
    ParameterLib::Parameter<double> const& fluid_viscosity;
    ParameterLib::Parameter<double> const& fluid_density;
    ParameterLib::Parameter<double> const& biot_coefficient;
    ParameterLib::Parameter<double> const& biot_modulus;
    ParameterLib::Parameter<double> const& drained_modulus;
    ParameterLib::Parameter<double> const& porosity;
    MeshLib::PropertyVector<double>* ele_grad_d = nullptr;
    MeshLib::PropertyVector<double>* ele_d = nullptr;
    MeshLib::PropertyVector<double>* ele_u_dot_grad_d = nullptr;
    MeshLib::PropertyVector<double>* width = nullptr;
    MeshLib::PropertyVector<double>* width_prev = nullptr;
    MeshLib::PropertyVector<double>* cum_grad_d = nullptr;
    FluidType::Fluid_Type const fluid_type;
    double const fluid_compressibility =
        std::numeric_limits<double>::quiet_NaN();
    double const specific_gas_constant =
        std::numeric_limits<double>::quiet_NaN();
    double const reference_temperature =
        std::numeric_limits<double>::quiet_NaN();
    Eigen::Vector3d const source_location;
    double const source = 0.0;
    double poroelastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
    double dt;
    double t;

    /// will be removed after linking with MPL
    double getFluidDensity(double const& t,
                           ParameterLib::SpatialPosition const& x_position,
                           double const& p_fr)
    {
        if (fluid_type == FluidType::Fluid_Type::INCOMPRESSIBLE_FLUID ||
            fluid_type == FluidType::Fluid_Type::COMPRESSIBLE_FLUID)
        {
            return fluid_density(t, x_position)[0];
        }
        if (fluid_type == FluidType::Fluid_Type::IDEAL_GAS)
        {
            return p_fr / (specific_gas_constant * reference_temperature);
        }
        OGS_FATAL("unknown fluid type %d", static_cast<int>(fluid_type));
    }

    /// will be removed after linking with MPL
    double getFluidCompressibility(double const& p_fr)
    {
        if (fluid_type == FluidType::Fluid_Type::INCOMPRESSIBLE_FLUID)
        {
            return 0.0;
        }
        if (fluid_type == FluidType::Fluid_Type::COMPRESSIBLE_FLUID)
        {
            return fluid_compressibility;
        }
        if (fluid_type == FluidType::Fluid_Type::IDEAL_GAS)
        {
            return 1.0 / p_fr;
        }
        OGS_FATAL("unknown fluid type %d", static_cast<int>(fluid_type));
    }
};

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib
