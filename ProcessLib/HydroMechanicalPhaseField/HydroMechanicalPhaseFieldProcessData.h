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
        int split_method_, int fperm_method_, int pf_scaling_method_, int fixed_strs_method_,
        double const reg_param_, double const pf_irrv_, double const li_disc_,
        double const cum_grad_d_CutOff_, int const at_param_,
        double const fixed_strs1_, double const fixed_strs2_,
        bool poroelastic_coupling_,
        ParameterLib::Parameter<double> const& intrinsic_permeability_,
        ParameterLib::Parameter<double> const& fluid_viscosity_,
        ParameterLib::Parameter<double> const& fluid_density_,
        ParameterLib::Parameter<double> const& grain_modulus_,
        ParameterLib::Parameter<double> const& drained_modulus_,
        ParameterLib::Parameter<double> const& porosity_,
        ParameterLib::Parameter<double> const& width_init_,
        double const geostatic_pressure_, double const fperm_param1_,
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
          fperm_method(fperm_method_),
          pf_scaling_method(pf_scaling_method_),
          fixed_strs_method(fixed_strs_method_),
          reg_param(reg_param_),
          pf_irrv(pf_irrv_),
          li_disc(li_disc_),
          cum_grad_d_CutOff(cum_grad_d_CutOff_),
          at_param(at_param_),
          fixed_strs1(fixed_strs1_),
          fixed_strs2(fixed_strs2_),
          poroelastic_coupling(poroelastic_coupling_),
          intrinsic_permeability(intrinsic_permeability_),
          fluid_viscosity(fluid_viscosity_),
          fluid_density(fluid_density_),
          grain_modulus(grain_modulus_),
          drained_modulus(drained_modulus_),
          porosity(porosity_),
          width_init(width_init_),
          geostatic_pressure(geostatic_pressure_),
          fperm_param1(fperm_param1_),
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
    int fperm_method = 0;
    int pf_scaling_method = 0;
    int fixed_strs_method = 0;
    double reg_param = 0.01;
    double const pf_irrv = 0.05;
    double const li_disc = 60;
    double cum_grad_d_CutOff = 0.5;
    int const at_param = 2;
    double const fixed_strs1 = 0.0;
    double const fixed_strs2 = 0.0;
    bool poroelastic_coupling = true;
    ParameterLib::Parameter<double> const& intrinsic_permeability;
    ParameterLib::Parameter<double> const& fluid_viscosity;
    ParameterLib::Parameter<double> const& fluid_density;
    ParameterLib::Parameter<double> const& grain_modulus;
    ParameterLib::Parameter<double> const& drained_modulus;
    ParameterLib::Parameter<double> const& porosity;
    ParameterLib::Parameter<double> const& width_init;
    MeshLib::PropertyVector<double>* ele_grad_d = nullptr;
    MeshLib::PropertyVector<double>* ele_d = nullptr;
    MeshLib::PropertyVector<double>* ele_u_dot_grad_d = nullptr;
    MeshLib::PropertyVector<double>* width = nullptr;
    MeshLib::PropertyVector<double>* width_nl_prev = nullptr;
    MeshLib::PropertyVector<double>* width_prev = nullptr;
    MeshLib::PropertyVector<double>* cum_grad_d = nullptr;
    MeshLib::PropertyVector<double>* frac_velocity = nullptr;
    double const geostatic_pressure = 0.0;
    double const fperm_param1 = 1.0;

    std::vector<bool> width_comp_visited;
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
