/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateDPHMPhaseFieldProcess.h"

#include "DPHMPhaseFieldProcess.h"
#include "DPHMPhaseFieldProcessData.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace DPHMPhaseField
{
template <int DisplacementDim>
std::unique_ptr<Process> createDPHMPhaseFieldProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "DPHM_PHASE_FIELD");
    DBUG("Create DPHMPhaseFieldProcess.");

    INFO(
        "Solve the coupling with the staggered scheme,"
        "which is the only option for DPHM Phasefield in the "
        "current code");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    int frac_hydro_process_id = 0;
    int pore_hydro_process_id = 1;
    int mechanics_related_process_id = 2;
    int phase_field_process_id = 3;

    auto process_variable_p_f = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__process_variables__frac_pressure}
         "frac_pressure"});
    process_variables.push_back(std::move(process_variable_p_f));
    ProcessVariable* variable_p_f =
        &process_variables[process_variables.size() - 1][0].get();

    auto process_variable_p_p = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__process_variables__pore_pressure}
         "pore_pressure"});
    process_variables.push_back(std::move(process_variable_p_p));
    ProcessVariable* variable_p_p =
        &process_variables[process_variables.size() - 1][0].get();

    auto process_variable_u = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__process_variables__displacement}
         "displacement"});
    process_variables.push_back(std::move(process_variable_u));
    ProcessVariable* variable_u =
        &process_variables[process_variables.size() - 1][0].get();

    auto process_variable_ph = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__process_variables__phasefield}
         "phasefield"});
    process_variables.push_back(std::move(process_variable_ph));
    ProcessVariable* variable_ph =
        &process_variables[process_variables.size() - 1][0].get();

    DBUG("Associate fracture pressure with process variable '%s'.",
         variable_p_f->getName().c_str());
    if (variable_p_f->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Fracture pressure process variable '%s' is not a scalar variable "
            "but "
            "has "
            "%d components.",
            variable_p_f->getName().c_str(),
            variable_p_f->getNumberOfComponents());
    }

    DBUG("Associate pore pressure with process variable '%s'.",
         variable_p_p->getName().c_str());
    if (variable_p_p->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pore pressure process variable '%s' is not a scalar variable but "
            "has "
            "%d components.",
            variable_p_p->getName().c_str(),
            variable_p_p->getNumberOfComponents());
    }

    DBUG("Associate displacement with process variable '%s'.",
         variable_u->getName().c_str());

    if (variable_u->getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            variable_u->getName().c_str(),
            variable_u->getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate phase field with process variable '%s'.",
         variable_ph->getName().c_str());
    if (variable_ph->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phasefield process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_ph->getName().c_str(),
            variable_ph->getNumberOfComponents());
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    auto const phasefield_parameters_config =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__phasefield_parameters}
        config.getConfigSubtree("phasefield_parameters");

    // Residual stiffness
    auto& residual_stiffness = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__phasefield_parameters__residual_stiffness}
        "residual_stiffness", parameters, 1);
    DBUG("Use '%s' as residual stiffness.", residual_stiffness.name.c_str());

    // Crack resistance
    auto& crack_resistance = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__phasefield_parameters__crack_resistance}
        "crack_resistance", parameters, 1);
    DBUG("Use '%s' as crack resistance.", crack_resistance.name.c_str());

    // Crack length scale
    auto& crack_length_scale = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__phasefield_parameters__crack_length_scale}
        "crack_length_scale", parameters, 1);
    DBUG("Use '%s' as crack length scale.", crack_length_scale.name.c_str());

    // Solid density
    auto& solid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__reference_solid_density}
        "solid_density", parameters, 1);
    DBUG("Use '%s' as solid density parameter.", solid_density.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (specific_body_force.size() != DisplacementDim)
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                specific_body_force.size(), DisplacementDim);

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    // Intrinsic permeability
    auto& intrinsic_permeability = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1);

    DBUG("Use '%s' as intrinsic permeability parameter.",
         intrinsic_permeability.name.c_str());

    // Fluid viscosity
    auto& fluid_viscosity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__fluid_viscosity}
        "fluid_viscosity", parameters, 1);
    DBUG("Use '%s' as fluid viscosity parameter.",
         fluid_viscosity.name.c_str());

    // Fluid density
    auto& fluid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__fluid_density}
        "fluid_density", parameters, 1);
    DBUG("Use '%s' as fluid density parameter.", fluid_density.name.c_str());

    // Grain modulus
    auto& grain_modulus = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__grain_modulus}
        "grain_modulus", parameters, 1);

    // drained modulus
    auto& drained_modulus = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__drained_modulus}
        "drained_modulus", parameters, 1);

    // Porosity
    auto& porosity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__DPHM_PHASE_FIELD__porosity}
        "porosity", parameters, 1);
    DBUG("Use '%s' as porosity parameter.", porosity.name.c_str());

    auto const fluid_type = FluidType::strToFluidType(
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__fluid_type}
        config.getConfigParameter<std::string>("fluid_type"));
    DBUG("Use 'fluid_type' as fluid type parameter.");
    // Reference temperature
    double const reference_temperature =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__reference_temperature}
        config.getConfigParameter<double>(
            "reference_temperature", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'reference_temperature' as reference temperature.");

    // Specific gas constant
    double const specific_gas_constant =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__specific_gas_constant}
        config.getConfigParameter<double>(
            "specific_gas_constant", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'specific_gas_constant' as specific gas constant.");

    // Fluid compressibility
    double const fluid_compressibility =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__fluid_compressibility}
        config.getConfigParameter<double>(
            "fluid_compressibility", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'fluid_compressibility' as fluid compressibility parameter.");

    Eigen::Vector3d const source_location = [&]() {
        auto const v =
            //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__source_location}
            config.getConfigParameter<std::vector<double>>("source_location");
        if (v.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the source location vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                v.size(), DisplacementDim);
        }
        Eigen::Vector3d vec3 = Eigen::Vector3d::Zero();
        std::copy_n(v.data(), DisplacementDim, vec3.data());
        return vec3;
    }();

    auto source_read =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__source}
        config.getConfigParameterOptional<double>("source");

    double source;
    if (source_read)
        source = *source_read;
    else
        source = 0.0;

    auto reg_param_read =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__reg_param}
        config.getConfigParameterOptional<double>("reg_param");

    double reg_param;
    if (reg_param_read)
        reg_param = *reg_param_read;
    else
        reg_param = 0.01;

    auto pf_irrv_read =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__pf_irrv}
        config.getConfigParameterOptional<double>("pf_irrv");

    double pf_irrv;
    if (pf_irrv_read)
        pf_irrv = *pf_irrv_read;
    else
        pf_irrv = 0.05;

    auto li_disc_read =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__li_disc}
        config.getConfigParameterOptional<double>("li_disc");

    double li_disc;
    if (li_disc_read)
        li_disc = *li_disc_read;
    else
        li_disc = 60;

    auto cut_off_read =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__cut_off}
        config.getConfigParameterOptional<double>("cut_off");

    double cut_off;
    if (cut_off_read)
        cut_off = *cut_off_read;
    else
        cut_off = 60;

    auto at_num =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__at_num}
        config.getConfigParameterOptional<int>("at_num");

    int at_param;
    if (at_num && (*at_num == 1))
        at_param = 1;
    else
        at_param = 2;

    auto split =
        //! \ogs_file_param{prj__processes__process__DPHM_PHASE_FIELD__split_method}
        config.getConfigParameterOptional<int>("split_method");

    int split_method;
    if (split && (*split == 1))
        split_method = 1;
    else if (split && (*split == 2))
        split_method = 2;
    else if (split && (*split == 3))
        split_method = 3;
    else
        split_method = 0;

    auto const leakoff_scheme =
        //! \ogs_file_param{prj__processes__process__DPHMPHASE_FIELD__leak_off_scheme}
        config.getConfigParameterOptional<std::string>("leakoff_scheme");
    if (leakoff_scheme &&
        ((*leakoff_scheme != "yes") && (*leakoff_scheme != "no")))
    {
        OGS_FATAL(
            "leakoff_scheme must be 'yes' or 'no' but "
            "'%s' was given",
            leakoff_scheme->c_str());
    }

    const bool leak_off =
        (leakoff_scheme && (*leakoff_scheme == "yes"));

    DPHMPhaseFieldProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(solid_constitutive_relations),
        residual_stiffness,
        crack_resistance,
        crack_length_scale,
        solid_density,
        specific_body_force,
        split_method,
        reg_param,
        pf_irrv,
        li_disc,
        cut_off,
        at_param,
        intrinsic_permeability,
        fluid_viscosity,
        fluid_density,
        grain_modulus,
        drained_modulus,
        porosity,
        fluid_type,
        fluid_compressibility,
        reference_temperature,
        specific_gas_constant,
        source_location,
        source};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"pressure_displacement_phasefield"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<DPHMPhaseFieldProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(named_function_caller), mechanics_related_process_id,
        phase_field_process_id, frac_hydro_process_id, pore_hydro_process_id);
}

template std::unique_ptr<Process> createDPHMPhaseFieldProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createDPHMPhaseFieldProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace DPHMPhaseField
}  // namespace ProcessLib