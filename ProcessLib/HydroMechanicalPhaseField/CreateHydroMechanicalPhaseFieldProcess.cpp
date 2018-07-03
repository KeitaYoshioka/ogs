/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHydroMechanicalPhaseFieldProcess.h"

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropicPhaseField.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"

#include "HydroMechanicalPhaseFieldProcess.h"
#include "HydroMechanicalPhaseFieldProcessData.h"

namespace ProcessLib
{
namespace HydroMechanicalPhaseField
{
template <int DisplacementDim>
std::unique_ptr<Process> createHydroMechanicalPhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HYDRO_MECHANICAL_PHASE_FIELD");
    DBUG("Create HydroMechanicalPhaseFieldProcess.");

    INFO(
        "Solve the coupling with the staggered scheme,"
        "which is the only option for TM-Phasefield in the current code");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    int hydro_process_id = 0;
    int mechanics_related_process_id = 1;
    int phase_field_process_id = 2;

    auto process_variable_P = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__process_variables__pressure}
         "pressure"});
    process_variables.push_back(std::move(process_variable_P));
    ProcessVariable* variable_p =
        &process_variables[process_variables.size() - 1][0].get();

    auto process_variable_u = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__process_variables__displacement}
         "displacement"});
    process_variables.push_back(std::move(process_variable_u));
    ProcessVariable* variable_u =
        &process_variables[process_variables.size() - 1][0].get();
    auto process_variable_ph = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__process_variables__phasefield}
         "phasefield"});
    process_variables.push_back(std::move(process_variable_ph));
    ProcessVariable* variable_ph =
        &process_variables[process_variables.size() - 1][0].get();

    DBUG("Associate displacement with process variable \'%s\'.",
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

    DBUG("Associate phase field with process variable \'%s\'.",
         variable_ph->getName().c_str());
    if (variable_ph->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phase field process variable '%s' is not a scalar variable but "
            "has "
            "%d components.",
            variable_ph->getName().c_str(),
            variable_ph->getNumberOfComponents());
    }

    DBUG("Associate temperature with process variable \'%s\'.",
         variable_p->getName().c_str());
    if (variable_p->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Temperature process variable '%s' is not a scalar variable but "
            "has "
            "%d components.",
            variable_p->getName().c_str(),
            variable_p->getNumberOfComponents());
    }

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const phasefield_parameters_config =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__phasefield_parameters}
        config.getConfigSubtree("phasefield_parameters");

    auto const type =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::PhaseFieldExtension<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropic")
    {
        auto elastic_model =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
                parameters, constitutive_relation_config);
        material = std::make_unique<
            MaterialLib::Solids::LinearElasticIsotropicPhaseField<
                DisplacementDim>>(
            std::move(elastic_model->getMaterialProperties()));
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Residual stiffness
    auto& residual_stiffness = findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__phasefield_parameters__residual_stiffness}
        "residual_stiffness", parameters, 1);
    DBUG("Use \'%s\' as residual stiffness.", residual_stiffness.name.c_str());

    // Crack resistance
    auto& crack_resistance = findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__phasefield_parameters__crack_resistance}
        "crack_resistance", parameters, 1);
    DBUG("Use \'%s\' as crack resistance.", crack_resistance.name.c_str());

    // Crack length scale
    auto& crack_length_scale = findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__phasefield_parameters__crack_length_scale}
        "crack_length_scale", parameters, 1);
    DBUG("Use \'%s\' as crack length scale.", crack_length_scale.name.c_str());


    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__reference_solid_density}
        "solid_density", parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.", solid_density.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICAL_PHASE_FIELD__specific_body_force}
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
    auto& intrinsic_permeability = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1);

    DBUG("Use \'%s\' as intrinsic permeability parameter.",
         intrinsic_permeability.name.c_str());

    // Fluid viscosity
    auto& fluid_viscosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__fluid_viscosity}
        "fluid_viscosity", parameters, 1);
    DBUG("Use \'%s\' as fluid viscosity parameter.",
         fluid_viscosity.name.c_str());

    // Fluid density
    auto& fluid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__fluid_density}
        "fluid_density", parameters, 1);
    DBUG("Use \'%s\' as fluid density parameter.", fluid_density.name.c_str());

    // Biot coefficient
    auto& biot_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__biot_coefficient}
        "biot_coefficient", parameters, 1);
    DBUG("Use \'%s\' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Biot's modulus
    auto& biot_modulus = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__biot_modulus}
        "biot_modulus", parameters, 1);

    // Porosity
    auto& porosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__porosity}
        "porosity", parameters, 1);
    DBUG("Use \'%s\' as porosity parameter.", porosity.name.c_str());


    auto pf_irrv_read =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD__pf_irrv}
        config.getConfigParameterOptional<double>("pf_irrv");

    double pf_irrv;
    if(pf_irrv_read)
        pf_irrv = *pf_irrv_read;
    else
        pf_irrv = 0.05;

    auto at_num =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD__at_num}
        config.getConfigParameterOptional<int>("at_num");

    int at_param;
    if(at_num && (*at_num == 1))
        at_param = 1;
    else
        at_param = 2;

    HydroMechanicalPhaseFieldProcessData<DisplacementDim> process_data{
        std::move(material),
        residual_stiffness,
        crack_resistance,
        crack_length_scale,
        solid_density,
        specific_body_force,
        pf_irrv,
        at_param,
        intrinsic_permeability,
        fluid_viscosity,
        fluid_density,
        biot_coefficient,
        biot_modulus,
        porosity};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"pressure_phasefield_displacement"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<HydroMechanicalPhaseFieldProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        mechanics_related_process_id, phase_field_process_id,
        hydro_process_id);
}

template std::unique_ptr<Process> createHydroMechanicalPhaseFieldProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createHydroMechanicalPhaseFieldProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib