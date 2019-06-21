/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreatePhaseFieldInSituProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "PhaseFieldInSituProcess.h"
#include "PhaseFieldInSituProcessData.h"

namespace ProcessLib
{
namespace PhaseFieldInSitu
{
template <int DisplacementDim>
std::unique_ptr<Process> createPhaseFieldInSituProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "PHASE_FIELD_INSITU");
    DBUG("Create PhaseFieldInSituProcess.");

    INFO(
        "Solve the coupling with the staggered scheme,"
        "which is the only option for Phasefield in the current code");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__PHASE_FIELD_INSITU__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    // u_p
    int mechanics_process0_id = 0;
    // u_s
    int mechanics_process1_id = 1;
    int phase_field_process_id = 2;

    auto process_variable_u0 = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_INSITU__process_variables__displacement0}
         "displacement0"});
    process_variables.push_back(std::move(process_variable_u0));
    ProcessVariable* variable_u0 =
        &process_variables[process_variables.size() - 1][0].get();

    auto process_variable_u1 = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_INSITU__process_variables__displacement1}
         "displacement1"});
    process_variables.push_back(std::move(process_variable_u1));
    ProcessVariable* variable_u1 =
        &process_variables[process_variables.size() - 1][0].get();

    auto process_variable_ph = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_INSITU__process_variables__phasefield}
         "phasefield"});
    process_variables.push_back(std::move(process_variable_ph));
    ProcessVariable* variable_ph =
        &process_variables[process_variables.size() - 1][0].get();

    DBUG("Associate displacement0 with process variable \'%s\'.",
         variable_u0->getName().c_str());
    if (variable_u0->getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement0 dimension: got %d, expected %d",
            variable_u0->getName().c_str(),
            variable_u0->getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate displacement1 with process variable \'%s\'.",
         variable_u1->getName().c_str());
    if (variable_u1->getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement1 dimension: got %d, expected %d",
            variable_u1->getName().c_str(),
            variable_u1->getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate phase field with process variable \'%s\'.",
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
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD_INSITU__phasefield_parameters}
        config.getConfigSubtree("phasefield_parameters");

    // Residual stiffness
    auto& residual_stiffness = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_INSITU__phasefield_parameters__residual_stiffness}
        "residual_stiffness", parameters, 1);
    DBUG("Use \'%s\' as residual stiffness.", residual_stiffness.name.c_str());

    // Crack resistance
    auto& crack_resistance = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_INSITU__phasefield_parameters__crack_resistance}
        "crack_resistance", parameters, 1);
    DBUG("Use \'%s\' as crack resistance.", crack_resistance.name.c_str());

    // Crack length scale
    auto& crack_length_scale = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_INSITU__phasefield_parameters__crack_length_scale}
        "crack_length_scale", parameters, 1);
    DBUG("Use \'%s\' as crack length scale.", crack_length_scale.name.c_str());

    // Solid density
    auto& solid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_INSITU__reference_solid_density}
        "solid_density", parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.", solid_density.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__PHASE_FIELD_INSITU__specific_body_force}
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

    auto pf_irrv_read =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD_INSITU__pf_irrv}
        config.getConfigParameterOptional<double>("pf_irrv");

    auto const crack_scheme =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD_INSITU__hydro_crack_scheme}
        config.getConfigParameterOptional<std::string>("hydro_crack_scheme");
    if (crack_scheme &&
        ((*crack_scheme != "propagating") && (*crack_scheme != "static")))
    {
        OGS_FATAL(
            "hydro_crack_scheme must be \"propagating\" or \"static\" but "
            "\"%s\" was given",
            crack_scheme->c_str());
    }

    const bool propagating_crack =
        (crack_scheme && (*crack_scheme == "propagating"));
    const bool crack_pressure =
        (crack_scheme &&
         ((*crack_scheme == "propagating") || (*crack_scheme == "static")));

    double pf_irrv;
    if (pf_irrv_read)
        pf_irrv = *pf_irrv_read;
    else
        pf_irrv = 0.05;

    auto at_num =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD__at_num}
        config.getConfigParameterOptional<int>("at_num");

    int at_param;
    if (at_num && (*at_num == 1))
        at_param = 1;
    else if (at_num && (*at_num == 3))
        at_param = 3;
    else
        at_param = 2;
    auto split =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD__split_method}
        config.getConfigParameterOptional<int>("split_method");

    int split_method;
    if (split && (*split == 1))
        split_method = 1;
    else
        split_method = 0;

    PhaseFieldInSituProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(solid_constitutive_relations),
        residual_stiffness,
        crack_resistance,
        crack_length_scale,
        solid_density,
        specific_body_force,
        propagating_crack,
        split_method,
        crack_pressure,
        pf_irrv,
        at_param};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"displacement0_displacement1_phasefield"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<PhaseFieldInSituProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(named_function_caller), mechanics_process0_id,
        mechanics_process1_id, phase_field_process_id);
}

template std::unique_ptr<Process> createPhaseFieldInSituProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createPhaseFieldInSituProcess<3>(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config);

}  // namespace PhaseFieldInSitu
}  // namespace ProcessLib
