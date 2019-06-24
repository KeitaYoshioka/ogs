/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreatePhaseFieldExternalProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "PhaseFieldExternalProcess.h"
#include "PhaseFieldExternalProcessData.h"

namespace ProcessLib
{
namespace PhaseFieldExternal
{
template <int DisplacementDim>
std::unique_ptr<Process> createPhaseFieldExternalProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "PHASE_FIELD_EXTERNAL");
    DBUG("Create PhaseFieldExternalProcess.");

    INFO(
        "Solve the coupling with the staggered scheme,"
        "which is the only option for TM-Phasefield in the current code");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__PHASE_FIELD_EXTERNAL__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    int mechanics_related_process_id = 0;
    int phase_field_process_id = 1;

    auto process_variable_u = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__process_variables__displacement}
         "displacement"});
    process_variables.push_back(std::move(process_variable_u));
    ProcessVariable* variable_u =
        &process_variables[process_variables.size() - 1][0].get();
    auto process_variable_ph = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__process_variables__phasefield}
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
            "Phasefield process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_ph->getName().c_str(),
            variable_ph->getNumberOfComponents());
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    auto const phasefield_parameters_config =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD_EXTERNAL__phasefield_parameters}
        config.getConfigSubtree("phasefield_parameters");

    // Residual stiffness
    auto& residual_stiffness = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__phasefield_parameters__residual_stiffness}
        "residual_stiffness", parameters, 1);
    DBUG("Use \'%s\' as residual stiffness.", residual_stiffness.name.c_str());

    // Crack resistance
    auto& crack_resistance = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__phasefield_parameters__crack_resistance}
        "crack_resistance", parameters, 1);
    DBUG("Use \'%s\' as crack resistance.", crack_resistance.name.c_str());

    // Crack length scale
    auto& crack_length_scale = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__phasefield_parameters__crack_length_scale}
        "crack_length_scale", parameters, 1);
    DBUG("Use \'%s\' as crack length scale.", crack_length_scale.name.c_str());

    // Solid density
    auto& solid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__reference_solid_density}
        "solid_density", parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.", solid_density.name.c_str());

    // Linear thermal expansion coefficient
    auto& linear_thermal_expansion_coefficient = ParameterLib::findParameter<
        double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__reference__linear_thermal_expansion_coefficient}
        "linear_thermal_expansion_coefficient", parameters, 1);
    DBUG("Use \'%s\' as linear thermal expansion coefficient.",
         linear_thermal_expansion_coefficient.name.c_str());

    // Pressure provided by external tool
    auto& pressure_ext = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__reference__pressure_ext}
        "pressure_ext", parameters, 1);
    DBUG("Use \'%s\' as pressure provided by extneral tool.",
         pressure_ext.name.c_str());

    // Temperature provided by external tool
    auto& temperature_ext = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__reference__temperature_ext}
        "temperature_ext", parameters, 1);
    DBUG("Use \'%s\' as temperature provided by extneral tool.",
         temperature_ext.name.c_str());

    // Biot's coefficient
    auto& biot_coefficient = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD_EXTERNAL__biot_coefficient}
        "biot_coefficient", parameters, 1);
    DBUG("Use \'%s\' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Reference temperature
    const double reference_temperature =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD_EXTERNAL__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__PHASE_FIELD_EXTERNAL__specific_body_force}
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

    auto reg_param_read =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD__reg_param}
        config.getConfigParameterOptional<double>("reg_param");

    double reg_param;
    if (reg_param_read)
        reg_param = *reg_param_read;
    else
        reg_param = 0.01;

    auto pf_irrv_read =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD__pf_irrv}
        config.getConfigParameterOptional<double>("pf_irrv");

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

    PhaseFieldExternalProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(solid_constitutive_relations),
        residual_stiffness,
        crack_resistance,
        crack_length_scale,
        solid_density,
        linear_thermal_expansion_coefficient,
        pressure_ext,
        temperature_ext,
        biot_coefficient,
        reference_temperature,
        specific_body_force,
        split_method,
        reg_param,
        pf_irrv,
        at_param};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"temperature_phasefield_displacement"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<PhaseFieldExternalProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(named_function_caller), mechanics_related_process_id,
        phase_field_process_id);
}

template std::unique_ptr<Process> createPhaseFieldExternalProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createPhaseFieldExternalProcess<3>(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order, BaseLib::ConfigTree const& config);

}  // namespace PhaseFieldExternal
}  // namespace ProcessLib
