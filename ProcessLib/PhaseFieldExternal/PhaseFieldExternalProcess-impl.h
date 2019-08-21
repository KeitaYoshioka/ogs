/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PhaseFieldExternalFEM.h"
#include "PhaseFieldExternalProcess.h"
#include "PhaseFieldExternalProcessData.h"

#include <cassert>

#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace PhaseFieldExternal
{
template <int DisplacementDim>
PhaseFieldExternalProcess<DisplacementDim>::PhaseFieldExternalProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    PhaseFieldExternalProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    int const mechanics_related_process_id,
    int const phase_field_process_id)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              false),
      _process_data(std::move(process_data)),
      _mechanics_related_process_id(mechanics_related_process_id),
      _phase_field_process_id(phase_field_process_id)
{
}

template <int DisplacementDim>
bool PhaseFieldExternalProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
PhaseFieldExternalProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and phase field process or heat conduction.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
PhaseFieldExternalProcess<DisplacementDim>::getDOFTable(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or heat conduction.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap&
PhaseFieldExternalProcess<DisplacementDim>::getDOFTableByProcessID(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or heat conduction.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map_single_component);

    // For displacement equation.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    getProcessVariables(_mechanics_related_process_id)[0]
                        .get()
                        .getNumberOfComponents(),
                    [&]() { return *_mesh_subset_all_nodes; });

    std::vector<int> const vec_n_components{DisplacementDim};
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    // For phase field equation or the heat conduction.
    _sparsity_pattern_with_single_component = NumLib::computeSparsityPattern(
        *_local_to_global_index_map_single_component, _mesh);
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, PhaseFieldExternalLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        _mechanics_related_process_id, _phase_field_process_id);

    _secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &PhaseFieldExternalLocalAssemblerInterface::getIntPtSigma));

    _secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &PhaseFieldExternalLocalAssemblerInterface::getIntPtEpsilon));
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<DisplacementDim>::initializeBoundaryConditions()
{
    // Staggered scheme:
    // for the equations of temperature-deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_mechanics_related_process_id),
        _mechanics_related_process_id);
    // for the phase field
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_phase_field_process_id),
        _phase_field_process_id);
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble the equations for PhaseFieldExternalProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the staggered scheme
    if (_coupled_solutions->process_id == _mechanics_related_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "temperature-deformation in "
            "PhaseFieldExternalProcess for "
            "the staggered scheme.");
    }

    if (_coupled_solutions->process_id == _phase_field_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of"
            "phase field in "
            "PhaseFieldExternalProcess for "
            "the staggered scheme.");
    }

    dof_tables.emplace_back(
        getDOFTableByProcessID(_mechanics_related_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_phase_field_process_id));

    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
        Jac, _coupled_solutions);
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x,
    double const t,
    double const dt,
    const int process_id)
{
    DBUG("PreTimestep PhaseFieldExternalProcess.");

    _process_data.dt = dt;
    _process_data.t = t;
    _x_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);

    if (process_id != _mechanics_related_process_id)
    {
        return;
    }
    GlobalExecutor::executeMemberOnDereferenced(
        &PhaseFieldExternalLocalAssemblerInterface::preTimestep,
        _local_assemblers, getDOFTable(process_id), x, t, dt);
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x,
    double const /*t*/,
    double const /*dt*/,
    int const process_id)
{
    DBUG("PostTimestep PhaseFieldExternalProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &PhaseFieldExternalLocalAssemblerInterface::postTimestep,
        _local_assemblers, getDOFTable(process_id), x);
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<
    DisplacementDim>::postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                                         const double t,
                                                         const int process_id)
{
    if (process_id != _mechanics_related_process_id)
    {
        return;
    }

    DBUG("PostNonLinearSolver PhaseFieldExternalProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    const bool use_monolithic_scheme = false;
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        getDOFTable(process_id), x, t, use_monolithic_scheme);
}

template <int DisplacementDim>
void PhaseFieldExternalProcess<DisplacementDim>::updateConstraints(
    GlobalVector& lower, GlobalVector& upper)
{
    lower.setZero();
    MathLib::LinAlg::setLocalAccessibleVector(*_x_previous_timestep);
    MathLib::LinAlg::copy(*_x_previous_timestep, upper);

    GlobalIndexType x_begin = _x_previous_timestep->getRangeBegin();
    GlobalIndexType x_end = _x_previous_timestep->getRangeEnd();

    for (GlobalIndexType i = x_begin; i < x_end; i++)
        if ((*_x_previous_timestep)[i] > _process_data.pf_irrv)
            upper.set(i, 1.0);
}

}  // namespace PhaseFieldExternal
}  // namespace ProcessLib