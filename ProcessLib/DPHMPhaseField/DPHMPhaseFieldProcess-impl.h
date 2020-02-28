/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "DPHMPhaseFieldFEM.h"
#include "DPHMPhaseFieldProcess.h"
#include "DPHMPhaseFieldProcessData.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace DPHMPhaseField
{
template <int DisplacementDim>
DPHMPhaseFieldProcess<DisplacementDim>::DPHMPhaseFieldProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    DPHMPhaseFieldProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    int const mechanics_related_process_id, int const phase_field_process_id,
    int const frac_hydro_process_id, int const pore_hydro_process_id)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              false),
      _process_data(std::move(process_data)),
      _mechanics_related_process_id(mechanics_related_process_id),
      _phase_field_process_id(phase_field_process_id),
      _frac_hydro_process_id(frac_hydro_process_id),
      _pore_hydro_process_id(pore_hydro_process_id)
{
}

template <int DisplacementDim>
bool DPHMPhaseFieldProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
DPHMPhaseFieldProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and phase field process or hydro.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
DPHMPhaseFieldProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or hydro.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap&
DPHMPhaseFieldProcess<DisplacementDim>::getDOFTableByProcessID(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or hydro.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::constructDofTable()
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
void DPHMPhaseFieldProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, DPHMPhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        _mechanics_related_process_id, _phase_field_process_id,
        _frac_hydro_process_id, _pore_hydro_process_id);

    _secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &DPHMPhaseFieldLocalAssemblerInterface::getIntPtSigma));

    _secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &DPHMPhaseFieldLocalAssemblerInterface::getIntPtEpsilon));

    _process_data.ele_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "damage", MeshLib::MeshItemType::Cell,
        1);

    _process_data.ele_u_dot_grad_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "u_dot_grad_d",
        MeshLib::MeshItemType::Cell, 1);

    _process_data.width = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "width", MeshLib::MeshItemType::Cell,
        1);

    _process_data.width_prev = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "width_prev",
        MeshLib::MeshItemType::Cell, 1);

    _process_data.cum_grad_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "cum_grad_d",
        MeshLib::MeshItemType::Cell, 1);

    _process_data.ele_grad_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "grad_damage",
        MeshLib::MeshItemType::Cell, DisplacementDim);
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions()
{
    // Staggered scheme:
    // for the equations of deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_mechanics_related_process_id),
        _mechanics_related_process_id);
    // for the phase field
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_phase_field_process_id),
        _phase_field_process_id);
    // for fracture flow
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_frac_hydro_process_id), _frac_hydro_process_id);

    // for porous flow
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_pore_hydro_process_id), _pore_hydro_process_id);
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, GlobalVector const& x,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for DPHMPhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, dt, x, process_id, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, GlobalVector const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the staggered scheme
    if (process_id == _mechanics_related_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "poroelastic-deformation in "
            "DPHMPhaseFieldProcess for "
            "the staggered scheme.");
    }
    else if (process_id == _phase_field_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "phase field in "
            "DPHMPhaseFieldProcess for "
            "the staggered scheme.");
    }
    else if (process_id == _frac_hydro_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "fracture fluid flow in "
            "DPHMPhaseFieldProcess for "
            "the staggered scheme.");
    }
    else
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "porous fluid flow in "
            "DPHMPhaseFieldProcess for "
            "the staggered scheme.");
    }
    setCoupledSolutionsOfPreviousTimeStep();

    dof_tables.emplace_back(
        getDOFTableByProcessID(_mechanics_related_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_phase_field_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_frac_hydro_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_pore_hydro_process_id));

    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, dt, x, xdot, dxdot_dx, dx_dx,
        process_id, M, K, b, Jac, _coupled_solutions);
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    double const t,
    double const dt,
    const int process_id)
{
    DBUG("PreTimestep DPHMPhaseFieldProcess.");

    _process_data.dt = dt;
    _process_data.t = t;
    _x_previous_timestep[process_id] =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(*x[process_id]);
    //    MathLib::LinAlg::setLocalAccessibleVector(*_x_previous_timestep);
    //    MathLib::LinAlg::copy(x,*_x_previous_timestep);

    if (process_id != _mechanics_related_process_id)
    {
        return;
    }
    GlobalExecutor::executeMemberOnDereferenced(
        &DPHMPhaseFieldLocalAssemblerInterface::preTimestep, _local_assemblers,
        getDOFTable(process_id), *x[process_id], t, dt);
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    double const t,
    double const dt,
    int const process_id)
{
    DBUG("PostTimestep DPHMPhaseFieldProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &DPHMPhaseFieldLocalAssemblerInterface::postTimestep, _local_assemblers,
        getDOFTable(process_id), *x[process_id], t, dt);
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x,
    const double t,
    double const dt,
    const int process_id)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    dof_tables.emplace_back(
        getDOFTableByProcessID(_mechanics_related_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_phase_field_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_frac_hydro_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_pore_hydro_process_id));

    const bool use_monolithic_scheme = false;
    if (process_id == _phase_field_process_id)
    {
        _process_data.width_comp_visited =
            std::vector(_mesh.getNumberOfElements(), false);
        INFO("Fracture width computation");
        GlobalExecutor::executeMemberOnDereferenced(
            &DPHMPhaseFieldLocalAssemblerInterface::computeFractureNormal,
            _local_assemblers, dof_tables, _coupled_solutions);
        GlobalExecutor::executeMemberOnDereferenced(
            &DPHMPhaseFieldLocalAssemblerInterface::computeFractureWidth,
            _local_assemblers, dof_tables, t, _coupled_solutions, _mesh);

        INFO("After pressure solution, pressure is copied for fixed stress");
        GlobalExecutor::executeMemberOnDereferenced(
            &DPHMPhaseFieldLocalAssemblerInterface::postNonLinearSolver,
            _local_assemblers, getDOFTable(process_id), x, t, dt,
            use_monolithic_scheme);
    }
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::updateConstraints(
    GlobalVector& lower, GlobalVector& upper)
{
    lower.setZero();
    MathLib::LinAlg::setLocalAccessibleVector(
        *_x_previous_timestep[_phase_field_process_id]);
    MathLib::LinAlg::copy(*_x_previous_timestep[_phase_field_process_id],
                          upper);

    GlobalIndexType x_size =
        _x_previous_timestep[_phase_field_process_id]->size();

    for (GlobalIndexType i = 0; i < x_size; i++)
        if ((*_x_previous_timestep[_phase_field_process_id])[i] >
            _process_data.pf_irrv)
            upper.set(i, 1.0);
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<DisplacementDim>::
    setCoupledSolutionsOfPreviousTimeStepPerProcess(const int process_id)
{
    MathLib::LinAlg::setLocalAccessibleVector(
        *_x_previous_timestep[process_id]);
    _coupled_solutions->coupled_xs_t0[process_id] =
        _x_previous_timestep[process_id].get();
}

template <int DisplacementDim>
void DPHMPhaseFieldProcess<
    DisplacementDim>::setCoupledSolutionsOfPreviousTimeStep()
{
    _coupled_solutions->coupled_xs_t0.resize(4);
    setCoupledSolutionsOfPreviousTimeStepPerProcess(
        _mechanics_related_process_id);
    setCoupledSolutionsOfPreviousTimeStepPerProcess(_phase_field_process_id);
    setCoupledSolutionsOfPreviousTimeStepPerProcess(_frac_hydro_process_id);
    setCoupledSolutionsOfPreviousTimeStepPerProcess(_pore_hydro_process_id);
}

}  // namespace DPHMPhaseField
}  // namespace ProcessLib
