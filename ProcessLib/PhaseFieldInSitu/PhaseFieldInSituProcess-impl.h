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

#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "PhaseFieldInSituFEM.h"
#include "PhaseFieldInSituProcess.h"
#include "PhaseFieldInSituProcessData.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace PhaseFieldInSitu
{
template <int DisplacementDim>
PhaseFieldInSituProcess<DisplacementDim>::PhaseFieldInSituProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    PhaseFieldInSituProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    int const mechanics_process0_id,
    int const mechanics_process1_id,
    int const phase_field_process_id)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              false),
      _process_data(std::move(process_data)),
      _mechanics_process0_id(mechanics_process0_id),
      _mechanics_process1_id(mechanics_process1_id),
      _phase_field_process_id(phase_field_process_id)
{
}

template <int DisplacementDim>
bool PhaseFieldInSituProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
PhaseFieldInSituProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    if (process_id == _mechanics_process0_id ||
        process_id == _mechanics_process1_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and phase field.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
PhaseFieldInSituProcess<DisplacementDim>::getDOFTable(
    const int process_id) const
{
    if (process_id == _mechanics_process0_id ||
        process_id == _mechanics_process1_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or hydro.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap&
PhaseFieldInSituProcess<DisplacementDim>::getDOFTableByProcessID(
    const int process_id) const
{
    if (process_id == _mechanics_process0_id ||
        process_id == _mechanics_process1_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or hydro.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
void PhaseFieldInSituProcess<DisplacementDim>::constructDofTable()
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

    // For displacement0 equation.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    getProcessVariables(_mechanics_process0_id)[0]
                        .get()
                        .getNumberOfComponents(),
                    [&]() { return *_mesh_subset_all_nodes; });

    std::vector<int> const vec_n_components{DisplacementDim};
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    // For displacement1 equation.
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    getProcessVariables(_mechanics_process1_id)[0]
                        .get()
                        .getNumberOfComponents(),
                    [&]() { return *_mesh_subset_all_nodes; });

    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    // For phase field equation
    _sparsity_pattern_with_single_component = NumLib::computeSparsityPattern(
        *_local_to_global_index_map_single_component, _mesh);
}

template <int DisplacementDim>
void PhaseFieldInSituProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, PhaseFieldInSituLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        _mechanics_process0_id, _mechanics_process1_id,
        _phase_field_process_id);

    _secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &PhaseFieldInSituLocalAssemblerInterface::getIntPtSigma));

    _secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &PhaseFieldInSituLocalAssemblerInterface::getIntPtEpsilon));
}

template <int DisplacementDim>
void PhaseFieldInSituProcess<DisplacementDim>::initializeBoundaryConditions()
{
    // Staggered scheme:
    // for displacement0.
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_mechanics_process0_id), _mechanics_process0_id);

    // for displacement1
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_mechanics_process1_id), _mechanics_process1_id);

    // for phase-field
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_phase_field_process_id),
        _phase_field_process_id);
}

template <int DisplacementDim>
void PhaseFieldInSituProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, GlobalVector const& x,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble the equations for PhaseFieldInSituProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, dt, x, process_id, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void PhaseFieldInSituProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, GlobalVector const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the staggered scheme
    if (process_id == _mechanics_process0_id)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "deformation0 in "
            "PhaseFieldInSituProcess for "
            "the staggered scheme.");
    }
    else if (process_id == _mechanics_process1_id)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "deformation1 in "
            "PhaseFieldInSituProcess for "
            "the staggered scheme.");
    }
    else
    {
        DBUG(
            "Assemble the Jacobian equations of z"
            "phase field in "
            "PhaseFieldInSituProcess for "
            "the staggered scheme.");
    }
    dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process0_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process1_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_phase_field_process_id));

    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, dt, x, xdot, dxdot_dx, dx_dx,
        process_id, M, K, b, Jac, _coupled_solutions);
}

template <int DisplacementDim>
void PhaseFieldInSituProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep PhaseFieldInSituProcess.");

    _process_data.dt = dt;
    _process_data.t = t;
    _process_data.injected_volume = _process_data.t;
    _x_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(*x[process_id]);

    GlobalExecutor::executeMemberOnDereferenced(
        &PhaseFieldInSituLocalAssemblerInterface::preTimestep,
        _local_assemblers, getDOFTable(process_id), *x[process_id], t, dt);
    /*    if (_coupled_solutions->process_id == _mechanics_process0_id)
        {
            std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
                dof_tables;

            dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process0_id));
            dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process1_id));
            dof_tables.emplace_back(
                getDOFTableByProcessID(_phase_field_process_id));

            auto& u_p =
                _coupled_solutions->coupled_xs[_mechanics_process0_id].get();
            auto& u_s =
                _coupled_solutions->coupled_xs[_mechanics_process1_id].get();
            // u_p = 1/p * u_p - 1/p * u_s
            MathLib::LinAlg::axpby(
                const_cast<GlobalVector&>(u_p), 1 / _process_data.pressure,
                -1 / _process_data.pressure, const_cast<GlobalVector&>(u_s));
            INFO("u_p is scaled back for a new time step");
        }*/
}

template <int DisplacementDim>
void PhaseFieldInSituProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    int const process_id)
{
    DBUG("PostTimestep PhaseFieldInSituProcess.");

    if (process_id == _phase_field_process_id)
    {
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeMemberOnDereferenced(
            &PhaseFieldInSituLocalAssemblerInterface::postTimestep,
            _local_assemblers, getDOFTable(process_id), *x[process_id], t, dt);

        _process_data.elastic_energy = 0.0;
        _process_data.surface_energy = 0.0;
        _process_data.pressure_work = 0.0;
        _process_data.crack_length = 0.0;

        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;
        dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process0_id));
        dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process1_id));
        dof_tables.emplace_back(
            getDOFTableByProcessID(_phase_field_process_id));
        /*
                auto& u_p =
                    _coupled_solutions->coupled_xs[_mechanics_process0_id].get();
                auto& u_s =
                    _coupled_solutions->coupled_xs[_mechanics_process1_id].get();

                // u_p holds the unscaled displacement
                // u_p = p*u_p + u_s
                MathLib::LinAlg::aypx(const_cast<GlobalVector&>(u_p),
                                      _process_data.pressure,
                                      const_cast<GlobalVector&>(u_s));
                INFO("u_p is superposed for output ");
                */
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &PhaseFieldInSituLocalAssemblerInterface::computeEnergy,
            _local_assemblers, pv.getActiveElementIDs(), dof_tables, *x[process_id], t, dt,
            _process_data.elastic_energy, _process_data.surface_energy,
            _process_data.crack_length, _process_data.pressure_work,
            _coupled_solutions);

#ifdef USE_PETSC
        double const elastic_energy = _process_data.elastic_energy;
        MPI_Allreduce(&elastic_energy, &_process_data.elastic_energy, 1,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        double const surface_energy = _process_data.surface_energy;
        MPI_Allreduce(&surface_energy, &_process_data.surface_energy, 1,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        double const crack_length = _process_data.crack_length;
        MPI_Allreduce(&crack_length, &_process_data.crack_length, 1, MPI_DOUBLE,
                      MPI_SUM, PETSC_COMM_WORLD);
        double const pressure_work = _process_data.pressure_work;
        MPI_Allreduce(&pressure_work, &_process_data.pressure_work, 1,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
#endif

        INFO("Elastic energy: %g Surface energy: %g Pressure work: %g ",
             _process_data.elastic_energy, _process_data.crack_length,
             _process_data.surface_energy);
    }
}

template <int DisplacementDim>
void PhaseFieldInSituProcess<
    DisplacementDim>::postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                                         const double t,
                                                         double const /*dt*/,
                                                         const int process_id)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    if (process_id == _mechanics_process0_id)
    {
        _process_data.crack_volume0 = 0.0;
        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;
        dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process0_id));
        dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process1_id));
        dof_tables.emplace_back(
            getDOFTableByProcessID(_phase_field_process_id));
        DBUG(
            "PostNonLinearSolver crack volume computation for "
            "mechanics_process0.");

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &PhaseFieldInSituLocalAssemblerInterface::computeCrackIntegral,
            _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, t,
            _process_data.crack_volume0, _coupled_solutions);
#ifdef USE_PETSC
        double const crack_volume = _process_data.crack_volume0;
        MPI_Allreduce(&crack_volume, &_process_data.crack_volume0, 1,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
#endif
        INFO("Integral of crack with u_p: %g", _process_data.crack_volume0);
    }
    else if (process_id == _mechanics_process1_id)
    {
        _process_data.crack_volume1 = 0.0;
        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;

        dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process0_id));
        dof_tables.emplace_back(getDOFTableByProcessID(_mechanics_process1_id));
        dof_tables.emplace_back(
            getDOFTableByProcessID(_phase_field_process_id));

        DBUG(
            "PostNonLinearSolver crack volume computation for "
            "mechanics_process0.");

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &PhaseFieldInSituLocalAssemblerInterface::computeCrackIntegral,
            _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, t,
            _process_data.crack_volume1, _coupled_solutions);
#ifdef USE_PETSC
        double const crack_volume = _process_data.crack_volume1;
        MPI_Allreduce(&crack_volume, &_process_data.crack_volume1, 1,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
#endif
        INFO("Integral of crack with u_s: %g", _process_data.crack_volume1);

        if (_process_data.propagating_crack)
        {
            /*           _process_data.crack_volume0 =
                           _process_data.crack_volume0 /
               _process_data.unity_pressure; */
            _process_data.pressure_old = _process_data.pressure;
            // p = (V_inj - V_fs)/V_fp
            // V_fs: stress only, V_fp: pressure only
            _process_data.pressure =
                (_process_data.injected_volume - _process_data.crack_volume1) /
                _process_data.crack_volume0;
            _process_data.pressure_error =
                std::fabs(_process_data.pressure_old - _process_data.pressure) /
                _process_data.pressure;
            INFO("Internal pressure: %g and Pressure error: %.4e",
                 _process_data.pressure, _process_data.pressure_error);
            GlobalVector u_p{
                *_coupled_solutions->coupled_xs[_mechanics_process0_id]};
            GlobalVector u_s{
                *_coupled_solutions->coupled_xs[_mechanics_process1_id]};

            // u_p holds the unscaled displacement
            // u_p = p*u_p + u_s
            MathLib::LinAlg::aypx(u_p, _process_data.pressure, u_s);
            MathLib::LinAlg::copy(
                u_p,
                const_cast<GlobalVector&>(
                    *_coupled_solutions->coupled_xs[_mechanics_process0_id]));
        }
    }
    else
    {
        if (_process_data.propagating_crack)
        {
            GlobalVector u_p{
                *_coupled_solutions->coupled_xs[_mechanics_process0_id]};
            GlobalVector u_s{
                *_coupled_solutions->coupled_xs[_mechanics_process1_id]};

            // u_p = 1/p * u_p - 1/p * u_s
            MathLib::LinAlg::axpby(u_p, 1 / _process_data.pressure,
                                   -1 / _process_data.pressure, u_s);
            MathLib::LinAlg::copy(
                u_p,
                const_cast<GlobalVector&>(
                    *_coupled_solutions->coupled_xs[_mechanics_process0_id]));
            INFO("u_p is scaled back for non-linear iteration");
        }
    }

}  // namespace PhaseFieldInSitu

template <int DisplacementDim>
void PhaseFieldInSituProcess<DisplacementDim>::updateConstraints(
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

}  // namespace PhaseFieldInSitu
}  // namespace ProcessLib
