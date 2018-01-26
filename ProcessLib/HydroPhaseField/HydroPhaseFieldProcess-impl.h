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

#include "NumLib/DOF/ComputeSparsityPattern.h"

#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

#include "HydroPhaseFieldFEM.h"

namespace ProcessLib
{
namespace HydroPhaseField
{
template <int DisplacementDim>
HydroPhaseFieldProcess<DisplacementDim>::HydroPhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HydroPhaseFieldProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme,
    int const mechanics_related_process_id,
    int const hydro_process_id,
    int const phase_field_process_id)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _process_data(std::move(process_data),
                    _mechanics_related_process_id(mechanics_related_process_id),
                    _hydro_process_id(hydro_process_id),
                    _phase_field_process_id(phase_field_process_id))
{
}

template <int DisplacementDim>
bool HydroPhaseFieldProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
HydroPhaseFieldProcess<DisplacementDim>::getMatrixSpecifications(
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
HydroPhaseFieldProcess<DisplacementDim>::getDOFTable(const int process_id) const
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
HydroPhaseFieldProcess<DisplacementDim>::getDOFTableByProcessID(
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
void HydroPhaseFieldProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());

    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map_single_component);
    if (_use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic schme is not implemented for HydroPhaseField process");
    }
    else
    {
        // For displacement equation.
        const int process_id = 0;
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(process_id)[0].get().getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);

        // For phase field equation or hydro.
        _sparsity_pattern_with_single_component =
            NumLib::computeSparsityPattern(
                *_local_to_global_index_map_single_component, _mesh);
    }
}

template <int DisplacementDim>
void HydroPhaseFieldProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, HydroPhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaXZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXX));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonYY));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtEpsilonYZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtEpsilonXZ));
    }
}

template <int DisplacementDim>
void HydroPhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic schme is not implemented for HydroPhaseField process");
    }
    else
    {
        // for the equations of temperature-deformation.
        initializeProcessBoundaryConditionsAndSourceTerms(
            getDOFTableByProcessID(_mechanics_related_process_id),
            _mechanics_related_process_id);

        // for the hydro (pressure) field
        initializeProcessBoundaryConditionsAndSourceTerms(
            getDOFTableByProcessID(_hydro_process_id), _hydro_process_id);

        // for the phase field
        initializeProcessBoundaryConditionsAndSourceTerms(
            getDOFTableByProcessID(_phase_field_process_id),
            _phase_field_process_id);
    }
}

template <int DisplacementDim>
void HydroPhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble the equations for HydroPhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void HydroPhaseFieldProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG(
            "HydroPhaseFieldProcess for the monolithic scheme is not "
            "implemented.");
    }
    else
    {
        // For the staggered scheme
        if (_coupled_solutions->process_id == _mechanics_related_process_id)
        {
            DBUG(
                "Assemble the equations of deformation in "
                "HydroPhaseFieldProcess for the staggered scheme.");
        }
        else if (_coupled_solutions->process_id == _hydro_process_id)
        {
            DBUG(
                "Assemble the equations of pressure in "
                "HydroPhaseFieldProcess for the staggered scheme.");
        }
        else if (_coupled_solutions->process_id == _phase_field_process_id)
        {
            DBUG(
                "Assemble the equations of phase-field in "
                "HydroPhaseFieldProcess for the staggered scheme.");
        }
        dof_tables.emplace_back(
            getDOFTableByProcessID(_mechanics_related_process_id));
        dof_tables.emplace_back(getDOFTableByProcessID(_hydro_process_id));
        dof_tables.emplace_back(
            getDOFTableByProcessID(_phase_field_process_id));
    }

    // Hydro process
    if (_coupled_solutions->process_id == _hydro_process_id)
    {
        // calculate pressure
        // copy pressure to solution vector
    }
    else
    {
        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
            Jac, _coupled_solutions);
    }
}

template <int DisplacementDim>
void HydroPhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep HydroPhaseFieldProcess.");

    _process_data.dt = dt;
    _process_data.t = t;
    _process_data.injected_volume = _process_data.flow_rate * dt;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        getDOFTable(process_id), x, t, dt);
}

template <int DisplacementDim>
void HydroPhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, int const process_id)
{
    DBUG("PostTimestep HydroPhaseFieldProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        getDOFTable(process_id), x);
}

template <int DisplacementDim>
void HydroPhaseFieldProcess<
    DisplacementDim>::postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                                         const double t,
                                                         const int process_id)
{
    if (process_id == _mechanics_related_process_id)
    {
        double integral = 0;

        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;

        dof_tables.emplace_back(
            getDOFTableByProcessID(_mechanics_related_process_id));
        dof_tables.emplace_back(getDOFTableByProcessID(_hydro_process_id));
        dof_tables.emplace_back(
            getDOFTableByProcessID(_phase_field_process_id));

        DBUG("PostNonLinearSolver crack volume computation.");

        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::computeCrackIntegral, _local_assemblers,
            dof_tables, x, t, integral, _use_monolithic_scheme,
            _coupled_solutions);

        INFO("Integral of crack: %g", integral);
    }
    else if (process_id == _hydro_process_id)
    {
        // u = pressure * u_1
    }
}

}  // namespace HydroPhaseField
}  // namespace ProcessLib
