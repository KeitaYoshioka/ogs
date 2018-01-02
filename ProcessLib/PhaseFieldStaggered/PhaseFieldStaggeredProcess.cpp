/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldStaggeredProcess.h"

#include <cassert>

#include "MeshLib/PropertyVector.h"

#include "CreateLocalAssemblers.h"
#include "PhaseFieldStaggeredLocalAssembler.h"
#include "PhaseFieldStaggeredProcessData.h"

namespace ProcessLib
{
namespace PhaseFieldStaggered
{
PhaseFieldStaggeredProcess::PhaseFieldStaggeredProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    PhaseFieldStaggeredProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    DBUG("Creat Phase field staggered process.")
}



void PhaseFieldStaggeredProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables()[0];
    createLocalAssemblers<PhaseFieldStaggeredLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        _coupling_term);
}

void PhaseFieldStaggeredProcess::assembleConcreteProcess(const double t,
                                                         GlobalVector const& x,
                                                         GlobalMatrix& M,
                                                         GlobalMatrix& K,
                                                         GlobalVector& b)
{
    DBUG("Assemble PhaseFieldStaggeredProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, _coupling_term);
}

void PhaseFieldStaggeredProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian PhaseFieldStaggeredProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, _coupling_term);
}

void PhaseFieldStaggeredProcess::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt)
{
    DBUG("PreTimestep PhaseFieldSmallDeformationProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        *_local_to_global_index_map, x, t, dt);
}

void PhaseFieldStaggeredProcess::postTimestepConcreteProcess(
    GlobalVector const& x)
{
    DBUG("PostTimestep PhaseFieldSmallDeformationProcess.");


    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        *_local_to_global_index_map, x);
}

void PhaseFieldStaggeredProcess::computeSecondaryVariableConcrete(
    const double t, GlobalVector const& x)
{
    //  DBUG("Compute damage gradient.");
    GlobalExecutor::executeMemberOnDereferenced(
        &PhaseFieldStaggeredLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, *_local_to_global_index_map, t, x, _coupling_term);
}

void PhaseFieldStaggeredProcess::setStaggeredCouplingTermToLocalAssemblers()
{
    GlobalExecutor::executeMemberOnDereferenced(
        &PhaseFieldStaggeredLocalAssemblerInterface::setStaggeredCouplingTerm,
        _local_assemblers, _coupling_term);
}

std::vector<double> PhaseFieldStaggeredProcess::getIntGradDamage(const std::size_t element_id) const
{
    return _local_assemblers[element_id]->getIntPtGradDamage();
}

}  // namespace PhaseFieldStaggered
}  // namespace ProcessLib
