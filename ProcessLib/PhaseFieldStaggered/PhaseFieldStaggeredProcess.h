/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"
#include "PhaseFieldStaggeredProcessData.h"

namespace ProcessLib {
namespace PhaseFieldStaggered {
class PhaseFieldStaggeredLocalAssemblerInterface;
struct PhaseFieldStaggeredProcessData;

class PhaseFieldStaggeredProcess final : public Process {
public:
  PhaseFieldStaggeredProcess(
      MeshLib::Mesh& mesh,
      std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
      std::vector<std::unique_ptr<ParameterBase>> const& parameters,
      unsigned const integration_order,
      std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
      PhaseFieldStaggeredProcessData&& process_data,
      SecondaryVariableCollection&& secondary_variables,
      NumLib::NamedFunctionCaller&& named_function_caller);

  bool isLinear() const override { return true; };

  void computeSecondaryVariableConcrete(double const t,
                                        GlobalVector const& x) override;
  void setStaggeredCouplingTermToLocalAssemblers() override;

private:
  void initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const& dof_table,
                                 MeshLib::Mesh const& mesh,
                                 unsigned const integration_order) override;

  void assembleConcreteProcess(const double t, GlobalVector const& x,
                               GlobalMatrix& M, GlobalMatrix& K,
                               GlobalVector& b) override;

  void assembleWithJacobianConcreteProcess(
      const double t, GlobalVector const& x, GlobalVector const& xdot,
      const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
      GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

  void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                  double const dt) override;

  PhaseFieldStaggeredProcessData _process_data;
//  const std::unique_ptr<PhaseFieldStaggeredProcessData> _process_data;
  std::vector<std::unique_ptr<PhaseFieldStaggeredLocalAssemblerInterface>>
      _local_assemblers;
};

} // namespace PhaseFieldStaggered
} // namespace ProcessLib
