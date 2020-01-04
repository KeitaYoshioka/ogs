/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "DPHMPhaseFieldProcessData.h"
#include "LocalAssemblerInterface.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace DPHMPhaseField
{
struct DPHMPhaseFieldLocalAssemblerInterface;

template <int DisplacementDim>
class DPHMPhaseFieldProcess final : public Process
{
public:
    DPHMPhaseFieldProcess(
        std::string name, MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        DPHMPhaseFieldProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        int const mechanics_related_process_id,
        int const phase_field_process_id, int const frac_hydro_process_id,
        int const pore_hydro_process);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

private:
    void constructDofTable() override;

    void initializeBoundaryConditions() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 GlobalVector const& x, int const process_id,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, GlobalVector const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
                                    const int process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     double const t, double const dt,
                                     int const process_id) override;

    void postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                            const double t,
                                            double const dt,
                                            int const process_id) override;
    void updateConstraints(GlobalVector& lower, GlobalVector& upper) override;

    // To be replaced.
    NumLib::LocalToGlobalIndexMap& getDOFTableByProcessID(
        const int process_id) const;

    void setCoupledSolutionsOfPreviousTimeStepPerProcess(const int process_id);

    /// Set the solutions of the previous time step to the coupled term.
    /// It is only for the staggered scheme, and it must be called within
    /// the coupling loop because that the coupling term is only created there.
    void setCoupledSolutionsOfPreviousTimeStep();

private:
    DPHMPhaseFieldProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<DPHMPhaseFieldLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;

    /// Previous time step solution used for the constraints.
    std::array<std::unique_ptr<GlobalVector>, 4> _x_previous_timestep;

    /// Sparsity pattern for the phase field equation, and it is initialized
    //  only if the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_single_component;

    /// ID of the processes that contains mechanical process.
    int const _mechanics_related_process_id;

    /// ID of phase field process.
    int const _phase_field_process_id;

    /// ID of hydro process in fracture.
    int const _frac_hydro_process_id;

    /// ID of hydro process in pore.
    int const _pore_hydro_process_id;
};

extern template class DPHMPhaseFieldProcess<2>;
extern template class DPHMPhaseFieldProcess<3>;

}  // namespace DPHMPhaseField
}  // namespace ProcessLib
