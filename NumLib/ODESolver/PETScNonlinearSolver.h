/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#ifdef USE_PETSC

#include <petscsnes.h>

#include "ConvergenceCriterion.h"
#include "NonlinearSolver.h"
#include "TimeDiscretizedODESystem.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{
class PETScNonlinearSolver final : public NonlinearSolverBase
{
public:
    //! Type of the nonlinear equation system to be solved.
    using System = NonlinearSystem<NonlinearSolverTag::Newton>;

    /*! Constructs a new instance.
     *
     * \param linear_solver the linear solver used by this nonlinear solver.
     */
    explicit PETScNonlinearSolver(GlobalLinearSolver& linear_solver);

    ~PETScNonlinearSolver() = default;

    //! Set the nonlinear equation system that will be solved.
    //! TODO doc
    void setEquationSystem(System& eq, ConvergenceCriterion& conv_crit);

    void assemble(GlobalVector const& x) const override;

    NonlinearSolverStatus solve(GlobalVector& x,
               std::function<void(int, GlobalVector const&)> const&
               /*postIterationCallback*/, int process_id) override;

private:
    /* TODO (naumov) remove or use inside
    GlobalLinearSolver& _linear_solver;
    */

    SNES _snes_solver;

    System* _equation_system = nullptr;
    ConvergenceCriterion* _convergence_criterion = nullptr;

    std::size_t _residual_id = 0u;  //!< ID of the residual vector.
    std::size_t _jacobian_id = 0u;  //!< ID of the Jacobian matrix.

    std::size_t _petsc_residual_id = 0u;
    std::size_t _petsc_x_id = 0u;
    std::size_t _petsc_jacobian_id = 0u;

    std::size_t _petsc_xl_id = 0u;
    std::size_t _petsc_xu_id = 0u;
};

//! @}

}  // namespace NumLib
#endif  // USE_PETSC
