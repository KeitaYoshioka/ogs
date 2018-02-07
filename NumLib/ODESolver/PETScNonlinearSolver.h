/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#ifdef USE_PETSC

#include <petscmat.h>
#include <petscsnes.h>
#include <petscvec.h>
#include <memory>

#include "BaseLib/RunTime.h"

#include "ConvergenceCriterion.h"
#include "NonlinearSolver.h"
#include "TimeDiscretizedODESystem.h"

namespace BaseLib
{
class ConfigTree;
}
namespace detail
{
struct PetscContext
{
    using System = NumLib::NonlinearSystem<NumLib::NonlinearSolverTag::Newton>;
    System* system;
    GlobalVector* x;
    GlobalVector* r;
    GlobalMatrix* J;
};
}  // namespace detail

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
    explicit PETScNonlinearSolver(GlobalLinearSolver& linear_solver)
        : _linear_solver(linear_solver)
    {
        SNESCreate(PETSC_COMM_WORLD, &_snes_solver);
        // SNESSetType(_snes_solver, "vi");
        SNESSetFromOptions(_snes_solver);
    }

    ~PETScNonlinearSolver() = default;

    //! Set the nonlinear equation system that will be solved.
    //! TODO doc
    void setEquationSystem(System& eq, ConvergenceCriterion& conv_crit)
    {
        _equation_system = &eq;
        _convergence_criterion = &conv_crit;
    }

    void assemble(GlobalVector const& x) const override
    {
        _equation_system->assemble(x);
    }

    bool solve(GlobalVector& x,
               std::function<void(unsigned, GlobalVector const&)> const&
               /*postIterationCallback*/) override
    {
        DBUG("PETScNonlinearSolver: solve()");
        using TimeDiscretizedSystem = TimeDiscretizedODESystem<
            ODESystemTag::FirstOrderImplicitQuasilinear,
            NonlinearSolverTag::Newton>;

        auto* system = static_cast<TimeDiscretizedSystem*>(_equation_system);

        DBUG("PETScNonlinearSolver: create vectors");
        // r and J on which the ogs assembly operates.
        auto& r = NumLib::GlobalVectorProvider::provider.getVector(
            system->getMatrixSpecifications(1), _residual_id);
        auto& J = NumLib::GlobalMatrixProvider::provider.getMatrix(
            system->getMatrixSpecifications(1), _jacobian_id);

        // temporary r and J for petsc operations. These are copies of r and J
        // after the assembly.
        auto& petsc_r = NumLib::GlobalVectorProvider::provider.getVector(
            system->getMatrixSpecifications(1), _petsc_residual_id);
        auto& petsc_x = NumLib::GlobalVectorProvider::provider.getVector(
            system->getMatrixSpecifications(1), _petsc_x_id);
        VecCopy(x.getRawVector(), petsc_x.getRawVector());
        auto& petsc_J = NumLib::GlobalMatrixProvider::provider.getMatrix(
            system->getMatrixSpecifications(1), _petsc_jacobian_id);

        ::detail::PetscContext petsc_context{_equation_system, &x, &r, &J};

        auto residual = [](SNES /*snes*/, Vec /*petsc_x*/, Vec petsc_r,
                           void* petsc_context) -> PetscErrorCode {
            auto context = static_cast<::detail::PetscContext*>(petsc_context);

            DBUG("PETScNonlinearSolver: residual callback called.");

            // TODO Maybe need to overwrite the ogs x with petsc_x.
            /*
            DBUG("BEFORE ASSEMBLY")
            DBUG("The ogs-x vector.")
            VecView(context->x->getRawVector(), PETSC_VIEWER_STDOUT_WORLD);
            DBUG("The petsc-x vector.")
            VecView(petsc_x, PETSC_VIEWER_STDOUT_WORLD);
            */

            // Assemble in ogs context.
            BaseLib::RunTime time_assembly;
            time_assembly.start();
            context->system->assemble(*context->x);

            INFO("[time] Assembly took %g s.", time_assembly.elapsed());
            context->system->getResidual(*context->x, *context->r);

            context->system->getJacobian(*context->J);
            context->system->applyKnownSolutionsNewton(*context->J, *context->r,
                                                       *context->x);
            /*
            DBUG("AFTER ASSEMBLY");
            DBUG("The ogs-x vector.")
            VecView(context->x->getRawVector(), PETSC_VIEWER_STDOUT_WORLD);
            DBUG("The ogs-r vector.")
            VecView(context->r->getRawVector(), PETSC_VIEWER_STDOUT_WORLD);
            */

            // Now copy the results into petsc vectors.
            VecCopy(context->r->getRawVector(), petsc_r);

            /*
            DBUG("The petsc-x vector.")
            VecView(petsc_x, PETSC_VIEWER_STDOUT_WORLD);
            DBUG("The petsc-r vector.")
            VecView(petsc_r, PETSC_VIEWER_STDOUT_WORLD);
            */
            return 0;
        };

        auto jacobian = [](SNES /*snes*/,
                           Vec /*petsc_x*/,
                           Mat petsc_J,
                           Mat /*petsc_B* Same as petsc_J*/,
                           void* petsc_context) {
            DBUG("PETScNonlinearSolver: jacobian callback called.")
            // Assume the system is already assembled.
            auto context = static_cast<::detail::PetscContext*>(petsc_context);
            context->system->getJacobian(*context->J);

            /*
            DBUG("The ogs-J matrix.")
            MatView(context->J->getRawMatrix(), PETSC_VIEWER_STDOUT_WORLD);
            */
            // Now copy the results into petsc vectors.
            MatCopy(context->J->getRawMatrix(), petsc_J,
                    DIFFERENT_NONZERO_PATTERN);
            /*
            DBUG("The petsc-J matrix.")
            MatView(petsc_J, PETSC_VIEWER_STDOUT_WORLD);
            */
            return 0;
        };

        DBUG("PETScNonlinearSolver: set function");
        SNESSetFunction(_snes_solver, petsc_r.getRawVector(), residual,
                        &petsc_context);

        DBUG("PETScNonlinearSolver: set jacobian");
        // The jacobian and the preconditioner matrices point to the same
        // location.
        SNESSetJacobian(_snes_solver, petsc_J.getRawMatrix(),
                        petsc_J.getRawMatrix(), jacobian, &petsc_context);

        DBUG("PETScNonlinearSolver: set constraints");
        // Constraints
        Vec xl, xu;
        VecDuplicate(x.getRawVector(), &xl);
        VecDuplicate(x.getRawVector(), &xu);
        VecSet(xl, 0.0);
        VecSet(xu, 1.0);

        SNESVISetVariableBounds(_snes_solver, xl, xu);

        DBUG("PETScNonlinearSolver: call SNESSolve");
        SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
        PetscInt iterations;

        SNESSolve(_snes_solver, nullptr, petsc_x.getRawVector());

        SNESGetConvergedReason(_snes_solver, &reason);
        SNESGetIterationNumber(_snes_solver, &iterations);
        INFO("PETScSNES used %d iterations.", iterations);
        INFO("PETSsSNES reason %d.", reason);

        // Copy back the solution.
        VecCopy(petsc_x.getRawVector(), x.getRawVector());
        // XXX is copying of r and J necessary?
        NumLib::GlobalMatrixProvider::provider.releaseMatrix(J);
        NumLib::GlobalVectorProvider::provider.releaseVector(r);
        NumLib::GlobalMatrixProvider::provider.releaseMatrix(petsc_J);
        NumLib::GlobalVectorProvider::provider.releaseVector(petsc_r);
        NumLib::GlobalVectorProvider::provider.releaseVector(petsc_x);

        VecDestroy(&xl);
        VecDestroy(&xu);

        return reason >= 0;
    }

private:
    GlobalLinearSolver& _linear_solver;

    SNES _snes_solver;

    System* _equation_system = nullptr;
    ConvergenceCriterion* _convergence_criterion = nullptr;

    std::size_t _residual_id = 0u;  //!< ID of the residual vector.
    std::size_t _jacobian_id = 0u;  //!< ID of the Jacobian matrix.

    std::size_t _petsc_residual_id = 0u;
    std::size_t _petsc_x_id = 0u;
    std::size_t _petsc_jacobian_id = 0u;
};

//! @}

}  // namespace NumLib
#endif  // USE_PETSC
