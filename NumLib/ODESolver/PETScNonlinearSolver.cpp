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

#ifdef USE_PETSC

#include "PETScNonlinearSolver.h"

#include <petscmat.h>
#include <petscvec.h>

#include "BaseLib/RunTime.h"

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
PETScNonlinearSolver::PETScNonlinearSolver(GlobalLinearSolver&)
{
    SNESCreate(PETSC_COMM_WORLD, &_snes_solver);
    // SNESSetType(_snes_solver, "vi");
    SNESSetFromOptions(_snes_solver);
}

void PETScNonlinearSolver::setEquationSystem(System& eq,
                                             ConvergenceCriterion& conv_crit)
{
    _equation_system = &eq;
    _convergence_criterion = &conv_crit;
}

void PETScNonlinearSolver::assemble(GlobalVector const& x) const
{
    _equation_system->assemble(x);
}

bool PETScNonlinearSolver::solve(
    GlobalVector& x,
    std::function<void(unsigned, GlobalVector const&)> const&
    /*postIterationCallback*/,
    int const process_id)
{
    DBUG("PETScNonlinearSolver: solve()");
    using TimeDiscretizedSystem =
        TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                                 NonlinearSolverTag::Newton>;

    auto* system = static_cast<TimeDiscretizedSystem*>(_equation_system);

    DBUG("PETScNonlinearSolver: create vectors");
    // r and J on which the ogs assembly operates.
    auto& r = NumLib::GlobalVectorProvider::provider.getVector(
        system->getMatrixSpecifications(process_id), _residual_id);
    auto& J = NumLib::GlobalMatrixProvider::provider.getMatrix(
        system->getMatrixSpecifications(process_id), _jacobian_id);

    BaseLib::RunTime timer_dirichlet;
    double time_dirichlet = 0.0;

    timer_dirichlet.start();
    system->computeKnownSolutions(x);
    system->applyKnownSolutions(x);
    time_dirichlet += timer_dirichlet.elapsed();
    INFO("[time] Applying Dirichlet BCs took %g s.", time_dirichlet);

    // temporary r and J for petsc operations. These are copies of r and J
    // after the assembly.
    auto& petsc_r = NumLib::GlobalVectorProvider::provider.getVector(
        system->getMatrixSpecifications(process_id), _petsc_residual_id);
    auto& petsc_x = NumLib::GlobalVectorProvider::provider.getVector(
        system->getMatrixSpecifications(process_id), _petsc_x_id);
    VecCopy(x.getRawVector(), petsc_x.getRawVector());
    auto& petsc_J = NumLib::GlobalMatrixProvider::provider.getMatrix(
        system->getMatrixSpecifications(process_id), _petsc_jacobian_id);

    ::detail::PetscContext petsc_context{_equation_system, &x, &r, &J};

    auto residual = [](SNES /*snes*/, Vec petsc_x, Vec petsc_r,
                       void* petsc_context) -> PetscErrorCode {
        auto context = static_cast<::detail::PetscContext*>(petsc_context);

        DBUG("PETScNonlinearSolver: residual callback called.");

        // TODO Maybe need to overwrite the ogs x with petsc_x.
        /*

        // TODO (naumov) context->sys->preIteration(iteration, x);

        DBUG("BEFORE ASSEMBLY")
        DBUG("The ogs-x vector.")
        VecView(context->x->getRawVector(), PETSC_VIEWER_STDOUT_WORLD);
        DBUG("The petsc-x vector.")
        VecView(petsc_x, PETSC_VIEWER_STDOUT_WORLD);
        */
        VecCopy(petsc_x, context->x->getRawVector());

        // Assemble in ogs context.
        BaseLib::RunTime time_assembly;
        time_assembly.start();
        context->system->assemble(*context->x);

        INFO("[time] Assembly took %g s.", time_assembly.elapsed());
        context->system->getResidual(*context->x, *context->r);

        context->system->getJacobian(*context->J);
        context->system->applyKnownSolutionsNewton(
            *context->J, *context->r, *context->x);
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
        // don't call getJacobian, because it makes a copy of the system's J and
        // overwrites the context->J
        // context->system->getJacobian(*context->J);

        /*
        DBUG("The ogs-J matrix.")
        MatView(context->J->getRawMatrix(), PETSC_VIEWER_STDOUT_WORLD);
        */
        // Now copy the results into petsc vectors.
        MatCopy(context->J->getRawMatrix(), petsc_J, DIFFERENT_NONZERO_PATTERN);
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

    // Constraints
    DBUG("PETScNonlinearSolver: set constraints");
    auto& xl = NumLib::GlobalVectorProvider::provider.getVector(
        system->getMatrixSpecifications(process_id), _petsc_xl_id);
    auto& xu = NumLib::GlobalVectorProvider::provider.getVector(
        system->getMatrixSpecifications(process_id), _petsc_xu_id);
    VecSet(xl.getRawVector(), 0.0);
    VecSet(xu.getRawVector(), 1.0);

    system->updateConstraints(xl, xu);
    MathLib::finalizeVectorAssembly(xl);
    MathLib::finalizeVectorAssembly(xu);

    SNESVISetVariableBounds(_snes_solver, xl.getRawVector(), xu.getRawVector());

    // Solve
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

    NumLib::GlobalVectorProvider::provider.releaseVector(xl);
    NumLib::GlobalVectorProvider::provider.releaseVector(xu);

    return reason >= 0;
}

}  // namespace NumLib
#endif  // USE_PETSC
