/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldDependentHydroBoundaryCondition.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>

namespace ProcessLib
{
void PhaseFieldDependentHydroBoundaryCondition::getEssentialBCValues(
    const double /*t*/, GlobalVector const& /*x*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    ParameterLib::SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(bc_values.ids.size() + _bc_values.ids.size());
    bc_values.values.reserve(bc_values.values.size() +
                             _bc_values.values.size());

    std::copy(_bc_values.ids.begin(), _bc_values.ids.end(),
              std::back_inserter(bc_values.ids));
    std::copy(_bc_values.values.begin(), _bc_values.values.end(),
              std::back_inserter(bc_values.values));
}
// update new values and corresponding indices.
void PhaseFieldDependentHydroBoundaryCondition::preTimestep(
    const double /*t*/, std::vector<GlobalVector*> const& x,
    int const process_id)
{
    // phase-field variable is considered irreversible if it loses more than 95%
    // of the stiffness, which is a widely used threshold.
    double activeDamage = 2.0;
    double currentValue = 0.0;

    _bc_values.ids.clear();
    _bc_values.values.clear();

    auto const mesh_id = _mesh.getID();
    auto const& nodes = _mesh.getNodes();
    for (auto const* n : nodes)
    {
        std::size_t node_id = n->getID();
        MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node, node_id);
        const auto g_idx =
            _dof_table.getGlobalIndex(l, _variable_id, _component_id);

        // For the DDC approach (e.g. with PETSc option), the negative
        // index of global_index means that the entry by that index is a ghost
        // one, which should be dropped. Especially for PETSc routines
        // MatZeroRows and MatZeroRowsColumns, which are called to apply the
        // Dirichlet BC, the negative index is not accepted like other matrix or
        // vector PETSc routines. Therefore, the following if-condition is
        // applied.
        if (g_idx < 0)
            continue;
        currentValue = (*x[process_id])[g_idx];

        if ((*x[3])[g_idx] > activeDamage)
        {
            _bc_values.ids.emplace_back(g_idx);
            _bc_values.values.emplace_back(currentValue);
        }
    }
}
// update new values and corresponding indices.
void PhaseFieldDependentHydroBoundaryCondition::postNonLinearSolver(
    const double /*t*/, GlobalVector const& x,
    int const process_id, CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{    

    double activeDamage = 2.0;
    double currentValue = 0.0;

    _bc_values.ids.clear();
    _bc_values.values.clear();

    auto const mesh_id = _mesh.getID();
    auto const& nodes = _mesh.getNodes();
    for (auto const* n : nodes)
    {
        std::size_t node_id = n->getID();
        MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node, node_id);
        const auto g_idx =
            _dof_table.getGlobalIndex(l, _variable_id, _component_id);

        // For the DDC approach (e.g. with PETSc option), the negative
        // index of global_index means that the entry by that index is a ghost
        // one, which should be dropped. Especially for PETSc routines
        // MatZeroRows and MatZeroRowsColumns, which are called to apply the
        // Dirichlet BC, the negative index is not accepted like other matrix or
        // vector PETSc routines. Therefore, the following if-condition is
        // applied.
        if (g_idx < 0)
            continue;
        currentValue = x[g_idx];

        if ((*cpl_xs->coupled_xs[3])[g_idx] > activeDamage)
        {
            _bc_values.ids.emplace_back(g_idx);
            _bc_values.values.emplace_back(currentValue);
        }
    }
}

std::unique_ptr<PhaseFieldDependentHydroBoundaryCondition>
createPhaseFieldDependentHydroBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    int const variable_id, int const component_id)
{
    DBUG(
        "Constructing PhaseFieldDependentHydroBoundaryCondition from "
        "config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter(
        "type", "PhaseFieldDependentHydroBoundaryCondition");

    return std::make_unique<
        PhaseFieldDependentHydroBoundaryCondition>(
        dof_table, mesh, variable_id, component_id);
}

}  // namespace ProcessLib
