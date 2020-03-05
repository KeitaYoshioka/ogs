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

#include "ParameterLib/Utils.h"

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

// void PhaseFieldDependentHydroBoundaryCondition::preTimestep(
//    const double /*t*/, std::vector<GlobalVector*> const& x,
//    int const process_id)
//{
//    // phase-field variable is considered irreversible if it loses more than
//    95%
//    // of the stiffness, which is a widely used threshold.

//    double currentValue = 0.0;

//    _bc_values.ids.clear();
//    _bc_values.values.clear();

//    std::vector<int> comp_nodes;
//    auto const& elems = _mesh.getElements();
//    auto const& width =
//        _mesh.getProperties().getPropertyVector<double>("width");

//    for (auto const* e : elems)
//    {
//        std::size_t elem_id = e->getID();
//        if ((*width)[elem_id] > 0.0)
//        {
//            for (int i = 0; i < e->getNumberOfNodes(); i++)
//            {
//                std::size_t node_id = e->getNode(i)->getID();
//                if (std::find(comp_nodes.begin(), comp_nodes.end(), node_id)
//                ==
//                    comp_nodes.end())
//                {
//                    comp_nodes.emplace_back(node_id);
//                }
//            }
//        }
//    }

//    for (auto const* e : elems)
//    {
//        std::size_t elem_id = e->getID();
//        if ((*width)[elem_id] < 1.e-20)
//        {
//            for (int i = 0; i < e->getNumberOfNodes(); i++)
//            {
//                std::size_t node_id = e->getNode(i)->getID();
//                currentValue = 1.e6;// (*x[process_id])[node_id];

//                if (std::find(_bc_values.ids.begin(), _bc_values.ids.end(),
//                              node_id) == _bc_values.ids.end() &&
//                    std::find(comp_nodes.begin(), comp_nodes.end(), node_id)
//                    ==
//                        comp_nodes.end())
//                {
//                    _bc_values.ids.emplace_back(node_id);
//                    _bc_values.values.emplace_back(currentValue);
//                }
//            }
//        }
//    }
//    /*
//        double activeDamage = 1.0;
//        auto const mesh_id = _mesh.getID();
//        auto const& nodes = _mesh.getNodes();
//        for (auto const* n : nodes)
//        {
//            std::size_t node_id = n->getID();
//            MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node,
//            node_id); const auto g_idx =
//                _dof_table.getGlobalIndex(l, _variable_id, _component_id);

//            if (g_idx < 0)
//                continue;
//            currentValue = (*x[process_id])[g_idx];

//            if ((*x[1])[g_idx] > activeDamage - 1.e-16)
//            {
//                _bc_values.ids.emplace_back(g_idx);
//                _bc_values.values.emplace_back(currentValue);
//            }
//        }*/
//    INFO("Pressure BC is updated.");
//}
// update new values and corresponding indices.
void PhaseFieldDependentHydroBoundaryCondition::postNonLinearSolver(
    const double t, GlobalVector const& /*x*/, int const /*process_id*/,
    CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    ParameterLib::SpatialPosition pos;
    double currentValue = 0.0;
    double activeDamage = 1.0;
    _bc_values.ids.clear();
    _bc_values.values.clear();
    std::vector<int> comp_nodes;
    auto const& elems = _mesh.getElements();
    auto const& width =
        _mesh.getProperties().getPropertyVector<double>("width");

    for (auto const* e : elems)
    {
        std::size_t elem_id = e->getID();
        if ((*width)[elem_id] > 0.0)
        {
            for (int i = 0; i < e->getNumberOfNodes(); i++)
            {
                std::size_t node_id = e->getNode(i)->getID();

                if (std::find(comp_nodes.begin(), comp_nodes.end(), node_id) ==
                        comp_nodes.end() &&
                    (*cpl_xs->coupled_xs[1])[node_id] < activeDamage - 1.e-16)
                    {
                        comp_nodes.emplace_back(node_id);
                    }
            }
        }
    }

    for (auto const* e : elems)
    {
        std::size_t elem_id = e->getID();
        if ((*width)[elem_id] < 1.e-20)
        {
            for (int i = 0; i < e->getNumberOfNodes(); i++)
            {
                std::size_t node_id = e->getNode(i)->getID();
                pos.setNodeID(node_id);
                currentValue = _parameter(t, pos)[0];

                if (std::find(_bc_values.ids.begin(), _bc_values.ids.end(),
                              node_id) == _bc_values.ids.end() &&
                    std::find(comp_nodes.begin(), comp_nodes.end(), node_id) ==
                        comp_nodes.end())
                {
                    _bc_values.ids.emplace_back(node_id);
                    _bc_values.values.emplace_back(currentValue);
                }
            }
        }
    }
    INFO("Pressure BC is updated.");

    /*
    double activeDamage = 1.0;
    auto const mesh_id = _mesh.getID();
    auto const& nodes = _mesh.getNodes();
    for (auto const* n : nodes)
    {
        std::size_t node_id = n->getID();
        MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node, node_id);
        const auto g_idx =
            _dof_table.getGlobalIndex(l, _variable_id, _component_id);

        if (g_idx < 0)
            continue;
        currentValue = x[g_idx];

        if ((*cpl_xs->coupled_xs[1])[g_idx] > activeDamage - 1.e-16)
        {
            _bc_values.ids.emplace_back(g_idx);
            _bc_values.values.emplace_back(currentValue);
        }
    }*/
}

std::unique_ptr<PhaseFieldDependentHydroBoundaryCondition>
createPhaseFieldDependentHydroBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    int const variable_id, int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters)
{
    DBUG(
        "Constructing PhaseFieldDependentHydroBoundaryCondition from "
        "config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type",
                                "PhaseFieldDependentHydroBoundaryCondition");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Dirichlet__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto& parameter =
        ParameterLib::findParameter<double>(param_name, parameters, 1, &mesh);
    return std::make_unique<PhaseFieldDependentHydroBoundaryCondition>(
        parameter, dof_table, mesh, variable_id, component_id);
}

}  // namespace ProcessLib
