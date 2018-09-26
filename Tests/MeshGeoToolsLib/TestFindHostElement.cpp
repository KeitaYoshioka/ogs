/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>
#include <ctime>
#include <memory>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"
#include "ProcessLib/HydroMechanicalPhaseField/HydroMechanicalPhaseFieldFEM.h"
#include "Tests/AutoCheckTools.h"

using namespace MeshLib;
using namespace ProcessLib::HydroMechanicalPhaseField;
namespace ac = autocheck;

struct FindHostElement : public ::testing::Test
{
    ac::IntervalGenerator<double> x_gen{0, 10};
    ac::IntervalGenerator<double> y_gen{0, 10};

    // ac::randomTupleGenerator<double, 2> tuple_gen;//(x_gen);

    ac::gtest_reporter gtest_reporter;
};

std::size_t elementIdAtPoint(Eigen::Vector3d const& pnt)
{
    auto const column = static_cast<int>(std::floor(pnt[0]));
    auto const row = static_cast<int>(std::floor(pnt[1]));
    return row * 10 + column;
}

TEST_F(FindHostElement, HostElementQuadSameElement)
{
    auto mesh =
        std::unique_ptr<Mesh>(MeshGenerator::generateRegularQuadMesh(10., 10));

    auto same_element_returned = [&mesh](double const x,
                                         double const y) -> bool {

        Eigen::Vector3d const pnt(x, y, 0);
        auto const ref_ele_id = elementIdAtPoint(pnt);

        double const probe_offset = 10.;

        MeshLib::Element const* ref_ele = mesh->getElement(ref_ele_id);
        assert(ref_ele != nullptr);
        MeshLib::Element const* host_ele = nullptr;
        findHostElement(*ref_ele, pnt, host_ele, probe_offset);
        std::size_t const expected_ele_id = ref_ele_id;
        assert(host_ele != nullptr);
        return expected_ele_id == host_ele->getID();
    };

    ac::check<double, double>(same_element_returned, 100,
                              ac::make_arbitrary(x_gen, y_gen), gtest_reporter);
}

TEST_F(FindHostElement, HostElementQuad)
{
    auto mesh =
        std::unique_ptr<Mesh>(MeshGenerator::generateRegularQuadMesh(10., 10));

    auto same_element_returned = [&mesh](double const x,
                                         double const y) -> bool {

        Eigen::Vector3d const pnt(x, y, 0);
        double const probe_offset = 10.;
        auto const expected_ele_id = elementIdAtPoint(pnt);
        MeshLib::Element const* expected_ele =
            mesh->getElement(expected_ele_id);

        std::size_t num_neighbor = expected_ele->getNumberOfNeighbors();
        for (std::size_t i = 0; i < num_neighbor; i++)
        {
            MeshLib::Element const* ref_ele = expected_ele->getNeighbor(i);
            if (ref_ele == nullptr)
                continue;
            MeshLib::Element const* host_ele = nullptr;
            findHostElement(*ref_ele, pnt, host_ele, probe_offset);
            assert(host_ele != nullptr);
            if (expected_ele_id != host_ele->getID())
                return false;
        }

        return true;
    };

    ac::check<double, double>(same_element_returned, 100,
                              ac::make_arbitrary(x_gen, y_gen), gtest_reporter);
}

TEST_F(FindHostElement, HostElementQuadCorner)
{
    auto mesh =
        std::unique_ptr<Mesh>(MeshGenerator::generateRegularQuadMesh(10., 10));

    auto same_element_returned = [&mesh](std::size_t const node_id) -> bool {

        auto const& node = *mesh->getNode(node_id);
        Eigen::Vector3d const corner_pnt(node[0], node[1], node[2]);
        double const probe_offset = 10.;

        for (auto const* element : node.getElements())
        {
            if (element == nullptr)
                continue;
            MeshLib::Element const* host_ele = nullptr;
            findHostElement(*element, corner_pnt, host_ele, probe_offset);
            assert(host_ele != nullptr);
            if (element->getID() != host_ele->getID())
                return false;
        }
        return true;
    };
    ac::IntervalGenerator<std::size_t> node_id_gen{0, mesh->getNumberOfNodes()};
    ac::check<std::size_t>(same_element_returned, 100,
                           ac::make_arbitrary(node_id_gen), gtest_reporter);
}

TEST_F(FindHostElement, HostElementQuadEdge)
{
    auto mesh =
        std::unique_ptr<Mesh>(MeshGenerator::generateRegularQuadMesh(10., 10));

    auto same_element_returned = [&mesh](std::size_t const ele_id) -> bool {

        double const probe_offset = 10.;
        auto const* ref_ele = mesh->getElement(ele_id);
        int num_edge = ref_ele->getNumberOfEdges();
        for (int i = 0; i < num_edge; i++)
        {
            auto edge_ele = ref_ele->getEdge(i);
            auto n0 = *edge_ele->getNode(0);
            auto n1 = *edge_ele->getNode(1);
            MeshLib::Element const* host_ele = nullptr;
            Eigen::Vector3d const edge_pnt0(0.5 * (n0[0] + n1[0]), n0[1],
                                             n0[2]);
            findHostElement(*ref_ele, edge_pnt0, host_ele, probe_offset);
            assert(host_ele != nullptr);
            if (ref_ele->getID() != host_ele->getID())
                return false;
            Eigen::Vector3d const edge_pnt1(n0[0], 0.5 *(n0[1] + n1[1]),
                                             n0[2]);
            findHostElement(*ref_ele, edge_pnt1, host_ele, probe_offset);
            assert(host_ele != nullptr);
            if (ref_ele->getID() != host_ele->getID())
                return false;
        }

        return true;
    };
    ac::IntervalGenerator<std::size_t> ele_id_gen{0,
                                                  mesh->getNumberOfElements()-1};

    ac::IntervalGenerator<double> frac_gen{0,1};
    ac::check<std::size_t>(same_element_returned, 100,
                           ac::make_arbitrary(ele_id_gen), gtest_reporter);
}
