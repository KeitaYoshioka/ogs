/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "GeoLib/AnalyticalGeometry.h"
#include "MeshLib/Mesh.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace HydroMechanicalPhaseField
{
struct HydroMechanicalPhaseFieldLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
    virtual std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual void computeFractureWidth(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        double const t, CoupledSolutionsForStaggeredScheme const* const cpl_xs,
        MeshLib::Mesh const& mesh) = 0;

    virtual void computeFractureNormal(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        CoupledSolutionsForStaggeredScheme const* const cpl_xs) = 0;

    /*    virtual void computeEnergy(
            std::size_t mesh_item_id,
            std::vector<
                std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
                dof_tables,
            GlobalVector const& x, double const t, double& elastic_energy,
            double& surface_energy, double& pressure_work,
            bool const use_monolithic_scheme,
            CoupledSolutionsForStaggeredScheme const* const cpl_xs) = 0;*/
};

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib
