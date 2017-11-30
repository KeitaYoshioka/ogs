/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowLocalAssembler.h
 *
 * Created on August 19, 2016, 2:28 PM
 */

#pragma once

#include <map>
#include <typeindex>
#include <unordered_map>
#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "PhaseFieldStaggeredProcessData.h"

namespace ProcessLib
{
struct StaggeredCouplingTerm;

namespace PhaseFieldStaggered
{
const unsigned NUM_NODAL_DOF = 1;

struct IntegrationPointData final
{
    double strain_energy_tensile = 0;
    double strain_energy_tensile_prev = 0;
    double delta_strain_energy_tensile = 0;
    double history_variable = 0;
    double history_variable_prev = 0;

    void pushBackState()
    {
        if (history_variable_prev < history_variable)
        {
            history_variable_prev = history_variable;
        }
        delta_strain_energy_tensile =
            strain_energy_tensile - strain_energy_tensile_prev;
        strain_energy_tensile_prev = strain_energy_tensile;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class PhaseFieldStaggeredLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    PhaseFieldStaggeredLocalAssemblerInterface(
        StaggeredCouplingTerm* const coupling_term)
        : _coupling_term(coupling_term)
    {
    }

    void setStaggeredCouplingTerm(std::size_t const /*mesh_item_id*/,
                                  StaggeredCouplingTerm* const coupling_term)
    {
        _coupling_term = coupling_term;
    }

#if 0
    virtual std::vector<double> const& getIntPtGradDamageX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtGradDamageY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtGradDamageZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;
#endif

protected:
    // TODO: remove _coupling_term or move integration point data from local
    // assembler class to a new class to make local assembler unique for each
    // process.
    /// Pointer that is set from a Process class.
    StaggeredCouplingTerm* _coupling_term;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class PhaseFieldStaggeredLocalAssembler
    : public PhaseFieldStaggeredLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

    using MatrixOfVelocityAtIntegrationPoints = Eigen::Map<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>;

    using RhsVector = typename ShapeMatricesType::template VectorType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * GlobalDim>;

public:
    PhaseFieldStaggeredLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric, unsigned const integration_order,
        PhaseFieldStaggeredProcessData& process_data,
        StaggeredCouplingTerm* coupling_term)
        : PhaseFieldStaggeredLocalAssemblerInterface(coupling_term),
          _element(element),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _process_data(process_data)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.resize(n_integration_points);
    }

    void assemble(double const, std::vector<double> const&,
                  std::vector<double>&, std::vector<double>&,
                  std::vector<double>&) override
    {
        OGS_FATAL(
            "PhaseFieldSmallDeformationLocalAssembler: assembly without "
            "jacobian is not "
            "implemented.");
    }

    void assembleWithCoupledTerm(
        double const t,
        std::vector<double> const& local_x,
        std::vector<double>& local_M_data,
        std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        LocalCouplingTerm const& coupled_term) override;

    void assembleWithJacobianAndCoupling(
        double const t, std::vector<double> const& local_x,
        std::vector<double> const& local_xdot, const double dxdot_dx,
        const double dx_dx, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data,
        LocalCouplingTerm const& coupling_term) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

#if 0
    std::vector<double> const& getIntPtGradDamageX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override;

    std::vector<double> const& getIntPtGradDamageY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override;

    std::vector<double> const& getIntPtGradDamageZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override;
#endif
private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;

    PhaseFieldStaggeredProcessData& _process_data;

    void assembleWithCoupledPhaseFieldStaggered(
        double const t, std::vector<double> const& local_x,
        std::vector<double> const& strain_energy_tensile_ips,
        std::vector<double>& local_Jac_data,
        std::vector<double>& local_rhs_data);
    std::vector<IntegrationPointData,
                Eigen::aligned_allocator<IntegrationPointData>>
        _ip_data;
};

}  // end of namespace
}  // end of namespace

#include "PhaseFieldStaggeredLocalAssembler-impl.h"
