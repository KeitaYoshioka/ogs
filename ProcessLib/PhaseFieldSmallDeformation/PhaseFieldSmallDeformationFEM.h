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
#include <vector>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/PhaseFieldExtension.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "PhaseFieldSmallDeformationProcessData.h"
#include "ProcessLib/PhaseFieldStaggered/PhaseFieldStaggeredProcess.h"

namespace ProcessLib
{
namespace PhaseFieldSmallDeformation
{
template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    typename BMatricesType::KelvinVectorType sigma_tensile, sigma_compressive,
        sigma_real_prev, sigma_real;
    double strain_energy_tensile;

    typename BMatricesType::KelvinVectorType eps, eps_prev;
    typename BMatricesType::KelvinMatrixType C_tensile, C_compressive;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    double d_ip, d_ip_prev;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_real_prev = sigma_real;
        material_state_variables->pushBackState();
        d_ip_prev = d_ip;
    }

    void updateConstitutiveRelation(double const t,
                                    SpatialPosition const& x_position,
                                    double const degradation)
    {
        static_cast<MaterialLib::Solids::PhaseFieldExtension<DisplacementDim>&>(
            solid_material)
            .calculateDegradedStress(t, x_position, eps, strain_energy_tensile,
                                     sigma_tensile, sigma_compressive,
                                     C_tensile, C_compressive, sigma_real,
                                     degradation);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class PhaseFieldSmallDeformationLocalAssembler
    : public PhaseFieldSmallDeformationLocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename BMatricesType::NodalForceVectorType;

    PhaseFieldSmallDeformationLocalAssembler(
        PhaseFieldSmallDeformationLocalAssembler const&) = delete;
    PhaseFieldSmallDeformationLocalAssembler(
        PhaseFieldSmallDeformationLocalAssembler&&) = delete;

    PhaseFieldSmallDeformationLocalAssembler(
        MeshLib::Element const& e, std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric, unsigned const integration_order,
        PhaseFieldSmallDeformationProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * sm.detJ;

            ip_data.N = sm.N;
            ip_data.dNdx = sm.dNdx;

            // Initialize current time step values
            ip_data.sigma_real.setZero(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.eps.setZero(KelvinVectorDimensions<DisplacementDim>::value);

            ip_data.C_tensile.setZero(
                KelvinVectorDimensions<DisplacementDim>::value,
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.C_compressive.setZero(
                KelvinVectorDimensions<DisplacementDim>::value,
                KelvinVectorDimensions<DisplacementDim>::value);

            ip_data.sigma_tensile.setZero(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.sigma_compressive.setZero(
                KelvinVectorDimensions<DisplacementDim>::value);

            ip_data.eps_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);

            _secondary_data.N[ip] = shape_matrices[ip].N;

            ip_data.d_ip = 1.0;
            ip_data.d_ip_prev = 1.0;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        OGS_FATAL(
            "PhaseFieldSmallDeformationLocalAssembler: assembly without "
            "jacobian is not "
            "implemented.");
    }

    void assembleWithJacobianAndCoupling(
        double const t, std::vector<double> const& local_x,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCouplingTerm const& coupled_term) override
    {
        for (auto const& coupled_process_pair : coupled_term.coupled_processes)
        {
            if (coupled_process_pair.first ==
                std::type_index(typeid(ProcessLib::PhaseFieldStaggered::
                                           PhaseFieldStaggeredProcess)))
            {
                auto const& pcs = static_cast<
                    ProcessLib::PhaseFieldStaggered::
                        PhaseFieldStaggeredProcess const&>(
                    coupled_process_pair.second);

                const auto local_d = coupled_term.local_coupled_xs.at(
                    coupled_process_pair.first);
                auto const grad_damage_ips = pcs.getIntGradDamage(_element.getID());

                assembleWithCoupledPhaseFieldJacobian(
                    t, local_x, local_Jac_data, local_b_data, local_d, grad_damage_ips);
            }
            else
            {
                OGS_FATAL(
                    "This coupled process is not presented for "
                    "PhaseFieldSmallDeformation process");
            }
        }
    }

    void assembleWithCoupledTerm(double const t,
                                 std::vector<double> const& local_x,
                                 std::vector<double>& /*local_M_data*/,
                                 std::vector<double>& local_K_data,
                                 std::vector<double>& local_b_data,
                                 LocalCouplingTerm const& coupled_term) override
    {
        for (auto const& coupled_process_pair : coupled_term.coupled_processes)
        {
            if (coupled_process_pair.first ==
                std::type_index(typeid(ProcessLib::PhaseFieldStaggered::
                                           PhaseFieldStaggeredProcess)))
            {
                auto const& pcs = static_cast<
                    ProcessLib::PhaseFieldStaggered::
                        PhaseFieldStaggeredProcess const&>(
                    coupled_process_pair.second);

                const auto local_d = coupled_term.local_coupled_xs.at(
                    coupled_process_pair.first);
                auto const grad_damage_ips = pcs.getIntGradDamage(_element.getID());

                assembleWithCoupledPhaseFieldJacobian(t, local_x, local_K_data,
                                                      local_b_data, local_d, grad_damage_ips);
            }
            else
            {
                OGS_FATAL(
                    "This coupled process is not presented for "
                    "PhaseFieldSmallDeformation process");
            }
        }
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        std::vector<double> local_dummy_d;
        local_dummy_d.resize(local_b_data.size());
        std::vector<double> dummy_grad_d;
        dummy_grad_d.resize(local_b_data.size());

        assembleWithCoupledPhaseFieldJacobian(t, local_x, local_Jac_data,
                                              local_b_data, local_dummy_d, dummy_grad_d);
    }

    void assembleWithCoupledPhaseFieldJacobian(
        double const t, std::vector<double> const& local_u,
        std::vector<double>& local_Jac_data, std::vector<double>& local_b_data,
        std::vector<double> const& local_d, std::vector<double> const& grad_damage_data)
    {
        auto const local_matrix_size = local_u.size();

        auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            local_b_data, local_matrix_size);

        auto local_grad_damage = Eigen::Matrix<double,ShapeFunction::NPOINTS, DisplacementDim>(
            grad_damage_data.data());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        assert(grad_damage_data.size() == n_integration_points*DisplacementDim);


        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

  //          GlobalDimVectorType const grad_damage = local_grad_damage.row(ip);

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u_op = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);
            for (int i = 0; i < DisplacementDim; ++i)
                N_u_op
                    .template block<1, displacement_size / DisplacementDim>(
                        i, i * displacement_size / DisplacementDim)
                    .noalias() = N;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);

            auto& eps = _ip_data[ip].eps;

            eps.noalias() =
                B *
                Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_u.data(), ShapeFunction::NPOINTS * DisplacementDim);
            // How to access to the material properties?
            double const k = 1.e-8;  //_process_data.residual_stiffness(t, x_position)[0];
            double d_ip = 0.;
            NumLib::shapeFunctionInterpolate(local_d, N, d_ip);
            double const degradation = d_ip * d_ip * (1 - k) + k;

            _ip_data[ip].updateConstitutiveRelation(t, x_position, /*local_u,*/
                                                    degradation);
            auto& sigma_real = _ip_data[ip].sigma_real;

            auto const& C_tensile = _ip_data[ip].C_tensile;
            auto const& C_compressive = _ip_data[ip].C_compressive;

            auto const rho = _process_data.solid_density(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;

            // TODO set these two in input file
            int const coupling_level = 1;
            bool const has_CrackPres = false;
            // TODO pressure is passed from outside of the process
            double const pres = 0.0;

            if (coupling_level == 1)
            {
                if(has_CrackPres)
                {
                    //int aa = N_u_op.transpose() * rho * b;
                    //int bb = N_u_op.transpose() * grad_damage;

                    local_b.noalias() -= (B.transpose() * sigma_real -
                                      N_u_op.transpose() * rho * b /*+ N_u_op.transpose() * grad_damage * pres*/) *
                                     w;
                }
                else{
                 local_b.noalias() -= (B.transpose() * sigma_real -
                                      N_u_op.transpose() * rho * b) *
                                     w;                   
                }

            }
            else if (coupling_level == 2)
            {
                auto& sigma_tensile = _ip_data[ip].sigma_tensile;
                _ip_data[ip].d_ip = d_ip;
                double d_ip_prev = _ip_data[ip].d_ip_prev;
                double delta_d_ip = d_ip - d_ip_prev;
                local_b.noalias() -=
                    (B.transpose() *
                         (sigma_real + 2 * d_ip * delta_d_ip * sigma_tensile) -
                     N_u_op.transpose() * rho * b) *
                    w;
            }

            local_Jac.noalias() += B.transpose() *
                                   (degradation * C_tensile + C_compressive) *
                                   B * w;
        }
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

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/, GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO make a general implementation for converting KelvinVectors
        // back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma_real;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size && component < 3;
                 ++component)
            {  // xx, yy, zz components
                cache_mat(component, ip) = sigma[component];
            }
            for (typename KelvinVectorType::Index component = 3;
                 component < kelvin_vector_size;
                 ++component)
            {  // mixed xy, yz, xz components
                cache_mat(component, ip) = sigma[component] / std::sqrt(2);
            }
        }

        return cache;
    }

    // TODO remove the component-wise methods.
    std::vector<double> const& getIntPtSigmaXX(
        const double /*t*/, GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 0);
    }

    std::vector<double> const& getIntPtSigmaYY(
        const double /*t*/, GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 1);
    }

    std::vector<double> const& getIntPtSigmaZZ(
        const double /*t*/, GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 2);
    }

    std::vector<double> const& getIntPtSigmaXY(
        const double /*t*/, GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 3);
    }

    std::vector<double> const& getIntPtSigmaYZ(
        const double /*t*/, GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 4);
    }

    std::vector<double> const& getIntPtSigmaXZ(
        const double /*t*/, GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 5);
    }

    std::vector<double> getIntPtStrainEnergyTensile() const override
    {
        auto const num_intpts = _ip_data.size();
        std::vector<double> st_result;

        st_result.resize(num_intpts);
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            st_result[ip] = _ip_data[ip].strain_energy_tensile;
        }

        return st_result;
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/, GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO make a general implementation for converting KelvinVectors
        // back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size && component < 3;
                 ++component)
            {  // xx, yy, zz components
                cache_mat(component, ip) = eps[component];
            }
            for (typename KelvinVectorType::Index component = 3;
                 component < kelvin_vector_size;
                 ++component)
            {  // mixed xy, yz, xz components
                cache_mat(component, ip) = eps[component] / std::sqrt(2);
            }
        }

        return cache;
    }

    // TODO remove the component-wise methods
    std::vector<double> const& getIntPtEpsilonXX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 0);
    }

    std::vector<double> const& getIntPtEpsilonYY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 1);
    }

    std::vector<double> const& getIntPtEpsilonZZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 2);
    }

    std::vector<double> const& getIntPtEpsilonXY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 3);
    }

    std::vector<double> const& getIntPtEpsilonYZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 4);
    }

    std::vector<double> const& getIntPtEpsilonXZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 5);
    }

    unsigned getNumberOfIntegrationPoints() const override
    {
        return _integration_method.getNumberOfPoints();
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override
    {
        return *_ip_data[integration_point].material_state_variables;
    }

private:
    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                             std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.sigma_real[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.sigma_real[component] / std::sqrt(2));
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        std::vector<double>& cache, std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.eps[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.eps[component] / std::sqrt(2));
        }

        return cache;
    }

    PhaseFieldSmallDeformationProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace PhaseFieldSmallDeformation
}  // namespace ProcessLib
