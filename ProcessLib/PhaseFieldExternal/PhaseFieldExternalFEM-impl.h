/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   PhaseFieldExternalFEM-impl.h
 */
#pragma once

#include "PhaseFieldExternalFEM.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace PhaseFieldExternal
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldExternalLocalAssembler<ShapeFunction, IntegrationMethod,
                                      DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    if (local_coupled_solutions.process_id == _phase_field_process_id)
    {
        assembleWithJacobianForPhaseFieldEquations(
            t, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
            local_b_data, local_Jac_data, local_coupled_solutions);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(
        t, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
        local_b_data, local_Jac_data, local_coupled_solutions);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldExternalLocalAssembler<ShapeFunction, IntegrationMethod,
                                      DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, std::vector<double> const& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    auto const& local_d =
        local_coupled_solutions.local_coupled_xs[_phase_field_process_id];
    auto const& local_u =
        local_coupled_solutions.local_coupled_xs[_mechanics_related_process_id];
    assert(local_d.size() == phasefield_size);
    assert(local_u.size() == displacement_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_u.data(), displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<displacement_size>>(
        local_b_data, displacement_size);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    double ele_d = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;
        ele_d += N * d;
    }

    ele_d = ele_d / n_integration_points;

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(_element,
                                                                     N);
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;

        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        double const k = _process_data.residual_stiffness(t, x_position)[0];
        double rho_sr = _process_data.solid_density(t, x_position)[0];
        double const alpha = _process_data.linear_thermal_expansion_coefficient(
            t, x_position)[0];
        double const beta = _process_data.biot_coefficient(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;
        double const T_ext = _process_data.temperature_ext(t, x_position)[0];
        double const T0 = _process_data.reference_temperature;
        double const P_ext = _process_data.pressure_ext(t, x_position)[0];
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        double const delta_T = T_ext - T0;
        // calculate real density
        double const rho_s = rho_sr / (1 + 3 * alpha * delta_T);

        double const degradation = ele_d * ele_d * (1 - k) + k;

        auto const C_eff = degradation * C_tensile + C_compressive;
        eps.noalias() = B * u;
        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u, alpha,
                                                delta_T, degradation);

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;

        local_Jac.noalias() += B.transpose() * C_eff * B * w;

        local_rhs.noalias() -=
            (B.transpose() * (sigma_eff - ele_d * beta * P_ext * identity2) -
             N_u.transpose() * rho_s * b - P_ext * N_u.transpose() * dNdx * d) *
            w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldExternalLocalAssembler<ShapeFunction, IntegrationMethod,
                                      DisplacementDim>::
    assembleWithJacobianForPhaseFieldEquations(
        double const t, std::vector<double> const& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    using DeformationVector =
        typename ShapeMatricesType::template VectorType<displacement_size>;
    using PhaseFieldVector =
        typename ShapeMatricesType::template VectorType<phasefield_size>;
    using PhaseFieldMatrix =
        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        phasefield_size>;

    auto const& local_u = local_coupled_solutions.local_coupled_xs[0];
    auto const& local_d = local_coupled_solutions.local_coupled_xs[1];
    assert(local_u.size() == displacement_size);
    assert(local_d.size() == phasefield_size);

    auto const local_matrix_size = local_d.size();
    auto d =
        Eigen::Map<PhaseFieldVector const>(local_d.data(), phasefield_size);
    auto u =
        Eigen::Map<DeformationVector const>(local_u.data(), displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, local_matrix_size);

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];
        double const P_ext = _process_data.pressure_ext(t, x_position)[0];

        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;

        auto& ip_data = _ip_data[ip];
        ip_data.strain_energy_tensile = strain_energy_tensile;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;

        // For AT2
        if (_process_data.at_param == 2)
        {
            local_Jac.noalias() +=
                (2 * N.transpose() * N * strain_energy_tensile +
                 gc * (N.transpose() * N / ls + dNdx.transpose() * dNdx * ls)) *
                w;

            local_rhs.noalias() -=
                (N.transpose() * N * d * 2 * strain_energy_tensile +
                 gc * ((N.transpose() * N / ls + dNdx.transpose() * dNdx * ls) *
                           d -
                       N.transpose() / ls) -
                 P_ext * dNdx.transpose() * N_u * u) *
                w;
        }
        // For AT1
        else
        {
            local_Jac.noalias() +=
                (2 * N.transpose() * N * strain_energy_tensile +
                 gc * (0.75 * dNdx.transpose() * dNdx * ls)) *
                w;

            local_rhs.noalias() -=
                (N.transpose() * N * d * 2 * strain_energy_tensile +
                 gc * (-0.375 * N.transpose() / ls +
                       0.75 * dNdx.transpose() * dNdx * ls * d) -
                 P_ext * dNdx.transpose() * N_u * u) *
                w;
        }
    }
}

}  // namespace PhaseFieldExternal
}  // namespace ProcessLib
