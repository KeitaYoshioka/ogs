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

#include "NumLib/DOF/DOFTableUtil.h"
#include "PhaseFieldExternalFEM.h"
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
        double const t, double const dt, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    if (process_id == _phase_field_process_id)
    {
        assembleWithJacobianForPhaseFieldEquations(
            t, dt, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
            local_b_data, local_Jac_data, local_coupled_solutions);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(
        t, dt, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
        local_b_data, local_Jac_data, local_coupled_solutions);
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldExternalLocalAssembler<ShapeFunction, IntegrationMethod,
                                      DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, double const dt,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    using DeformationMatrix =
        typename ShapeMatricesType::template MatrixType<_displacement_size,
                                                        _displacement_size>;

    assert(local_coupled_solutions.local_coupled_xs.size() ==
           _phasefield_size + _displacement_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[_phasefield_index],
        _phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_solutions.local_coupled_xs[_displacement_index],
        _displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<_displacement_size,
                                                        _displacement_size>>(
        local_Jac_data, _displacement_size, _displacement_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<_displacement_size>>(
        local_b_data, _displacement_size);

    double const& reg_param = _process_data.reg_param;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const d_ip = N.dot(d);

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(_element,
                                                                     N);
        auto const& B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx, N, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;

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

        double const degradation = d_ip * d_ip * (1 - k) + k;

        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, alpha, delta_T, degradation,
            _process_data.split_method, reg_param);

        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        auto const& sigma = _ip_data[ip].sigma;
        auto const C_eff = degradation * C_tensile + C_compressive;
        eps.noalias() = B * u;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        _displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, _displacement_size>::Zero(DisplacementDim,
                                                           _displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, _displacement_size / DisplacementDim>(
                   i, i * _displacement_size / DisplacementDim)
                .noalias() = N;

        local_Jac.noalias() += B.transpose() * C_eff * B * w;

        local_rhs.noalias() -=
            (B.transpose() * (sigma - d_ip * beta * P_ext * identity2) -
             N_u.transpose() * rho_s * b - P_ext * N_u.transpose() * dNdx * d) *
            w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldExternalLocalAssembler<ShapeFunction, IntegrationMethod,
                                      DisplacementDim>::
    assembleWithJacobianForPhaseFieldEquations(
        double const t, double const /*dt*/,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[_phasefield_index],
        _phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_solutions.local_coupled_xs[_displacement_index],
        _displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, _phasefield_size, _phasefield_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, _phasefield_size);

    ParameterLib::SpatialPosition x_position;
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
                                                        _displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, _displacement_size>::Zero(DisplacementDim,
                                                           _displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, _displacement_size / DisplacementDim>(
                   i, i * _displacement_size / DisplacementDim)
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
