/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on January 8, 2018, 3:00 PM
 */
#pragma once

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace PhaseField
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    // For the equations with phase field.
    if (process_id == 1)
    {
        assembleWithJacobianPhaseFiledEquations(
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
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
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
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>;

    assert(local_coupled_solutions.local_coupled_xs.size() ==
           phasefield_size + displacement_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[phasefield_index],
        phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_solutions.local_coupled_xs[displacement_index],
        displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<DeformationMatrix>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<DeformationVector>(
        local_b_data, displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    double const& reg_param = _process_data.reg_param;

    auto local_pressure = 0.0;
    if (_process_data.crack_pressure)
    {
        local_pressure = _process_data.unity_pressure;
    }

    int const n_integration_points = _integration_method.getNumberOfPoints();

    double ele_d = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;
        ele_d += N.dot(d);
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
        eps.noalias() = B * u;
        double const k = _process_data.residual_stiffness(t, x_position)[0];

        double degradation;
        // KKL
        if (_process_data.at_param == 3)
            degradation = (4 * pow(ele_d, 3) - 3 * pow(ele_d, 4)) * (1 - k) + k;
        // ATn
        else
            degradation = ele_d * ele_d * (1 - k) + k;

        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation, _process_data.split_method, reg_param);

        auto& sigma = _ip_data[ip].sigma;
        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;

        local_rhs.noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho_sr * b -
             local_pressure * N_u.transpose() * dNdx * d) *
            w;

        local_Jac.noalias() +=
            B.transpose() * (degradation * C_tensile + C_compressive) * B * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianPhaseFiledEquations(
        double const t, double const dt,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[phasefield_index],
        phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_solutions.local_coupled_xs[displacement_index],
        displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, phasefield_size, phasefield_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, phasefield_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    double const& reg_param = _process_data.reg_param;

    auto local_pressure = 0.0;
    if (_process_data.crack_pressure)
    {
        local_pressure = _process_data.unity_pressure;
    }
    else if (_process_data.propagating_crack)
    {
        local_pressure = _process_data.pressure;
    }

    int const n_integration_points = _integration_method.getNumberOfPoints();
    double ele_d = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;
        ele_d += N.dot(d);
    }
    ele_d = ele_d / n_integration_points;

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];

        // for propagating crack, u is rescaled.
        if (_process_data.propagating_crack)
        {
            double const k = _process_data.residual_stiffness(t, x_position)[0];
            double degradation;
            // KKL
            if (_process_data.at_param == 3)
                degradation =
                    (4 * pow(ele_d, 3) - 3 * pow(ele_d, 4)) * (1 - k) + k;
            // ATn
            else
                degradation = ele_d * ele_d * (1 - k) + k;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const& B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);

            auto& eps = _ip_data[ip].eps;
            eps.noalias() = B * u;
            _ip_data[ip].updateConstitutiveRelation(
                t, x_position, dt, u, degradation, _process_data.split_method, reg_param);
        }

        auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;

        auto& ip_data = _ip_data[ip];
        ip_data.strain_energy_tensile = strain_energy_tensile;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

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
                 local_pressure * dNdx.transpose() * N_u * u) *
                w;
        }
        // For AT1
        else if (_process_data.at_param == 1)
        {
            local_Jac.noalias() +=
                (2 * N.transpose() * N * strain_energy_tensile +
                 gc * (0.75 * dNdx.transpose() * dNdx * ls)) *
                w;

            local_rhs.noalias() -=
                (N.transpose() * N * d * 2 * strain_energy_tensile +
                 gc * (-0.375 * N.transpose() / ls +
                       0.75 * dNdx.transpose() * dNdx * ls * d) -
                 local_pressure * dNdx.transpose() * N_u * u) *
                w;
        }
        // For KKL
        else
        {
            double f, dfdv;
            double Cv = 0.7165753;
            double coef = strain_energy_tensile - gc / (4 * Cv * ls);

            f = 12 * ele_d * ele_d * (1 - ele_d) * coef;

            dfdv = 12 * ele_d * (2 - 3 * ele_d) * coef;

            local_Jac.noalias() +=
                (dfdv * N.transpose() * N +
                 ls * gc / (2 * Cv) * (dNdx.transpose() * dNdx)) *
                w;

            local_rhs.noalias() -=
                ((gc * ls / (2 * Cv) * dNdx.transpose() * dNdx) * d +
                 f * N.transpose() -
                 local_pressure * dNdx.transpose() * N_u * u) *
                w;
        }
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    computeCrackIntegral(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        GlobalVector const& /*x*/, double const /*t*/, double& crack_volume,
        CoupledSolutionsForStaggeredScheme const* const cpl_xs) const
{
    assert(cpl_xs != nullptr);

    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    std::transform(dof_tables.begin(), dof_tables.end(),
                   std::back_inserter(indices_of_processes),
                   [&](NumLib::LocalToGlobalIndexMap const& dof_table) {
                       return NumLib::getIndices(mesh_item_id, dof_table);
                   });

    auto local_coupled_xs =
        getCoupledLocalSolutions(cpl_xs->coupled_xs, indices_of_processes);
    assert(local_coupled_xs.size() == displacement_size + phasefield_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    double ele_crack_vol = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

        ele_crack_vol += (N_u * u).dot(dNdx * d) * w;
    }
#ifdef USE_PETSC
    int const n_all_nodes = indices_of_processes[1].size();
    int const n_regular_nodes = std::count_if(
        begin(indices_of_processes[1]), end(indices_of_processes[1]),
        [](GlobalIndexType const& index) { return index >= 0; });
    if (n_all_nodes != n_regular_nodes)
    {
        ele_crack_vol *= static_cast<double>(n_regular_nodes) / n_all_nodes;
    }
#endif  // USE_PETSC
    crack_volume += ele_crack_vol;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    computeEnergy(std::size_t mesh_item_id,
                  std::vector<std::reference_wrapper<
                      NumLib::LocalToGlobalIndexMap>> const& dof_tables,
                  GlobalVector const& /*x*/, double const t,
                  double& elastic_energy, double& surface_energy,
                  double& pressure_work,
                  CoupledSolutionsForStaggeredScheme const* const cpl_xs)
{
    assert(cpl_xs != nullptr);

    std::vector<std::vector<GlobalIndexType>> indices_of_processes;
    indices_of_processes.reserve(dof_tables.size());
    std::transform(dof_tables.begin(), dof_tables.end(),
                   std::back_inserter(indices_of_processes),
                   [&](NumLib::LocalToGlobalIndexMap const& dof_table) {
                       return NumLib::getIndices(mesh_item_id, dof_table);
                   });

    auto const local_coupled_xs =
        getCoupledLocalSolutions(cpl_xs->coupled_xs, indices_of_processes);
    assert(local_coupled_xs.size() == displacement_size + phasefield_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

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
        auto pressure_ip = _process_data.pressure;
        auto u_corrected = pressure_ip * u;
        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
        {
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;
        }

        elastic_energy += _ip_data[ip].elastic_energy * w;

        // For AT2
        if (_process_data.at_param == 2)
        {
            surface_energy += 0.5 * gc *
                              ((1 - d_ip) * (1 - d_ip) / ls +
                               (dNdx * d).dot((dNdx * d)) * ls) *
                              w;
        }
        // For AT1
        else
        {
            surface_energy +=
                0.375 * gc *
                ((1 - d_ip) / ls + (dNdx * d).dot((dNdx * d)) * ls) * w;
        }

        if (_process_data.crack_pressure)
        {
            pressure_work +=
                pressure_ip * (N_u * u_corrected).dot(dNdx * d) * w;
        }
    }
}
}  // namespace PhaseField
}  // namespace ProcessLib
