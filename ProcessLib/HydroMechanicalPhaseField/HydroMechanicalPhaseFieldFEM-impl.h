/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "GeoLib/AnalyticalGeometry.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/Vector3.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace HydroMechanicalPhaseField
{
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
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

    if (process_id == _hydro_process_id)
    {
        assembleWithJacobianForHydroProcessEquations(
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
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        double const t, double const dt,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    assert(local_coupled_solutions.local_coupled_xs.size() ==
           _phasefield_size + _displacement_size + _pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[_phasefield_index],
        _phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_solutions.local_coupled_xs[_displacement_index],
        _displacement_size);
    auto const p = Eigen::Map<
        typename ShapeMatricesType::template VectorType<_pressure_size> const>(
        &local_coupled_solutions.local_coupled_xs[_pressure_index],
        _pressure_size);
    double const p_geo = _process_data.geostatic_pressure;

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

    double const k = _process_data.residual_stiffness(t, x_position)[0];
    double const Kd = _process_data.drained_modulus(t, x_position)[0];
    double const Ks = _process_data.grain_modulus(t, x_position)[0];
    auto const porosity = _process_data.porosity(t, x_position)[0];
    auto rho_sr = _process_data.solid_density(t, x_position)[0];

    double alpha = 0.0;
    if (_process_data.poroelastic_coupling)
        alpha = (1 - Kd / Ks);

    for (int ip = 0; ip < n_integration_points; ip++)
    {
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

        double const p_ip = N.dot(p);
        double const delta_p = p_ip - p_geo;
        double const p_fr =
            (_process_data.fluid_type == FluidType::Fluid_Type::IDEAL_GAS)
                ? p_ip
                : std::numeric_limits<double>::quiet_NaN();
        double const rho_fr =
            _process_data.getFluidDensity(t, x_position, p_fr);
        double const beta_p = _process_data.getFluidCompressibility(p_fr);

        double const rho = rho_sr * (1 - porosity) + porosity * rho_fr;
        double const degradation = d_ip * d_ip * (1 - k) + k;
        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation, _process_data.split_method,
            reg_param, alpha, delta_p);

        auto const& sigma = _ip_data[ip].sigma;
        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        _displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, _displacement_size>::Zero(DisplacementDim,
                                                           _displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, _displacement_size / DisplacementDim>(
                   i, i * _displacement_size / DisplacementDim)
                .noalias() = N;

        auto const& b = _process_data.specific_body_force;

        auto const C_eff = degradation * C_tensile + C_compressive;

        // Check the dimension
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;
        // For debugging purpose
        if (_element.getID() == 973 && ip == 0)
            DBUG("dbg in hydro-mech");
        local_rhs.noalias() -=
            (B.transpose() * sigma - N_u.transpose() * rho * b -
             delta_p * N_u.transpose() * dNdx * d) *
            w;

        local_Jac.noalias() += B.transpose() * C_eff * B * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    assembleWithJacobianForHydroProcessEquations(
        double const t, double const dt, std::vector<double> const& local_xdot,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    assert(local_coupled_solutions.local_coupled_xs.size() ==
           _phasefield_size + _displacement_size + _pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[_phasefield_index],
        _phasefield_size);
    auto const p = Eigen::Map<
        typename ShapeMatricesType::template VectorType<_pressure_size> const>(
        &local_coupled_solutions.local_coupled_xs[_pressure_index],
        _pressure_size);

    auto const p0 = Eigen::Map<
        typename ShapeMatricesType::template VectorType<_pressure_size> const>(
        &local_coupled_solutions.local_coupled_xs0[_pressure_index],
        _pressure_size);

    auto p_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<_pressure_size> const>(
        local_xdot.data(), _pressure_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<_pressure_size,
                                                        _pressure_size>>(
        local_Jac_data, _pressure_size, _pressure_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<_pressure_size>>(
        local_b_data, _pressure_size);

    typename ShapeMatricesType::NodalMatrixType mass =
        ShapeMatricesType::NodalMatrixType::Zero(_pressure_size,
                                                 _pressure_size);

    typename ShapeMatricesType::NodalMatrixType laplace =
        ShapeMatricesType::NodalMatrixType::Zero(_pressure_size,
                                                 _pressure_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    int const n_integration_points = _integration_method.getNumberOfPoints();

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material =
        *_process_data.solid_materials[0];

    auto const bulk_modulus = solid_material.getBulkModulus(t, x_position);
    double const Kd = _process_data.drained_modulus(t, x_position)[0];
    double const Ks = _process_data.grain_modulus(t, x_position)[0];
    double const perm = _process_data.intrinsic_permeability(t, x_position)[0];
    double const mu = _process_data.fluid_viscosity(t, x_position)[0];
    double const ls = _process_data.crack_length_scale(t, x_position)[0];
    double const porosity = _process_data.porosity(t, x_position)[0];
    double width = (*_process_data.width)[_element.getID()];
    double const width_nl_prev =
        (*_process_data.width_nl_prev)[_element.getID()];
    double const width_prev = (*_process_data.width_prev)[_element.getID()];

    double alpha = 0.0;
    if (_process_data.poroelastic_coupling)
        alpha = (1 - Kd / Ks);

    double const fixed_strs1 = _process_data.fixed_strs1;
    double const fixed_strs2 = _process_data.fixed_strs2;

    double ele_d = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;
        ele_d += N.dot(d);
    }
    ele_d = ele_d / n_integration_points;

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const d_ip = N.dot(d);
        double const p0_ip = N.dot(p0);
        double const p_ip = N.dot(p);

        auto& pressure = _ip_data[ip].pressure;
        auto& pressureNL = _ip_data[ip].pressureNL;

        auto& reg_source = _ip_data[ip].reg_source;

        pressure = N.dot(p);

        double const p_fr =
            (_process_data.fluid_type == FluidType::Fluid_Type::IDEAL_GAS)
                ? pressure
                : std::numeric_limits<double>::quiet_NaN();
        double const rho_fr =
            _process_data.getFluidDensity(t, x_position, p_fr);
        double const beta_p = _process_data.getFluidCompressibility(p_fr);
        double m_inv = porosity * beta_p + (alpha - porosity) / Ks;

        auto const vol_strain = Invariants::trace(_ip_data[ip].eps);
        auto const vol_strain_prev = Invariants::trace(_ip_data[ip].eps_prev);
        double const dv_dt = (vol_strain - vol_strain_prev) / dt;
        //        double const dp_dt = (p_ip - p0_ip) / dt;
        double const dp_dt = (pressureNL - p0_ip) / dt;

        // fixed_strs1 = 0.0 -> no pf modified fixed stress
        //               = 1.0 -> pf modified fixed stress

        double grad_d_norm = (dNdx * d).norm();
        if (d_ip < 0.01)
            grad_d_norm = 1 / ls;

        double const modulus_rm = /*(1 - d_ip * d_ip) * fixed_strs1 * m_inv +
                                  d_ip * d_ip * */
            alpha * alpha / Kd;

        mass.noalias() += (m_inv + alpha * alpha / Kd) * N.transpose() * N * w;
        //            ((fixed_strs1 + (1 - fixed_strs1) * d_ip * d_ip) * m_inv +
        //              d_ip * d_ip *  alpha * alpha / Kd) *
        //            N.transpose() * N * w;

        laplace.noalias() += (perm / mu * dNdx.transpose() * dNdx) * w;

        local_rhs.noalias() +=
            (modulus_rm * dp_dt - d_ip * d_ip * alpha * dv_dt) * N * w;

        local_rhs.noalias() += reg_source * N * w;

        double pf_scaling;

        if (_process_data.pf_scaling_method == 0)
        {
            pf_scaling = 3 * (1 - d_ip) / (4 * ls);
        }
        else if (_process_data.pf_scaling_method == 1)
        {
            pf_scaling = grad_d_norm;
        }
        else if (_process_data.pf_scaling_method == 2)
        {
            pf_scaling =
                3 / 8 * ((1 - d_ip) / ls + ls * grad_d_norm * grad_d_norm);
        }
        else if (_process_data.pf_scaling_method == 3)
        {
            pf_scaling = (1 - d_ip) * (1 - d_ip) / ls;
        }

        // For debugging purpose
        if (_element.getID() == 2 && ip == 0)
            DBUG("dbg in hydro");
        if (_element.getID() == 2353 && ip == 0)
            DBUG("dbg in hydro");
        if (d_ip < 1.0)
        {
            double w1, w2;
            double relax1 = fixed_strs1;
            double relax2 = fixed_strs2;

            w1 = relax1 * width + (1 - relax1) * width_nl_prev;
            w2 = relax2 * width + (1 - relax2) * width_nl_prev;

            double const dw_dt = (w1 - width_prev) / dt;

            double const frac_trans = 4 * pow(w2, 3) / (12 * mu);

            auto norm_gamma = (dNdx * d).normalized();
            decltype(dNdx) const dNdx_gamma =
                (dNdx - norm_gamma * norm_gamma.transpose() * dNdx).eval();

            mass.noalias() += width * beta_p / rho_fr *
                              (1 + d_ip * d_ip / m_inv / Kd * fixed_strs2) *
                              pf_scaling * N.transpose() * N * w;

            laplace.noalias() += (frac_trans * dNdx_gamma.transpose() *
                                  dNdx_gamma * pf_scaling) *
                                 w;
            local_rhs.noalias() +=
                -(dw_dt + dp_dt * d_ip * d_ip / m_inv / Kd * width * beta_p /
                              rho_fr * fixed_strs2) *
                pf_scaling * N * w;
        }
    }
    local_Jac.noalias() = laplace + mass / dt;

    local_rhs.noalias() -= laplace * p + mass * p_dot;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    assembleWithJacobianForPhaseFieldEquations(
        double const t, double const dt,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    assert(local_coupled_solutions.local_coupled_xs.size() ==
           _phasefield_size + _displacement_size + _pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[_phasefield_index],
        _phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_solutions.local_coupled_xs[_displacement_index],
        _displacement_size);
    auto const p = Eigen::Map<
        typename ShapeMatricesType::template VectorType<_pressure_size> const>(
        &local_coupled_solutions.local_coupled_xs[_pressure_index],
        _pressure_size);
    double const p_geo = _process_data.geostatic_pressure;

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, _phasefield_size, _phasefield_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, _phasefield_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    double const& reg_param = _process_data.reg_param;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    double ele_d = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;
        ele_d += N.dot(d);
    }
    ele_d = ele_d / n_integration_points;

    double const gc = _process_data.crack_resistance(t, x_position)[0];
    double const ls = _process_data.crack_length_scale(t, x_position)[0];

    double const k = _process_data.residual_stiffness(t, x_position)[0];
    double const Kd = _process_data.drained_modulus(t, x_position)[0];
    double const Ks = _process_data.grain_modulus(t, x_position)[0];
    auto const alpha = (1 - Kd / Ks);

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const p_ip = N.dot(p);

        double const delta_p = (p_ip - p_geo);

        double const degradation = ele_d * ele_d * (1 - k) + k;
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
        _ip_data[ip].updateConstitutiveRelation(
            t, x_position, dt, u, degradation, _process_data.split_method,
            reg_param, alpha, delta_p);

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

        // For debugging purpose
        if (_element.getID() == 973 && ip == 0)
            DBUG("dbg in pf");

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
                 delta_p * dNdx.transpose() * N_u * u) *
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
                 delta_p * dNdx.transpose() * N_u * u) *
                w;
        }
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    computeFractureVelocity(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
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

    auto local_coupled_xs =
        getCoupledLocalSolutions(cpl_xs->coupled_xs, indices_of_processes);
    assert(local_coupled_xs.size() ==
           _phasefield_size + _displacement_size + _pressure_size);

    auto const p = Eigen::Map<PressureVector const>(
        &local_coupled_xs[_pressure_index], _pressure_size);

    int const n_integration_points = _integration_method.getNumberOfPoints();

    GlobalDimVectorType frac_velocity =
        GlobalDimVectorType::Zero(DisplacementDim);
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    double damage = (*_process_data.ele_d)[_element.getID()];
    double width = (*_process_data.width)[_element.getID()];
    double const mu =
        _process_data.fluid_viscosity(_process_data.t, x_position)[0];
    double const perm =
        _process_data.intrinsic_permeability(_process_data.t, x_position)[0];

    if (damage < 0.05)
    {
        for (int ip = 0; ip < n_integration_points; ip++)
        {
            auto const& dNdx = _ip_data[ip].dNdx;
            double trans = (perm + pow(width, 3) / 3) / mu;
            frac_velocity += -trans * dNdx * p;
        }
        frac_velocity = frac_velocity / n_integration_points;
    }

    for (int i = 0; i < DisplacementDim; ++i)
        _process_data.frac_velocity->getComponent(_element.getID(), i) =
            frac_velocity[i];
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    computeFractureNormal(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
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

    auto local_coupled_xs =
        getCoupledLocalSolutions(cpl_xs->coupled_xs, indices_of_processes);
    assert(local_coupled_xs.size() ==
           _phasefield_size + _displacement_size + _pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[_phasefield_index], _phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[_displacement_index], _displacement_size);

    int const n_integration_points = _integration_method.getNumberOfPoints();

    (*_process_data.width_nl_prev)[_element.getID()] =
        (*_process_data.width)[_element.getID()];
    (*_process_data.width)[_element.getID()] = 0.0;

    double ele_d = 0.0;
    double ele_u_dot_grad_d = 0.0;
    GlobalDimVectorType ele_grad_d = GlobalDimVectorType::Zero(DisplacementDim);

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;

        ele_d += N.dot(d);
    }
    ele_d = ele_d / n_integration_points;

    if (ele_d > 0.0 && ele_d < 1.0)
    {
        for (int ip = 0; ip < n_integration_points; ip++)
        {
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            _displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    _displacement_size>::Zero(DisplacementDim,
                                              _displacement_size);

            for (int i = 0; i < DisplacementDim; ++i)
                N_u.template block<1, _displacement_size / DisplacementDim>(
                       i, i * _displacement_size / DisplacementDim)
                    .noalias() = N;

            ele_grad_d += dNdx * d;
            ele_u_dot_grad_d += (N_u * u).dot(dNdx * d);
        }
        ele_grad_d = ele_grad_d / n_integration_points;
        ele_u_dot_grad_d = ele_u_dot_grad_d / n_integration_points;
    }
    else
    {
        ele_grad_d.setZero();
        ele_u_dot_grad_d = 0.0;
    }
    (*_process_data.ele_d)[_element.getID()] = ele_d;
    (*_process_data.ele_u_dot_grad_d)[_element.getID()] = ele_u_dot_grad_d;
    for (int i = 0; i < DisplacementDim; ++i)
        _process_data.ele_grad_d->getComponent(_element.getID(), i) =
            ele_grad_d[i];
}

bool isPointAtCorner(Eigen::Vector3d pnt_end, GeoLib::Point p0,
                     GeoLib::Point p1)
{
    double eps = std::numeric_limits<double>::epsilon();
    if (((abs(pnt_end[0] - p0[0]) < eps) && (abs(pnt_end[1] - p0[1]) < eps) &&
         (abs(pnt_end[2] - p0[2]) < eps)) ||
        ((abs(pnt_end[0] - p1[0]) < eps) && (abs(pnt_end[1] - p1[1]) < eps) &&
         (abs(pnt_end[2] - p1[2]) < eps)))
        return true;
    else
        return false;
}

bool isPointOnEdge(Eigen::Vector3d pnt_end, GeoLib::Point p0, GeoLib::Point p1)
{
    double eps = std::numeric_limits<double>::epsilon();

    // is it on the line?
    if (abs((p0[1] - p1[1]) * (pnt_end[0] - p0[0]) -
            (p0[0] - p1[0]) * (pnt_end[1] - p0[1])) < eps)
    {
        // is it within the range?
        if (((pnt_end[0] >= std::min(p0[0], p1[0])) &&
             (pnt_end[0] <= std::max(p0[0], p1[0]))) &&
            ((pnt_end[1] >= std::min(p0[1], p1[1])) &&
             (pnt_end[1] <= std::max(p0[1], p1[1]))))
            return true;
        else
            return false;
    }
    else
        return false;
}

void findHostElement(MeshLib::Element const& current_ele,
                     Eigen::Vector3d pnt_end,
                     MeshLib::Element const*& neighbor_ele,
                     double const probe_offset)
{
    // first check if the destination point is in current_ele
    int intersection_count = 0;
    GeoLib::Point intersection_point;
    GeoLib::Point seg_start(pnt_end[0], pnt_end[1], pnt_end[2]);
    GeoLib::Point seg_end(pnt_end[0] + probe_offset, pnt_end[1], pnt_end[2]);
    GeoLib::LineSegment probe_line(&seg_start, &seg_end);
    int num_edge = current_ele.getNumberOfEdges();

    for (int i = 0; i < num_edge; i++)
    {
        auto edge_ele = current_ele.getEdge(i);
        auto n0 = *edge_ele->getNode(0);
        auto n1 = *edge_ele->getNode(1);
        GeoLib::Point point_0(n0[0], n0[1], n0[2]);
        GeoLib::Point point_1(n1[0], n1[1], n1[2]);

        // check if pnt_end lies on the corners or the edge
        if (isPointAtCorner(pnt_end, point_0, point_1) ||
            isPointOnEdge(pnt_end, point_0, point_1))
        {
            neighbor_ele = &current_ele;
            return;
        }

        GeoLib::LineSegment seg0(&point_0, &point_1);
        if (GeoLib::lineSegmentIntersect(seg0, probe_line, intersection_point))
            intersection_count++;
    }
    if (intersection_count == 1)
    {
        neighbor_ele = &current_ele;
        return;
    }

    // The point is not in curren_ele. Find the nearest node to perform the
    // neighbor search
    double distance = 0.0;
    int nearest_node_id;
    Eigen::Vector3d node;
    int num_node = current_ele.getNumberOfNodes();
    for (int i = 0; i < num_node; i++)
    {
        node = Eigen::Map<Eigen::Vector3d const>(
            current_ele.getNode(i)->getCoords(), 3);
        if (i == 0)
        {
            distance = (node - pnt_end).norm();
            nearest_node_id = i;
        }
        else if (distance > (node - pnt_end).norm())
        {
            distance = (node - pnt_end).norm();
            nearest_node_id = i;
        }
    }

    // Loop over the neighbor elements that share the nearest node, to find the
    // hosting element for pnt_end

    MeshLib::Node const* nearest_node = current_ele.getNode(nearest_node_id);
    int num_search_ele = nearest_node->getNumberOfElements();
    for (int i = 0; i < num_search_ele; i++)
    {
        intersection_count = 0;
        auto candidate_ele = nearest_node->getElement(i);
        if (current_ele.getID() == candidate_ele->getID())
            goto next_i;
        num_edge = candidate_ele->getNumberOfEdges();
        for (int j = 0; j < num_edge; j++)
        {
            auto edge_ele = candidate_ele->getEdge(j);
            auto n0 = *edge_ele->getNode(0);
            auto n1 = *edge_ele->getNode(1);
            GeoLib::Point point_0(n0[0], n0[1], n0[2]);
            GeoLib::Point point_1(n1[0], n1[1], n1[2]);
            GeoLib::LineSegment seg0(&point_0, &point_1);

            if (GeoLib::lineSegmentIntersect(seg0, probe_line,
                                             intersection_point))
                intersection_count++;
        }
    next_i:
        if (intersection_count == 1)
        {
            neighbor_ele = candidate_ele;
            return;
        }
    }
}
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    computeFractureWidth(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        double const t,
        CoupledSolutionsForStaggeredScheme const* const /*cpl_xs*/,
        MeshLib::Mesh const& /*mesh*/)
{
    double width = 0.0;
    double cumul_grad_d = 0.0;
    double elem_d = (*_process_data.ele_d)[_element.getID()];
    std::vector<int> elem_list;
    double temporal = -1.0;
    double frac_d = 0.05;

    if (_element.getID() == 477)
        DBUG("dbg in width comp");

    if (0.0 < elem_d && elem_d < 1.0 &&
        _process_data.width_comp_visited[_element.getID()] == false)
    {
        std::vector<std::vector<GlobalIndexType>> indices_of_processes;
        indices_of_processes.reserve(dof_tables.size());
        std::transform(dof_tables.begin(), dof_tables.end(),
                       std::back_inserter(indices_of_processes),
                       [&](NumLib::LocalToGlobalIndexMap const& dof_table) {
                           return NumLib::getIndices(mesh_item_id, dof_table);
                       });

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        double const li_inc =
            _process_data.crack_length_scale(t, x_position)[0] /
            _process_data.li_disc;
        double const probe_offset =
            _process_data.crack_length_scale(t, x_position)[0] *
            _process_data.li_disc;

        double CutOff = _process_data.cum_grad_d_CutOff;

        double deviation = 1.0;
        double cod_start = 0.0, cod_end = 0.0;
        double search_dir = 1.0;

        auto node_ref = _element.getCenterOfGravity();
        auto ref_ele_grad_d_head =
            Eigen::Map<typename ShapeMatricesType::template VectorType<
                DisplacementDim> const>(
                &_process_data.ele_grad_d->getComponent(_element.getID(), 0),
                DisplacementDim);
        Eigen::Vector3d ref_ele_grad_d = Eigen::Vector3d::Zero();
        ref_ele_grad_d.head(DisplacementDim) = ref_ele_grad_d_head;
        auto current_ele_grad_d = ref_ele_grad_d;
        Eigen::Vector3d cumul_ele_grad_d = Eigen::Vector3d::Zero();
        auto current_norm = ref_ele_grad_d.normalized();

        Eigen::Vector3d pnt_start, pnt_end;
        MeshLib::Element const* current_ele;
        MeshLib::Element const* neighbor_ele;
        Eigen::Vector3d delta_l = ref_ele_grad_d.normalized() * li_inc;
        double dist = delta_l.norm();

        // integral in positive direction
        pnt_start = Eigen::Map<Eigen::Vector3d const>(node_ref.getCoords(), 3);
        current_ele = &_element;
        current_ele_grad_d = ref_ele_grad_d;
        int count_i = 0;
        int count_frac_elem = 0;
        elem_list.push_back(current_ele->getID());

        while (elem_d < 1.0 && deviation >= 0.0)
        {
            // find the host element at the end of integral
            pnt_end = pnt_start + delta_l;

            const MathLib::Point3d pnt_end_copy{
                {pnt_end[0], pnt_end[1], pnt_end[2]}};

            neighbor_ele =
                current_ele->findElementInNeighboursWithPoint(pnt_end_copy);
            if (neighbor_ele == nullptr)
            {
                DBUG("neighbor not found");
                break;
            }
            if (current_ele->getID() == neighbor_ele->getID())
                count_i++;
            else
                count_i = 1;

            if (count_i > _process_data.li_disc)
            {
                DBUG("count exceeded");
                break;
            }
            elem_d = (*_process_data.ele_d)[neighbor_ele->getID()];
            if (std::find(elem_list.begin(), elem_list.end(),
                          neighbor_ele->getID()) == elem_list.end() &&
                elem_d < 1.0)
                elem_list.push_back(neighbor_ele->getID());

            // check the normal vector
            auto old_norm = current_norm;
            auto old_ele_grad_d = current_ele_grad_d;
            auto current_ele_grad_d_head =
                Eigen::Map<typename ShapeMatricesType::template VectorType<
                    DisplacementDim> const>(
                    &_process_data.ele_grad_d->getComponent(
                        neighbor_ele->getID(), 0),
                    DisplacementDim);
            current_ele_grad_d.head(DisplacementDim) = current_ele_grad_d_head;
            current_norm = current_ele_grad_d.normalized();
            if (current_ele_grad_d.norm() == 0.0 || elem_d < frac_d)
            {
                current_norm = old_norm;
                count_frac_elem++;
                if (count_i == 1)
                    count_frac_elem++;
                if (count_frac_elem > 10)
                    break;
            }

            // line integral
            cod_start = (*_process_data.ele_u_dot_grad_d)[current_ele->getID()];
            cod_end = (*_process_data.ele_u_dot_grad_d)[neighbor_ele->getID()];
            width += 0.5 * dist * (cod_start + cod_end);

            // for next element search
            current_ele = neighbor_ele;
            pnt_start = pnt_end;
            if (current_norm.dot(old_norm) < 0.0)
            {
                search_dir = -1.0 * search_dir;
                ref_ele_grad_d = -1.0 * ref_ele_grad_d;
            }
            delta_l = search_dir * current_norm * li_inc;
            deviation = (ref_ele_grad_d.normalized()).dot(current_norm);
            //  temporal = std::min(abs(deviation), temporal);

            cumul_ele_grad_d =
                cumul_ele_grad_d +
                0.5 * dist * (old_ele_grad_d + current_ele_grad_d);
            //                cumul_ele_grad_d + 0.5 * dist *
            //                                       (old_ele_grad_d.normalized()
            //                                       +
            //                                        current_ele_grad_d.normalized());
        }

        // integral in negative direction

        pnt_start = Eigen::Map<Eigen::Vector3d const>(node_ref.getCoords(), 3);
        current_ele = &_element;
        elem_d = (*_process_data.ele_d)[_element.getID()];
        current_ele_grad_d = ref_ele_grad_d;
        current_norm = current_ele_grad_d.normalized();
        ref_ele_grad_d = -1.0 * ref_ele_grad_d;
        delta_l = ref_ele_grad_d.normalized() * li_inc;
        deviation = -1.0;
        search_dir = -1.0;

        count_i = 0;
        count_frac_elem = 0;
        while (elem_d < 1.0 && deviation <= 0.0)
        {
            // find the host element at the end of integral
            pnt_end = pnt_start + delta_l;
            const MathLib::Point3d pnt_end_copy{
                {pnt_end[0], pnt_end[1], pnt_end[2]}};

            neighbor_ele =
                current_ele->findElementInNeighboursWithPoint(pnt_end_copy);
            if (neighbor_ele == nullptr)
            {
                DBUG("neighbor not found");
                break;
            }
            if (current_ele->getID() == neighbor_ele->getID())
                count_i++;
            else
                count_i = 1;
            if (count_i > _process_data.li_disc)
            {
                DBUG("count exceeded");
                break;
            }
            elem_d = (*_process_data.ele_d)[neighbor_ele->getID()];
            if (std::find(elem_list.begin(), elem_list.end(),
                          neighbor_ele->getID()) == elem_list.end() &&
                elem_d < 1.0)
                elem_list.push_back(neighbor_ele->getID());

            // check the normal vector
            auto old_norm = current_norm;
            auto old_ele_grad_d = current_ele_grad_d;
            auto current_ele_grad_d_head =
                Eigen::Map<typename ShapeMatricesType::template VectorType<
                    DisplacementDim> const>(
                    &_process_data.ele_grad_d->getComponent(
                        neighbor_ele->getID(), 0),
                    DisplacementDim);
            current_ele_grad_d.head(DisplacementDim) = current_ele_grad_d_head;
            current_norm = current_ele_grad_d.normalized();
            if (current_ele_grad_d.norm() == 0.0 || elem_d < frac_d)
            {
                current_norm = old_norm;
                if (count_i == 1)
                    count_frac_elem++;
                if (count_frac_elem > 10)
                    break;
            }

            // line integral
            cod_start = (*_process_data.ele_u_dot_grad_d)[current_ele->getID()];
            cod_end = (*_process_data.ele_u_dot_grad_d)[neighbor_ele->getID()];
            width += 0.5 * dist * (cod_start + cod_end);

            // for next element search
            current_ele = neighbor_ele;
            pnt_start = pnt_end;
            if (current_norm.dot(old_norm) < 0.0)
            {
                search_dir = -1.0 * search_dir;
                ref_ele_grad_d = -1.0 * ref_ele_grad_d;
            }

            delta_l = search_dir * current_norm * li_inc;
            deviation = (ref_ele_grad_d.normalized()).dot(current_norm);
            temporal = std::max(deviation, temporal);

            cumul_ele_grad_d =
                cumul_ele_grad_d +
                0.5 * dist * (old_ele_grad_d + current_ele_grad_d);
        }

        cumul_grad_d = cumul_ele_grad_d.norm();

        if (width < 0.0 || cumul_grad_d > CutOff ||
            /* temporal > -0.4 || */ count_frac_elem > 10 ||
            count_frac_elem == 0)
        {
            /* if ((*_process_data.width_nl_prev)[_element.getID()] > 0)
                 width = width * 0.5;
             else */
            width = 0.0;
        }

        if (count_frac_elem <= 10 && width >= 0.0)
        {
            for (std::size_t i = 0; i < elem_list.size(); i++)
            {
                if ((*_process_data.ele_d)[elem_list[i]] < 0.05)
                {
                    if (elem_list[i] == 2388)
                        DBUG("dbg in width frac width");
                    //                    if
                    //                    ((*_process_data.width)[elem_list[i]]
                    //                    < 1.e-10)
                    //                        (*_process_data.width)[elem_list[i]]
                    //                        = width;
                    //                    else if (width <
                    //                    (*_process_data.width)[elem_list[i]])
                    //                        (*_process_data.width)[elem_list[i]]
                    //                        = width;
                    //                    (*_process_data.cum_grad_d)[elem_list[i]]
                    //                    = temporal;
                }
                //                _process_data.width_comp_visited[elem_list[i]]
            }
        }

        (*_process_data.width)[_element.getID()] = width;
        (*_process_data.cum_grad_d)[_element.getID()] =
            cumul_grad_d;  // temporal;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    approximateFractureWidth(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        double const t,
        CoupledSolutionsForStaggeredScheme const* const /*cpl_xs*/,
        MeshLib::Mesh const& /*mesh*/)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    double const ls = _process_data.crack_length_scale(t, x_position)[0];
    int const n_integration_points = _integration_method.getNumberOfPoints();
    // For debugging purpose
    if (_element.getID() == 2555)
        DBUG("dbg in width approx");
    double width = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const vol_strain = Invariants::trace(_ip_data[ip].eps);

        width = vol_strain * ls / _process_data.fperm_param1;
    }
    if (width < 0)
        (*_process_data.width)[_element.getID()] = 0.0;
    else if ((*_process_data.ele_d)[_element.getID()] < 0.05)
        (*_process_data.width)[_element.getID()] = width / n_integration_points;
}
}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib
