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
           phasefield_size + displacement_size + pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[phasefield_index],
        phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_solutions.local_coupled_xs[displacement_index],
        displacement_size);
    auto const p = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        &local_coupled_solutions.local_coupled_xs[pressure_index],
        pressure_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        displacement_size>>(
        local_Jac_data, displacement_size, displacement_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<displacement_size>>(
        local_b_data, displacement_size);

    double const& reg_param = _process_data.reg_param;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    double ele_d = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N = _ip_data[ip].N;
        ele_d += N.dot(d);
    }
    ele_d = ele_d / n_integration_points;

    double const k = _process_data.residual_stiffness(t, x_position)[0];
    double const Kd = _process_data.drained_modulus(t, x_position)[0];
    double const Ks = _process_data.grain_modulus(t, x_position)[0];
    auto const porosity = _process_data.porosity(t, x_position)[0];
    auto rho_sr = _process_data.solid_density(t, x_position)[0];

    auto const alpha = (1 - Kd / Ks);

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

        double const p_ip = N.dot(p);
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
            reg_param);

        auto const& sigma = _ip_data[ip].sigma;
        auto const& C_tensile = _ip_data[ip].C_tensile;
        auto const& C_compressive = _ip_data[ip].C_compressive;

        typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                        displacement_size>
            N_u = ShapeMatricesType::template MatrixType<
                DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                          displacement_size);

        for (int i = 0; i < DisplacementDim; ++i)
            N_u.template block<1, displacement_size / DisplacementDim>(
                   i, i * displacement_size / DisplacementDim)
                .noalias() = N;

        auto const& b = _process_data.specific_body_force;

        auto const C_eff = degradation * C_tensile + C_compressive;

        // Check the dimension
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        local_rhs.noalias() -=
            (B.transpose() * (sigma - d_ip * alpha * p_ip * identity2) -
             N_u.transpose() * rho * b - p_ip * N_u.transpose() * dNdx * d) *
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
           phasefield_size + displacement_size + pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[phasefield_index],
        phasefield_size);
    auto const p = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        &local_coupled_solutions.local_coupled_xs[pressure_index],
        pressure_size);

    auto const p0 = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        &local_coupled_solutions.local_coupled_xs0[pressure_index], pressure_size);

    auto p_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_xdot.data(), pressure_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<pressure_size,
                                                        pressure_size>>(
        local_Jac_data, pressure_size, pressure_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<pressure_size>>(
        local_b_data, pressure_size);

    typename ShapeMatricesType::NodalMatrixType mass =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType laplace =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    double width = (*_process_data.width)[_element.getID()];
    double width_prev = (*_process_data.width_prev)[_element.getID()];
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    int const n_integration_points = _integration_method.getNumberOfPoints();

    //    double ele_d = 0.0;
    //    double ele_grad_d_norm = 0.0;
    double ele_source = 0.0;
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        //        auto const& N = _ip_data[ip].N;
        //        auto const& dNdx = _ip_data[ip].dNdx;
        //       ele_d += N.dot(d);
        //       ele_grad_d_norm += (dNdx * d).norm();
        ele_source += _ip_data[ip].reg_source;
    }
    //  ele_d = ele_d / n_integration_points;
    //  ele_grad_d_norm = ele_grad_d_norm / n_integration_points;
    ele_source = ele_source / n_integration_points;

    double const Kd = _process_data.drained_modulus(t, x_position)[0];
    double const Ks = _process_data.grain_modulus(t, x_position)[0];
    double const perm = _process_data.intrinsic_permeability(t, x_position)[0];
    double const mu = _process_data.fluid_viscosity(t, x_position)[0];

    auto const porosity = _process_data.porosity(t, x_position)[0];

    double const alpha = (1 - Kd / Ks);

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const d_ip = N.dot(d);
        double const p0_ip = N.dot(p0);

        auto& pressure = _ip_data[ip].pressure;
        auto& pressureNL = _ip_data[ip].pressureNL;
        pressure = N.dot(p);

        double const p_fr =
            (_process_data.fluid_type == FluidType::Fluid_Type::IDEAL_GAS)
                ? pressure * 2.e6
                : std::numeric_limits<double>::quiet_NaN();
        double const rho_fr =
            _process_data.getFluidDensity(t, x_position, p_fr);
        double const beta_p = _process_data.getFluidCompressibility(p_fr);
        double m_inv =
            porosity * beta_p + (alpha - porosity) * (1 - alpha) / Ks;

        auto const vol_strain = Invariants::trace(_ip_data[ip].eps);
        auto const vol_strain_prev = Invariants::trace(_ip_data[ip].eps_prev);
        double const dv_dt = (vol_strain - vol_strain_prev) / dt;
        double const dp_dt = (pressureNL - p0_ip) / dt;
        // pf_fixed_strs = 1.0 -> no pf fixed stress
        //               = 0.0 -> pf_fixed stress
        double const pf_fixed_strs = 0.0;

        double const grad_d_norm = (dNdx * d).norm();

        double const modulus_rm =
            alpha * alpha / Kd * d_ip * d_ip +
            m_inv * (1 - d_ip * d_ip) * (1 - pf_fixed_strs);

        local_rhs.noalias() +=
            (-modulus_rm * dp_dt + d_ip * d_ip * alpha * dv_dt) * N * w;

        mass.noalias() += ((1 + pf_fixed_strs * (d_ip * d_ip - 1)) * m_inv +
                           d_ip * d_ip * alpha * alpha / Kd) *
                          N.transpose() * N * w;

        local_rhs.noalias() += ele_source * grad_d_norm * N * w;

        laplace.noalias() += (perm / mu * dNdx.transpose() * dNdx) * w;

        if (d_ip > 0.0 && d_ip < 0.99)
        {
            double const dw_dt = (width - width_prev) / dt;
            auto norm_gamma = (dNdx * d).normalized();
            decltype(dNdx) const dNdx_gamma =
                (dNdx - norm_gamma * norm_gamma.transpose() * dNdx).eval();

            local_rhs.noalias() -= (dw_dt * grad_d_norm) * N * w;
            double const frac_trans = 4 * pow(width, 3) / (12 * mu);

            laplace.noalias() += (frac_trans * dNdx_gamma.transpose() *
                                  dNdx_gamma * grad_d_norm) *
                                 w;
            mass.noalias() +=
                (width * beta_p / rho_fr * grad_d_norm) * N.transpose() * N * w;
        }

        // For debugging purpose
        if (_element.getID() == 1 && ip == 0)
            DBUG("something");
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
           phasefield_size + displacement_size + pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_solutions.local_coupled_xs[phasefield_index],
        phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_solutions.local_coupled_xs[displacement_index],
        displacement_size);
    auto const p = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        &local_coupled_solutions.local_coupled_xs[pressure_index],
        pressure_size);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, phasefield_size, phasefield_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, phasefield_size);

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
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const p_ip = N.dot(p);

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
            reg_param);

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
                 p_ip * dNdx.transpose() * N_u * u -
                 N.transpose() * alpha * p_ip) *
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
                 p_ip * dNdx.transpose() * N_u * u -
                 N.transpose() * alpha * p_ip) *
                w;
        }
    }
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
           phasefield_size + displacement_size + pressure_size);

    auto const d = Eigen::Map<PhaseFieldVector const>(
        &local_coupled_xs[phasefield_index], phasefield_size);
    auto const u = Eigen::Map<DeformationVector const>(
        &local_coupled_xs[displacement_index], displacement_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();

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
            x_position.setIntegrationPoint(ip);
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);

            for (int i = 0; i < DisplacementDim; ++i)
                N_u.template block<1, displacement_size / DisplacementDim>(
                       i, i * displacement_size / DisplacementDim)
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
    double temporal = 0.0;
    if (0.0 < elem_d && elem_d < 0.99)
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

        /*        if (_element.getID() == 1921)
                    DBUG("something");
        */
        // integral in positive direction
        pnt_start = Eigen::Map<Eigen::Vector3d const>(node_ref.getCoords(), 3);
        current_ele = &_element;
        current_ele_grad_d = ref_ele_grad_d;
        int count_i = 0;
        int count_frac_elem = 0;
        while (elem_d < 0.99 && deviation >= 0.0)
        {
            // find the host element at the end of integral
            pnt_end = pnt_start + delta_l;
            findHostElement(*current_ele, pnt_end, neighbor_ele, probe_offset);
            if (current_ele->getID() == neighbor_ele->getID())
                count_i++;
            else
                count_i = 1;

            if (count_i > _process_data.li_disc)
            {
                DBUG("count exceeded");
                break;
            }
            //            if (mesh_item_id == 2652 && count_i == 1)
            //                DBUG("here");

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
            if (current_ele_grad_d.norm() == 0.0)
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
            elem_d = (*_process_data.ele_d)[neighbor_ele->getID()];

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
        while (elem_d < 0.99 && deviation <= 0.0)
        {
            // find the host element at the end of integral
            pnt_end = pnt_start + delta_l;
            findHostElement(*current_ele, pnt_end, neighbor_ele, probe_offset);
            if (current_ele->getID() == neighbor_ele->getID())
                count_i++;
            else
                count_i = 1;
            if (count_i > _process_data.li_disc)
            {
                DBUG("count exceeded");
                break;
            }

            //            if (mesh_item_id == 2652 && count_i == 1)
            //                DBUG("here");

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
            if (current_ele_grad_d.norm() == 0.0)
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
            elem_d = (*_process_data.ele_d)[neighbor_ele->getID()];

            cumul_ele_grad_d =
                cumul_ele_grad_d +
                0.5 * dist * (old_ele_grad_d + current_ele_grad_d);
            //                cumul_ele_grad_d + 0.5 * dist *
            //                                       (old_ele_grad_d.normalized()
            //                                       +
            //                                        current_ele_grad_d.normalized());
        }
        if (width < 0.0 || cumul_ele_grad_d.norm() > CutOff ||
            count_frac_elem > 10)
            width = 0.0;
        cumul_grad_d = cumul_ele_grad_d.norm();
        temporal = deviation;
    }

    (*_process_data.width)[_element.getID()] = width;
    (*_process_data.cum_grad_d)[_element.getID()] = cumul_grad_d;
}

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib
