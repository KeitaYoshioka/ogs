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

    if (local_coupled_solutions.process_id == _hydro_process_id)
    {
        assembleWithJacobianForHydroProcessEquations(
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
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
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
        double const alpha = _process_data.biot_coefficient(t, x_position)[0];
        double const d_ip = N.dot(d);
        double const degradation = d_ip * d_ip * (1 - k) + k;
        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                degradation);

        auto const& sigma_eff = _ip_data[ip].sigma_eff;
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

        auto rho_s = _process_data.solid_density(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;

        double const p_ip = _ip_data[ip].pressure;
        auto const C_eff = degradation * C_tensile + C_compressive;

        // Check the dimension
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        local_rhs.noalias() -=
            (B.transpose() * (sigma_eff - d_ip * alpha * p_ip * identity2) -
             N_u.transpose() * rho_s * b - p_ip * N_u.transpose() * dNdx * d) *
            w;

        local_Jac.noalias() += B.transpose() * C_eff * B * w;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    assembleWithJacobianForHydroProcessEquations(
        double const t, std::vector<double> const& local_xdot,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    auto const& local_d =
        local_coupled_solutions.local_coupled_xs[_phase_field_process_id];
    auto const& local_p =
        local_coupled_solutions.local_coupled_xs[_hydro_process_id];
    assert(local_p.size() == pressure_size);
    assert(local_d.size() == phasefield_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto p = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_p.data(), pressure_size);

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

    /*    typename ShapeMatricesType::NodalMatrixType fixed_stress =
            ShapeMatricesType::NodalMatrixType::Zero(pressure_size,
       pressure_size);

        typename ShapeMatricesType::NodalMatrixType source =
            ShapeMatricesType::NodalMatrixType::Zero(pressure_size,
       pressure_size);
    */
    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    double width = (*_process_data.width)[_element.getID()];
    double width_prev = (*_process_data.width_prev)[_element.getID()];

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const& eps = _ip_data[ip].eps;
        auto const& eps_prev = _ip_data[ip].eps_prev;

        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const m_inv = 1 / _process_data.biot_modulus(t, x_position)[0];
        auto const kappa = _process_data.drained_modulus(t, x_position)[0];

        using Invariants = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>;

        double const d_ip = N.dot(d);
        auto& pressure = _ip_data[ip].pressure;
        auto& pressure_prev = _ip_data[ip].pressure_prev;

        pressure = N.dot(p);
        double const perm =
            _process_data.intrinsic_permeability(t, x_position)[0];
        double const mu = _process_data.fluid_viscosity(t, x_position)[0];

        if (d_ip > 0.0 && d_ip < 1.0)
        {
            double const dw_dt = (width - width_prev) / dt;
            double const grad_d_norm = (dNdx * d).norm();
            auto norm_gamma = (dNdx * d).normalized();

            decltype(dNdx) const dNdx_gamma =
                (dNdx - norm_gamma * norm_gamma.transpose() * dNdx).eval();

            double const frac_trans = 4 * pow(width, 3) / (12 * mu);
            laplace.noalias() += (frac_trans * dNdx_gamma.transpose() *
                                  dNdx_gamma * grad_d_norm) *
                                 w;
            local_rhs.noalias() -= (dw_dt * grad_d_norm) * N * w;
        }
        double const modulus_rm =
            alpha * alpha / kappa - m_inv * (1 - d_ip * d_ip);
        // TODO
        //       double const source =
        double const dp_dt = (pressure - pressure_prev) / dt;

        double const dv_dt =
            (Invariants::trace(eps) - Invariants::trace(eps_prev)) / dt;

        laplace.noalias() += (perm / mu * dNdx.transpose() * dNdx) * w;

        mass.noalias() += (m_inv + d_ip * d_ip * alpha * alpha / kappa) *
                          N.transpose() * N * w;
        local_rhs.noalias() -=
            (modulus_rm * dp_dt + d_ip * d_ip * alpha * dv_dt) * N * w;
    }
    local_Jac.noalias() = laplace + mass / dt;

    local_rhs.noalias() -= laplace * p + mass * p_dot;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
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

    auto const& local_u =
        local_coupled_solutions.local_coupled_xs[_mechanics_related_process_id];
    auto const& local_d =
        local_coupled_solutions.local_coupled_xs[_phase_field_process_id];
    auto const& local_p =
        local_coupled_solutions.local_coupled_xs[_hydro_process_id];

    assert(local_u.size() == displacement_size);
    assert(local_d.size() == phasefield_size);
    assert(local_p.size() == pressure_size);

    auto const local_matrix_size = local_d.size();
    auto d =
        Eigen::Map<PhaseFieldVector const>(local_d.data(), phasefield_size);
    auto u =
        Eigen::Map<DeformationVector const>(local_u.data(), displacement_size);

    auto p = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_p.data(), pressure_size);

    auto local_Jac = MathLib::createZeroedMatrix<PhaseFieldMatrix>(
        local_Jac_data, local_matrix_size, local_matrix_size);
    auto local_rhs = MathLib::createZeroedVector<PhaseFieldVector>(
        local_b_data, local_matrix_size);

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    double const& dt = _process_data.dt;

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const gc = _process_data.crack_resistance(t, x_position)[0];
        double const ls = _process_data.crack_length_scale(t, x_position)[0];

        double const k = _process_data.residual_stiffness(t, x_position)[0];
        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        double const d_ip = N.dot(d);
        double const p_ip = N.dot(p);

        double const degradation = d_ip * d_ip * (1 - k) + k;
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
        _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u,
                                                degradation);

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
        getCurrentLocalSolutions(*cpl_xs, indices_of_processes);
    assert(local_coupled_xs.size() == 3);

    auto const& local_u = local_coupled_xs[_mechanics_related_process_id];
    auto const& local_d = local_coupled_xs[_phase_field_process_id];

    assert(local_d.size() == phasefield_size);

    auto d = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_d.data(), phasefield_size);

    auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
        displacement_size> const>(local_u.data(), displacement_size);

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();

    double ele_d = 0.0;
    double ele_u_dot_grad_d = 0.0;
    GlobalDimVectorType ele_grad_d = GlobalDimVectorType::Zero(DisplacementDim);

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N = _ip_data[ip].N;

        ele_d += N * d;
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

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<
    ShapeFunction, IntegrationMethod,
    DisplacementDim>::findNeighborElement(MeshLib::Element const& current_ele,
                                          GeoLib::LineSegment& LIntegral,
                                          int neighbor_ele_id)
{
    auto e0 = current_ele.getEdge(0);
    auto e1 = current_ele.getEdge(1);
    auto p0 = e0->getNode(0);
    auto p1 = e0->getNode(1);
    auto p2 = e1->getNode(1);
    GeoLib::LineSegment seg0(
        dynamic_cast<GeoLib::Point*>(const_cast<MeshLib::Node*>(p0)),
        dynamic_cast<GeoLib::Point*>(const_cast<MeshLib::Node*>(p1)));
    GeoLib::LineSegment seg1(
        dynamic_cast<GeoLib::Point*>(const_cast<MeshLib::Node*>(p1)),
        dynamic_cast<GeoLib::Point*>(const_cast<MeshLib::Node*>(p2)));
    GeoLib::LineSegment seg2(
        dynamic_cast<GeoLib::Point*>(const_cast<MeshLib::Node*>(p2)),
        dynamic_cast<GeoLib::Point*>(const_cast<MeshLib::Node*>(p0)));

    // Find neighbor element
    GeoLib::Point intersectionpoint;
    if (GeoLib::lineSegmentIntersect(seg0, LIntegral, intersectionpoint))
        neighbor_ele_id = 0;
    else if (GeoLib::lineSegmentIntersect(seg1, LIntegral, intersectionpoint))
        neighbor_ele_id = 1;
    else if (GeoLib::lineSegmentIntersect(seg2, LIntegral, intersectionpoint))
        neighbor_ele_id = 2;
    else
        neighbor_ele_id = -1;
}

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
void HydroMechanicalPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                             DisplacementDim>::
    computeFractureWidth(std::size_t mesh_item_id,
                         std::vector<std::reference_wrapper<
                             NumLib::LocalToGlobalIndexMap>> const& dof_tables,
                         CoupledSolutionsForStaggeredScheme const* const cpl_xs,
                         MeshLib::Mesh const& mesh)
{
    double width = 0.0;
    double elem_d = (*_process_data.ele_d)[_element.getID()];
    if (0.0 < elem_d && elem_d < 1.0)
    {
        assert(cpl_xs != nullptr);

        std::vector<std::vector<GlobalIndexType>> indices_of_processes;
        indices_of_processes.reserve(dof_tables.size());
        std::transform(dof_tables.begin(), dof_tables.end(),
                       std::back_inserter(indices_of_processes),
                       [&](NumLib::LocalToGlobalIndexMap const& dof_table) {
                           return NumLib::getIndices(mesh_item_id, dof_table);
                       });

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        double const ls =
            _process_data.crack_length_scale(_process_data.t, x_position)[0];

        // Define a line using the normal vector

        std::vector<GeoLib::Point*> pnts;

        auto node = _element.getCenterOfGravity();
        auto pnt_start = Eigen::Map<Eigen::Vector3d const>(node.getCoords(), 3);

        auto ele_grad_d_head =
            Eigen::Map<typename ShapeMatricesType::template VectorType<
                DisplacementDim> const>(
                &_process_data.ele_grad_d->getComponent(_element.getID(), 0),
                DisplacementDim);
        Eigen::Vector3d ele_grad_d = Eigen::Vector3d::Zero();
        ele_grad_d.head(DisplacementDim) = ele_grad_d_head;

        Eigen::Vector3d const pnt_end =
            pnt_start + ele_grad_d.normalized() * ls *
                            (_process_data.at_param == 1 ? 5 : 10);

        // Line segment for line integral
        GeoLib::Point LIntegral_end(pnt_end[0], pnt_end[1], pnt_end[2]);
        GeoLib::Point LIntegral_start(pnt_start[0], pnt_start[1], pnt_start[2]);
        GeoLib::LineSegment LIntegral(&LIntegral_start, &LIntegral_end);

        // Find the neighbor
        int neighbor_ele_id = -1;
        findNeighborElement(_element, LIntegral, neighbor_ele_id);

        INFO("neighbor_ele_id %d ",neighbor_ele_id);
        /*
                pnts.push_back(
                    new GeoLib::Point(pnt_start.x(), pnt_start.y(),
           pnt_start.z())); pnts.push_back( new GeoLib::Point(pnt_end.x(),
           pnt_end.y(), pnt_end.z())); GeoLib::Polyline ply(pnts);

                // Find a list of elements that are interseted by the line
                MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(
                    mesh,
                    std::make_unique<MeshGeoToolsLib::SearchLength>(),
                    MeshGeoToolsLib::SearchAllNodes::Yes);

                MeshGeoToolsLib::BoundaryElementsSearcher
           boundary_element_searcher( mesh, mesh_node_searcher);
                std::vector<MeshLib::Element*> const& intersected_elements(
                    boundary_element_searcher.getBoundaryElements(ply));

                // Perform a line integral along the line

                MathLib::Vector3 pnt0, pnt1;
                //        MeshLib::Node pnt0(0, 0, 0), pnt1(0, 0, 0);

                double dist, u_dot_grad_d_0, u_dot_grad_d_1;

                for (int i = 0; i < intersected_elements.size() - 1; ++i)
                {
                    pnt0 = intersected_elements[i]->getCenterOfGravity();
                    pnt1 = intersected_elements[i + 1]->getCenterOfGravity();
                    dist = sqrt(pow(pnt0[0] - pnt1[0], 2) + pow(pnt0[1] -
           pnt1[1], 2) + pow(pnt0[2] - pnt1[2], 2)); u_dot_grad_d_0 =
                        (*_process_data.ele_grad_d)[intersected_elements[i]->getID()];
                    u_dot_grad_d_1 =
                        (*_process_data
                              .ele_grad_d)[intersected_elements[i +
           1]->getID()]; width += 0.5 * dist * (u_dot_grad_d_0 +
           u_dot_grad_d_1);
                }
        */
        // Repeat the procedure for the negative normal direction
    }

    (*_process_data.width)[_element.getID()] = width;
}

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib
