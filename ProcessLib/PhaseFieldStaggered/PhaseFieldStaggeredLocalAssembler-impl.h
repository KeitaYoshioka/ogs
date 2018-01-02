/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 19, 2016, 2:28 PM
 */

#pragma once

#include "PhaseFieldStaggeredLocalAssembler.h"

#include "NumLib/Function/Interpolation.h"

#include "ProcessLib/StaggeredCouplingTerm.h"

#include "ProcessLib/PhaseFieldSmallDeformation/PhaseFieldSmallDeformationProcess.h"

namespace ProcessLib {
namespace PhaseFieldStaggered {
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void PhaseFieldStaggeredLocalAssembler<ShapeFunction, IntegrationMethod,
                                       GlobalDim>::
    assembleWithCoupledTerm(double const /*t*/,
                            std::vector<double> const & /*local_x*/,
                            std::vector<double> & /*local_M_data*/,
                            std::vector<double> & /*local_K_data*/,
                            std::vector<double> & /*local_b_data*/,
                            LocalCouplingTerm const & /*coupled_term*/) {
  OGS_FATAL("assembleWithCoupledTerm in PFStaggered is not implemented");
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void PhaseFieldStaggeredLocalAssembler<ShapeFunction, IntegrationMethod,
                                       GlobalDim>::
    assembleWithJacobianAndCoupling(
        double const t, std::vector<double> const &local_x,
        std::vector<double> const & /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double> & /*local_M_data*/,
        std::vector<double> & /*local_K_data*/,
        std::vector<double> &local_b_data, std::vector<double> &local_Jac_data,
        LocalCouplingTerm const &coupled_term) {
  SpatialPosition pos;
  pos.setElementID(_element.getID());

  for (auto const &coupled_process_pair : coupled_term.coupled_processes) {
    if (coupled_process_pair.first ==
        std::type_index(
            typeid(ProcessLib::PhaseFieldSmallDeformation::
                       PhaseFieldSmallDeformationProcess<GlobalDim>))) {
      auto const &pcs =
          static_cast<ProcessLib::PhaseFieldSmallDeformation::
                          PhaseFieldSmallDeformationProcess<GlobalDim> const &>(
              coupled_process_pair.second);

      const auto local_u =
          coupled_term.local_coupled_xs.at(coupled_process_pair.first);
      auto const strain_energy_tensile_ips =
          pcs.getIntStrainEnergyTensile(_element.getID());

      SpatialPosition pos;
      pos.setElementID(_element.getID());

      assembleWithCoupledPhaseFieldStaggered(t, local_x,
                                             strain_energy_tensile_ips, local_u,
                                             local_Jac_data, local_b_data);
    } else {
      OGS_FATAL("This coupled process is not presented for "
                "PhaseFieldStaggered process");
    }
  }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void PhaseFieldStaggeredLocalAssembler<ShapeFunction, IntegrationMethod,
                                       GlobalDim>::
    assembleWithCoupledPhaseFieldStaggered(
        double const t, std::vector<double> const &local_x,
        std::vector<double> const &strain_energy_tensile_ips,
        std::vector<double> const &local_u, std::vector<double> &local_Jac_data,
        std::vector<double> &local_rhs_data) {
  auto const local_matrix_size = local_x.size();
  // This assertion is valid only if all nodal d.o.f. use the same shape
  // matrices.
  assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

  auto local_Jac = MathLib::createZeroedMatrix<NodalMatrixType>(
      local_Jac_data, local_matrix_size, local_matrix_size);

  auto local_rhs =
      MathLib::createZeroedVector<RhsVector>(local_rhs_data, local_matrix_size);

  auto d = MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);
  // Multiply local_matrix_size with displacement dim
  auto u = MathLib::toVector<NodalVectorType>(local_u,
                                              local_matrix_size * GlobalDim);

  unsigned const n_integration_points = _integration_method.getNumberOfPoints();

  SpatialPosition x_position;
  assert(strain_energy_tensile_ips.size() == n_integration_points);

  for (unsigned ip = 0; ip < n_integration_points; ip++) {
    auto const w = _integration_method.getWeightedPoint(ip).getWeight() *
                   _shape_matrices[ip].detJ *
                   _shape_matrices[ip].integralMeasure;
    auto const &N = _shape_matrices[ip].N;
    auto const &dNdx = _shape_matrices[ip].dNdx;
    double const d_ip = N.dot(d);

    double const gc = _process_data.crack_resistance(t, x_position)[0];
    double const ls = _process_data.crack_length_scale(t, x_position)[0];
    double strain_energy_tensile = strain_energy_tensile_ips[ip];

    auto &ip_data = _ip_data[ip];
    ip_data.strain_energy_tensile = strain_energy_tensile;

    // TODO: following 3 variables should be set from input file
    bool use_history_variable = false;
    int coupling_level = 1;
    int const atnum = 2; // atnum = 1 will require Variational-Inequality (VI)
                         // non-linear solver e.g.
    bool const has_CrackPres = true;
    double const pres = 0.0;

    if (use_history_variable) {
      double history_variable = _ip_data[ip].history_variable;
      double history_variable_prev = _ip_data[ip].history_variable_prev;

      if (history_variable < strain_energy_tensile) {
        history_variable = strain_energy_tensile;
      } else {
        history_variable = history_variable_prev;
      }

      strain_energy_tensile = history_variable;
    }

    // int a = u.transpose() * dNdx;
    if (atnum == 2) {
      local_Jac.noalias() +=
          (2 * N.transpose() * N * strain_energy_tensile +
           gc * (N.transpose() * N / ls + dNdx.transpose() * dNdx * ls)) *
          w;

      if (coupling_level == 1) {
        if (has_CrackPres) {
          local_rhs.noalias() -=
              (N.transpose() * d_ip * 2 * strain_energy_tensile +
               gc * ((d_ip - 1) * N.transpose() / ls +
                     dNdx.transpose() * dNdx * ls * d)
               /* + pres * u.transpose() * dNdx*/) *
              w;
        } else {
          local_rhs.noalias() -=
              (N.transpose() * d_ip * 2 * strain_energy_tensile +
               gc * ((d_ip - 1) * N.transpose() / ls +
                     dNdx.transpose() * dNdx * ls * d)) *
              w;
        }

      } else if (coupling_level == 2) {
        double const delta_strain_energy_tensile =
            ip_data.delta_strain_energy_tensile;

        local_rhs.noalias() -=
            (N.transpose() * d_ip * 2 *
                 (strain_energy_tensile + delta_strain_energy_tensile) +
             gc * ((d_ip - 1) * N.transpose() / ls +
                   dNdx.transpose() * dNdx * ls * d)) *
            w;
      }
    } else if (atnum == 1) {
      local_Jac.noalias() += (2 * N.transpose() * N * strain_energy_tensile +
                              gc * (0.75 * dNdx.transpose() * dNdx * ls)) *
                             w;

      local_rhs.noalias() -= (N.transpose() * d_ip * 2 * strain_energy_tensile +
                              gc * (0.375 * N.transpose() / ls +
                                    dNdx.transpose() * dNdx * ls * d)) *
                             w;
    }
  }
}

} // namespace PhaseFieldStaggered
} // namespace ProcessLib
