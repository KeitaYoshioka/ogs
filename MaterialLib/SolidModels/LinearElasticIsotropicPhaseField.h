/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MathLib/KelvinVector.h"

namespace MaterialLib
{
namespace Solids
{
namespace Phasefield
{
// Heaviside step function
inline double Heaviside(double v)
{
    return (v < 0.0) ? 0.0 : 1.0;
}

// Regularized Heaviside step function
inline double Heaviside_reg(double v, double reg_param)
{
    return std::exp(v / reg_param) / (1 + std::exp(v / reg_param));
}

// Macaulay brackets: positive for tensile and negative for compressive
inline double Macaulay_pos(double v)
{
    return v * Heaviside(v);
}
inline double Macaulay_neg(double v)
{
    return v * Heaviside(-1.0 * v);
}

// Regularized Macaulay bracket: positive for tensile and negative for
// compressive
inline double Macaulay_pos_reg(double v, double reg_param)
{
    return v * Heaviside_reg(v, reg_param);
}
inline double Macaulay_neg_reg(double v, double reg_param)
{
    return v * Heaviside_reg(-1.0 * v, reg_param);
}

/** Decompose the stiffness into tensile and compressive part following AMor et
 * al.'s decomposition
 */
template <int DisplacementDim>
std::tuple<MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_real */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_compressive */,
           double /* strain_energy_tensile */, double /* elastic_energy */
           >
calculateDegradedStressAmor(
    double const degradation,
    double const bulk_modulus,
    double const mu,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    // calculation of deviatoric parts
    auto const& P_dev = Invariants::deviatoric_projection;
    KelvinVector const epsd_curr = P_dev * eps;

    // Hydrostatic part for the stress and the tangent.
    double const eps_curr_trace = Invariants::trace(eps);

    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    KelvinMatrix C_compressive = KelvinMatrix::Zero();

    double const strain_energy_tensile = bulk_modulus / 2 *
                                             Macaulay_pos(eps_curr_trace) *
                                             Macaulay_pos(eps_curr_trace) +
                                         mu * epsd_curr.transpose() * epsd_curr;

    double const strain_energy_compressive = bulk_modulus / 2 *
                                             Macaulay_neg(eps_curr_trace) *
                                             Macaulay_neg(eps_curr_trace);

    KelvinVector const sigma_tensile =
        bulk_modulus * Macaulay_pos(eps_curr_trace) * Invariants::identity2 +
        2 * mu * epsd_curr;

    KelvinVector const sigma_compressive =
        bulk_modulus * Macaulay_neg(eps_curr_trace) * Invariants::identity2;

    C_tensile.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * Heaviside(eps_curr_trace));
    C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();

    C_compressive.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * Heaviside(eps_curr_trace));

    double const elastic_energy =
        degradation * strain_energy_tensile + strain_energy_compressive;

    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;

    return std::make_tuple(sigma_real, sigma_tensile, C_tensile, C_compressive,
                           strain_energy_tensile, elastic_energy);
}

template <int DisplacementDim>
std::tuple<MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_real */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_compressive */,
           double /* strain_energy_tensile */, double /* elastic_energy */
           >
calculateDegradedStressAmor_reg(
    double const degradation,
    double const bulk_modulus,
    double const mu,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps, double const reg_param)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    // calculation of deviatoric parts
    auto const& P_dev = Invariants::deviatoric_projection;
    KelvinVector const epsd_curr = P_dev * eps;

    // Hydrostatic part for the stress and the tangent.
    double const eps_curr_trace = Invariants::trace(eps);

    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    KelvinMatrix C_compressive = KelvinMatrix::Zero();

    double const strain_energy_tensile = bulk_modulus / 2 *
                                             Macaulay_pos_reg(eps_curr_trace, reg_param) *
                                             Macaulay_pos_reg(eps_curr_trace, reg_param) +
                                         mu * epsd_curr.transpose() * epsd_curr;

    double const strain_energy_compressive = bulk_modulus / 2 *
                                             Macaulay_neg_reg(eps_curr_trace, reg_param) *
                                             Macaulay_neg_reg(eps_curr_trace, reg_param);

    KelvinVector const sigma_tensile =
        bulk_modulus * Macaulay_pos_reg(eps_curr_trace, reg_param) * Invariants::identity2 +
        2 * mu * epsd_curr;

    KelvinVector const sigma_compressive =
        bulk_modulus * Macaulay_neg_reg(eps_curr_trace, reg_param) * Invariants::identity2;

    C_tensile.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * Heaviside_reg(eps_curr_trace,reg_param));
    C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();

    C_compressive.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * Heaviside_reg(eps_curr_trace,reg_param));

    double const elastic_energy =
        degradation * strain_energy_tensile + strain_energy_compressive;

    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;

    return std::make_tuple(sigma_real, sigma_tensile, C_tensile, C_compressive,
                           strain_energy_tensile, elastic_energy);
}

template <int DisplacementDim>
std::tuple<
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> /* sigma_real */,
    MathLib::KelvinVector::KelvinVectorType<
        DisplacementDim> /* sigma_tensile */,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> /* C_tensile */,
    double /* strain_energy_tensile */, double /* elastic_energy */>
calculateIsotropicDegradedStress(
    double const degradation,
    double const bulk_modulus,
    double const mu,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    // calculation of deviatoric parts
    auto const& P_dev = Invariants::deviatoric_projection;
    KelvinVector const epsd_curr = P_dev * eps;

    // Hydrostatic part for the stress and the tangent.
    double const eps_curr_trace = Invariants::trace(eps);

    double const strain_energy_tensile =
        bulk_modulus / 2 * eps_curr_trace * eps_curr_trace +
        mu * epsd_curr.transpose() * epsd_curr;
    double const elastic_energy = degradation * strain_energy_tensile;
    KelvinVector const sigma_tensile =
        bulk_modulus * eps_curr_trace * Invariants::identity2 +
        2 * mu * epsd_curr;
    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    C_tensile.template topLeftCorner<3, 3>().setConstant(bulk_modulus);
    C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();
    KelvinVector const sigma_real = degradation * sigma_tensile;

    return {sigma_real, sigma_tensile, C_tensile, strain_energy_tensile,
            elastic_energy};
}

}  // namespace Phasefield
}  // namespace Solids
}  // namespace MaterialLib
