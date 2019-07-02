/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Eigenvalues>
#include <boost/math/special_functions/pow.hpp>
#include "MathLib/KelvinVector.h"

namespace MaterialLib
{
namespace Solids
{
namespace Phasefield
{
template <int DisplacementDim>
MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> aOdotB(
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& A,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& B);

// Heaviside step function
inline double Heaviside(double v)
{
    return (v < 0) ? 0.0 : 1.0;
}

// Regularized Heaviside step function
inline double Heaviside_reg(double v, double reg_param)
{
    return 1 / (1 + std::exp(-(v + 10 * reg_param) / reg_param));
}

// Macaulay brackets: positive for tensile and negative for compressive
inline double Macaulay_pos(double v)
{
    return v * Heaviside(v);
}
inline double Macaulay_neg(double v)
{
    return v * (1 - Heaviside(v));
}

// Regularized Macaulay bracket: positive for tensile and negative for
// compressive
inline double Macaulay_pos_reg(double v, double reg_param)
{
    return v * Heaviside_reg(v, reg_param);
}
inline double Macaulay_neg_reg(double v, double reg_param)
{
    return v * (1 - Heaviside_reg(v, reg_param));
}

inline double evaluateHPos(int const i, int const j,
                           Eigen::Matrix<double, 3, 1> const& principal_strain)
{
    if (i == j)
    {
        return 0.0;
    }
    if (fabs(principal_strain[i] - principal_strain[j]) <
        std::numeric_limits<double>::epsilon())
    {
        return 2 * Heaviside(principal_strain[i]);
    }
    return 2 *
           (Macaulay_pos(principal_strain[i]) -
            Macaulay_pos(principal_strain[j])) /
           (principal_strain[i] - principal_strain[j]);
}

inline double evaluateHNeg(int const i, int const j,
                           Eigen::Matrix<double, 3, 1> const& principal_strain)
{
    if (i == j)
    {
        return 0.0;
    }
    if (fabs(principal_strain[i] - principal_strain[j]) <
        std::numeric_limits<double>::epsilon())
    {
        return 2 * (1 - Heaviside(principal_strain[i]));
    }
    return 2 *
           (Macaulay_neg(principal_strain[i]) -
            Macaulay_neg(principal_strain[j])) /
           (principal_strain[i] - principal_strain[j]);
}

/** Decompose the stiffness into tensile and compressive part following Miehe et
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
calculateDegradedStressMiehe(
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
    auto const& identity2 = Invariants::identity2;

    // Hydrostatic part for the stress and the tangent.
    double const eps_curr_trace = Invariants::trace(eps);

    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    KelvinMatrix C_compressive = KelvinMatrix::Zero();

    // non-const for Eigen solver.
    auto eps_tensor = MathLib::KelvinVector::kelvinVectorToTensor(eps);

    Eigen::EigenSolver<decltype(eps_tensor)> eigen_solver(eps_tensor);
    Eigen::Matrix<double, 3, 1> const principal_strain =
        eigen_solver.eigenvalues().real();
    //    auto const eigen_vectors = eigen_solver.eigenvectors().real();
    //        auto eigen_vectors =
    //        eigen_solver.eigenvectors().real().col().normalized();
    //        eigen_vectors.col(0).normalize();
    double const sum_strain = principal_strain.sum();

    std::array<KelvinVector, 3> M_kelvin;

    for (int i = 0; i < 3; i++)
    {
        M_kelvin[i] = MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
            eigen_solver.eigenvectors().real().col(i).normalized() *
            eigen_solver.eigenvectors().real().col(i).normalized().transpose());
    }

    auto strain_energy_computation = [&](auto&& macaulay) {
        auto macaulay_squared = [&macaulay](double x) {
            return boost::math::pow<2>(macaulay(x));
        };
        return bulk_modulus / 2 * macaulay_squared(sum_strain) +
               mu * principal_strain.unaryExpr(macaulay_squared).sum();
    };

    double const strain_energy_tensile =
        strain_energy_computation(Macaulay_pos);

    double const strain_energy_compressive =
        strain_energy_computation(Macaulay_neg);

    auto stress_computation = [&](auto&& macaulay) {
        KelvinVector stress = bulk_modulus * macaulay(sum_strain) * identity2;
        for (int i = 0; i < DisplacementDim; i++)
            stress += mu * principal_strain[i] * M_kelvin[i];
        return stress;
    };

    KelvinVector const sigma_tensile = stress_computation(Macaulay_pos);

    KelvinVector const sigma_compressive = stress_computation(Macaulay_neg);

    C_tensile.template topLeftCorner<3, 3>().setConstant(bulk_modulus *
                                                         Heaviside(sum_strain));
    for (int i = 0; i < DisplacementDim; i++)
    {
        C_tensile.noalias() += 2 * mu * Heaviside(principal_strain[i]) *
                               M_kelvin[i] * M_kelvin[i].transpose();
        for (int j = 0; j < DisplacementDim; j++)
            C_tensile.noalias() +=
                mu * evaluateHPos(i, j, principal_strain) *
                aOdotB<DisplacementDim>(M_kelvin[i], M_kelvin[j]);
    }

    C_compressive.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * (1 - Heaviside(eps_curr_trace)));
    for (int i = 0; i < DisplacementDim; i++)
    {
        C_compressive.noalias() += 2 * mu *
                                   (1 - Heaviside(principal_strain[i])) *
                                   M_kelvin[i] * M_kelvin[i].transpose();
        for (int j = 0; j < DisplacementDim; j++)
            C_compressive.noalias() +=
                mu * evaluateHNeg(i, j, principal_strain) *
                aOdotB<DisplacementDim>(M_kelvin[i], M_kelvin[j]);
    }

    double const elastic_energy =
        degradation * strain_energy_tensile + strain_energy_compressive;

    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;

    return std::make_tuple(sigma_real, sigma_tensile, C_tensile, C_compressive,
                           strain_energy_tensile, elastic_energy);
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

    double const strain_energy_tensile =
        bulk_modulus / 2 * boost::math::pow<2>(Macaulay_pos(eps_curr_trace)) +
        mu * epsd_curr.transpose() * epsd_curr;

    double const strain_energy_compressive =
        bulk_modulus / 2 * boost::math::pow<2>(Macaulay_neg(eps_curr_trace));

    KelvinVector const sigma_tensile =
        bulk_modulus * Macaulay_pos(eps_curr_trace) * Invariants::identity2 +
        2 * mu * epsd_curr;

    KelvinVector const sigma_compressive =
        bulk_modulus * Macaulay_neg(eps_curr_trace) * Invariants::identity2;

    C_tensile.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * Heaviside(eps_curr_trace));
    C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();

    C_compressive.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * (1 - Heaviside(eps_curr_trace)));

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
    double const degradation, double const bulk_modulus, double const mu,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps,
    double const reg_param)
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

    double const strain_energy_tensile =
        bulk_modulus / 2 *
            boost::math::pow<2>(Macaulay_pos_reg(eps_curr_trace, reg_param)) +
        mu * epsd_curr.transpose() * epsd_curr;

    double const strain_energy_compressive =
        bulk_modulus / 2 *
        boost::math::pow<2>(Macaulay_neg_reg(eps_curr_trace, reg_param));

    KelvinVector const sigma_tensile =
        bulk_modulus * Macaulay_pos_reg(eps_curr_trace, reg_param) *
            Invariants::identity2 +
        2 * mu * epsd_curr;

    KelvinVector const sigma_compressive =
        bulk_modulus * Macaulay_neg_reg(eps_curr_trace, reg_param) *
        Invariants::identity2;

    C_tensile.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * Heaviside_reg(eps_curr_trace, reg_param));
    C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();

    C_compressive.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * (1 - Heaviside_reg(eps_curr_trace, reg_param)));

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
