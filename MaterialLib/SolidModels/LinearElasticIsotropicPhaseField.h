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
#include <array>
#include <boost/math/special_functions/pow.hpp>
#include "BaseLib/quicksort.h"
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

inline double Heaviside(double const v)
{
    return (v < 0) ? 0.0 : 1.0;
}

// Regularized Heaviside step function
inline double Heaviside_reg(double const v, double const reg_param)
{
    if (v < -reg_param)
        return 0.0;
    else if (-reg_param <= v && v < -reg_param / 2)
        return 2 * ((v / reg_param) + 1) * ((v / reg_param) + 1);
    else if (-reg_param / 2 <= v && v < 0.0)
        return 1 - 2 * (v / reg_param) * (v / reg_param);
    else
        return 1.0;
}

// Macaulay brackets: strain is positive for tensile and negative for
// compressive
inline double Macaulay_tens(double const v)
{
    return v * Heaviside(v);
}
inline double Macaulay_comp(double v)
{
    return v * (1 - Heaviside(v));
}

// Macaulay brackets: strain is positive for tensile and negative for
// compressive
inline double Macaulay_tens_reg(double v, double reg_param)
{
    return v * Heaviside_reg(v, reg_param);
}
inline double Macaulay_comp_reg(double v, double reg_param)
{
    return v * (1 - Heaviside_reg(v, reg_param));
}

inline double evaluateHTensMiehe(
    int const i, int const j,
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
           (Macaulay_tens(principal_strain[i]) -
            Macaulay_tens(principal_strain[j])) /
           (principal_strain[i] - principal_strain[j]);
}

inline double evaluateHCompMiehe(
    int const i, int const j,
    Eigen::Matrix<double, 3, 1> const& principal_strain)
{
    if (i == j)
    {
        return 0.0;
    }
    double num_zero = 1.e-14;

    if (fabs(principal_strain[i] - principal_strain[j]) < num_zero)
    {
        return 2 * (1 - Heaviside(principal_strain[i]));
    }
    return 2 *
           (Macaulay_comp(principal_strain[i]) -
            Macaulay_comp(principal_strain[j])) /
           (principal_strain[i] - principal_strain[j]);
}

inline double evaluateHMasonry(
    int const i, int const j, double lambda, double mu,
    Eigen::Matrix<double, 3, 1> const& principal_strain,
    Eigen::Matrix<double, 3, 1> const& principal_strain_tp,
    std::array<std::array<double, 3>, 3> const& alpha)
{
    if (i == j)
    {
        return 0.0;
    }

    double h_ij = 0.0;
    double Tr_a = principal_strain_tp.sum();
    double denom = fabs(principal_strain[0]);
    double num_zero = 1.e-3;
    double eval = fabs(principal_strain[i] - principal_strain[j]);
    if (denom > std::numeric_limits<double>::epsilon())
        eval = eval / denom;
    else
        eval = 0.0;

    if (eval < num_zero)
    {
        for (int k = 0; k < 3; k++)
        {
            for (int l = 0; l < 3; l++)
            {
                h_ij += lambda * alpha[k][i] * (alpha[l][i] - alpha[l][j]);
                if (k == l)
                    h_ij += 2 * mu * alpha[k][i] * (alpha[l][i] - alpha[l][j]);
            }
        }
    }
    else
    {
        for (int k = 0; k < 3; k++)
        {
            h_ij += (lambda * Tr_a + 2 * mu * principal_strain_tp[k]) *
                    (alpha[k][i] - alpha[k][j]) /
                    (principal_strain[i] - principal_strain[j]);
        }
    }
    return h_ij;
}

/** Decompose the stiffness into tensile and compressive part following Masonry
 * decomposition
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
calculateDegradedStressMasonry(
    double const degradation, double const lambda, double const mu,
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

    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    KelvinMatrix C_compressive = KelvinMatrix::Zero();

    auto eps_tensor = MathLib::KelvinVector::kelvinVectorToTensor(eps);

    Eigen::EigenSolver<decltype(eps_tensor)> eigen_solver(eps_tensor);
    Eigen::Matrix<double, 3, 1> principal_strain =
        eigen_solver.eigenvalues().real();
    Eigen::Matrix<double, 3, 1> principal_strain_tensile,
        principal_strain_compressive;
    std::array<std::array<double, 3>, 3> alpha = {
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    std::array<std::array<double, 3>, 3> beta = {
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};

 /*   Eigen::Matrix<double, 3, 1> principal_strain_temp,
        principal_strain_tensile_temp, principal_strain_compressive_temp;
    std::array<std::array<double, 3>, 3> alpha_temp = {
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    std::array<std::array<double, 3>, 3> beta_temp = {
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
*/
    double nu = lambda / (2 * (lambda + mu));
    KelvinVector sigma_tensile = KelvinVector::Zero();
    KelvinVector sigma_compressive = KelvinVector::Zero();

    double num_zero = 1.e-14;  // std::numeric_limits<double>::epsilon();
    for (int i = 0; i < 3; i++)
    {
        if (fabs(principal_strain[i]) < num_zero)
            principal_strain[i] = 0.0;
    }

//    principal_strain_temp = principal_strain;

    // sort eigenvalues and eigenvectors
    std::array<std::size_t, 3> permutation = {0, 1, 2};
    principal_strain = -1 * principal_strain;
    BaseLib::quicksort(principal_strain.data(), 0, DisplacementDim,
                       permutation.data());
    principal_strain = -1 * principal_strain;

    if (principal_strain[DisplacementDim - 1] >= 0.0)
    {
        principal_strain_tensile = principal_strain;
        alpha = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
        beta = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
        //        INFO("all tens");
    }
    else if (principal_strain[1] + nu * principal_strain[2] >= 0.0)
    {
        principal_strain_tensile = {
            principal_strain[0] + nu * principal_strain[2],
            principal_strain[1] + nu * principal_strain[2], 0.0};
        alpha = {{{1.0, 0.0, nu}, {0.0, 1.0, nu}, {0.0, 0.0, 0.0}}};
        beta = {{{0.0, 0.0, -nu}, {0.0, 0.0, -nu}, {0.0, 0.0, 1.0}}};
        //        INFO("not UT");
    }
    else if ((1 - nu) * principal_strain[0] +
                 nu * (principal_strain[1] + principal_strain[2]) >=
             0.0)
    {
        principal_strain_tensile = {
            principal_strain[0] +
                nu / (1 - nu) *
                    (principal_strain[1] + principal_strain[2]),
            0.0, 0.0};
        alpha = {{{1.0, nu / (1 - nu), nu / (1 - nu)},
                       {0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0}}};
        beta = {{{0.0, -nu / (1 - nu), -nu / (1 - nu)},
                      {0.0, 1.0, 0.0},
                      {0.0, 0.0, 1.0}}};
        //        INFO("Uniaxial tension");
    }
    else
    {
        principal_strain_tensile = {0.0, 0.0, 0.0};
        alpha = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
        beta = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
        //        INFO("all comp");
    }

    if (fabs(principal_strain[0]) < num_zero and
        fabs(principal_strain[1]) < num_zero and
        fabs(principal_strain[2]) < num_zero)
    {
        alpha = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
        beta = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
        //        INFO("zero strain");
    }



    principal_strain_compressive =
        principal_strain - principal_strain_tensile;

    //    for (int i = 0; i < 3; i++)
    //    {
    //        for (int j = 0; j < 3; j++)
    //        {
    //            alpha[permutation[i]][j] = alpha_temp[i][j];
    //            beta[permutation[i]][j] = beta_temp[i][j];
    //        }
    //        principal_strain_tensile[permutation[i]] =
    //            principal_strain_tensile_temp[i];
    //        principal_strain_compressive[permutation[i]] =
    //            principal_strain_compressive_temp[i];
    //    }

//    for (int i = 0; i < 3; i++)
//    {
//        principal_strain_tensile[i] = principal_strain_tensile_temp[i];
//        principal_strain_compressive[i] = principal_strain_compressive_temp[i];
//        principal_strain[i] = principal_strain_temp[i];
//    }
    double const sum_strain_tensile = principal_strain_tensile.sum();
    double const sum_strain_compressive = principal_strain_compressive.sum();

    std::array<KelvinVector, 3> M_kelvin;

    for (int i = 0; i < 3; i++)
    {
        M_kelvin[permutation[i]] = MathLib::KelvinVector::tensorToKelvin<
            DisplacementDim>(
            eigen_solver.eigenvectors().real().col(i).normalized() *
            eigen_solver.eigenvectors().real().col(i).normalized().transpose());
    }

    double const strain_energy_tensile =
        lambda / 2 * boost::math::pow<2>(sum_strain_tensile) +
        mu * principal_strain_tensile.squaredNorm();

    double const strain_energy_compressive =
        lambda / 2 * boost::math::pow<2>(sum_strain_compressive) +
        mu * principal_strain_compressive.squaredNorm();

    for (int i = 0; i < 3; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            sigma_tensile += (lambda * sum_strain_tensile +
                              2 * mu * principal_strain_tensile[k]) *
                             alpha[k][i] * M_kelvin[i];
            sigma_compressive += (lambda * sum_strain_compressive +
                                  2 * mu * principal_strain_compressive[k]) *
                                 beta[k][i] * M_kelvin[i];
        }
    }
    KelvinMatrix C_temp = KelvinMatrix::Zero();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    C_tensile.noalias() += lambda * alpha[l][i] * alpha[k][j] *
                                           M_kelvin[i] *
                                           M_kelvin[j].transpose();
                    C_temp.noalias() = M_kelvin[i] * M_kelvin[j].transpose();

                    C_compressive.noalias() += lambda * beta[l][i] *
                                               beta[k][j] * M_kelvin[i] *
                                               M_kelvin[j].transpose();
                    if (k == l)
                    {
                        C_tensile.noalias() += 2 * mu * alpha[l][i] *
                                               alpha[k][j] * M_kelvin[i] *
                                               M_kelvin[j].transpose();
                        C_compressive.noalias() += 2 * mu * beta[l][i] *
                                                   beta[k][j] * M_kelvin[i] *
                                                   M_kelvin[j].transpose();
                    }
                }
            }
            C_tensile.noalias() +=
                evaluateHMasonry(i, j, lambda, mu, principal_strain,
                                 principal_strain_tensile, alpha) *
                aOdotB<DisplacementDim>(M_kelvin[i], M_kelvin[j]);
            C_compressive.noalias() +=
                evaluateHMasonry(i, j, lambda, mu, principal_strain,
                                 principal_strain_compressive, beta) *
                aOdotB<DisplacementDim>(M_kelvin[i], M_kelvin[j]);
        }
    }

    KelvinMatrix C_eff = KelvinMatrix::Zero();
    C_eff = C_tensile + C_compressive;
    double const elastic_energy =
        degradation * strain_energy_tensile + strain_energy_compressive;

    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;
    /*

        double bulk_modulus = lambda + 2 / 3 * mu;
        auto const& P_dev = Invariants::deviatoric_projection;
        KelvinVector const epsd_curr = P_dev * eps;
        double const eps_curr_trace = Invariants::trace(eps);
        sigma_tensile = bulk_modulus * eps_curr_trace * Invariants::identity2 +
                        2 * mu * epsd_curr;
        C_tensile = KelvinMatrix::Zero();
        C_tensile.template topLeftCorner<3, 3>().setConstant(bulk_modulus);
        C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();
        KelvinVector const sigma_real = degradation * sigma_tensile;
        double const strain_energy_tensile =
            bulk_modulus / 2 * eps_curr_trace * eps_curr_trace +
            mu * epsd_curr.transpose() * epsd_curr;
        double const strain_energy_compressive = 0.0;
        double const elastic_energy =
            degradation * strain_energy_tensile + strain_energy_compressive;
    */
    return std::make_tuple(sigma_real, sigma_tensile, C_tensile, C_compressive,
                           strain_energy_tensile, elastic_energy);
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
    double const degradation, double const lambda, double const mu,
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
    auto const& identity2 = Invariants::identity2;

    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    KelvinMatrix C_compressive = KelvinMatrix::Zero();

    // non-const for Eigen solver.
    auto eps_tensor = MathLib::KelvinVector::kelvinVectorToTensor(eps);

    Eigen::EigenSolver<decltype(eps_tensor)> eigen_solver(eps_tensor);
    Eigen::Matrix<double, 3, 1> const principal_strain =
        eigen_solver.eigenvalues().real();
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
        return lambda / 2 * macaulay_squared(sum_strain) +
               mu * principal_strain.unaryExpr(macaulay_squared).sum();
    };

    auto stress_computation = [&](auto&& macaulay) {
        KelvinVector stress = lambda * macaulay(sum_strain) * identity2;
        for (int i = 0; i < 3; i++)
            stress += 2 * mu * macaulay(principal_strain[i]) * M_kelvin[i];
        return stress;
    };

    auto hs = [&](double const v) {
        return reg_param < 0.0 ? Heaviside(v) : Heaviside_reg(v, reg_param);
    };

    auto m_tens = [&](double const v) {
        return reg_param < 0.0 ? Macaulay_tens(v)
                               : Macaulay_tens_reg(v, reg_param);
    };

    auto m_comp = [&](double const v) {
        return reg_param < 0.0 ? Macaulay_comp(v)
                               : Macaulay_comp_reg(v, reg_param);
    };

    double const strain_energy_tensile = strain_energy_computation(m_tens);

    double const strain_energy_compressive = strain_energy_computation(m_comp);

    KelvinVector const sigma_tensile = stress_computation(m_tens);

    KelvinVector const sigma_compressive = stress_computation(m_comp);

    C_tensile.template topLeftCorner<3, 3>().setConstant(lambda *
                                                         hs(sum_strain));
    for (int i = 0; i < 3; i++)
    {
        C_tensile.noalias() += 2 * mu * hs(principal_strain[i]) * M_kelvin[i] *
                               M_kelvin[i].transpose();
        for (int j = 0; j < 3; j++)
        {
            C_tensile.noalias() +=
                mu * evaluateHTensMiehe(i, j, principal_strain) *
                aOdotB<DisplacementDim>(M_kelvin[i], M_kelvin[j]);
        }
    }

    C_compressive.template topLeftCorner<3, 3>().setConstant(
        lambda * (1 - hs(sum_strain)));
    KelvinMatrix C_temp = KelvinMatrix::Zero();
    for (int i = 0; i < 3; i++)
    {
        C_compressive.noalias() += 2 * mu * (1 - hs(principal_strain[i])) *
                                   M_kelvin[i] * M_kelvin[i].transpose();
        C_temp.noalias() = M_kelvin[i] * M_kelvin[i].transpose();
        for (int j = 0; j < 3; j++)
            C_compressive.noalias() +=
                mu * evaluateHCompMiehe(i, j, principal_strain) *
                aOdotB<DisplacementDim>(M_kelvin[i], M_kelvin[j]);
    }
    KelvinMatrix C_eff = KelvinMatrix::Zero();
    C_eff = C_tensile + C_compressive;
    double const elastic_energy =
        degradation * strain_energy_tensile + strain_energy_compressive;

    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;

    return std::make_tuple(sigma_real, sigma_tensile, C_tensile, C_compressive,
                           strain_energy_tensile, elastic_energy);
}

/** Decompose the stiffness into tensile and compressive part following Amor et
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
    auto const& identity2 = Invariants::identity2;

    // Hydrostatic part for the stress and the tangent.
    double const eps_curr_trace = Invariants::trace(eps);

    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    KelvinMatrix C_compressive = KelvinMatrix::Zero();

    auto strain_energy_computation_vol = [&](auto&& macaulay) {
        auto macaulay_squared = [&macaulay](double x) {
            return boost::math::pow<2>(macaulay(x));
        };
        return bulk_modulus / 2 * macaulay_squared(eps_curr_trace);
    };

    auto stress_computation_vol = [&](auto&& macaulay) {
        return bulk_modulus * macaulay(eps_curr_trace) * identity2;
    };

    auto hs = [&](double const v) {
        return reg_param < 0.0 ? Heaviside(v) : Heaviside_reg(v, reg_param);
    };

    auto m_tens = [&](double const v) {
        return reg_param < 0.0 ? Macaulay_tens(v)
                               : Macaulay_tens_reg(v, reg_param);
    };

    auto m_comp = [&](double const v) {
        return reg_param < 0.0 ? Macaulay_comp(v)
                               : Macaulay_comp_reg(v, reg_param);
    };

    double const strain_energy_tensile = strain_energy_computation_vol(m_tens) +
                                         mu * epsd_curr.transpose() * epsd_curr;

    double const strain_energy_compressive =
        strain_energy_computation_vol(m_comp);

    KelvinVector const sigma_tensile =
        stress_computation_vol(m_tens) + 2 * mu * epsd_curr;

    KelvinVector const sigma_compressive = stress_computation_vol(m_comp);

    C_tensile.template topLeftCorner<3, 3>().setConstant(bulk_modulus *
                                                         hs(eps_curr_trace));
    C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();

    C_compressive.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * (1 - hs(eps_curr_trace)));

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
