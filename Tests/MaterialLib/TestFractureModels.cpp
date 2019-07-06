/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <limits>

#include <gtest/gtest.h>

#include <Eigen/Eigen>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/FractureModels/LinearElasticIsotropic.h"
#include "MaterialLib/FractureModels/MohrCoulomb.h"

#include "ParameterLib/ConstantParameter.h"

using namespace MaterialLib::Fracture;

static const double eps_sigma = 1e6*1e-5;
static const double eps_C = 1e10*1e-5;

static NumLib::NewtonRaphsonSolverParameters const nonlinear_solver_parameters{
    1000, 1e-16};

TEST(MaterialLib_Fracture, LinearElasticIsotropic)
{
    ParameterLib::ConstantParameter<double> const kn("", 1e11);
    ParameterLib::ConstantParameter<double> const ks("", 1e9);
    LinearElasticIsotropic<2>::MaterialProperties const mp{kn, ks};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    LinearElasticIsotropic<2> fractureModel{penalty_aperture_cutoff,
                                            tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<2>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector2d const w_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma0 = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma_prev = sigma0;
    Eigen::Vector2d w(-1e-5, -1e-5);

    // Result vectors, not initialized.
    Eigen::Vector2d sigma;
    Eigen::Matrix2d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    EXPECT_NEAR(-1e4, sigma[0], eps_sigma);
    EXPECT_NEAR(-1e6, sigma[1], eps_sigma);
    EXPECT_NEAR(1e9, C(0, 0), eps_C);
    EXPECT_NEAR(0, C(0, 1), eps_C);
    EXPECT_NEAR(0, C(1, 0), eps_C);
    EXPECT_NEAR(1e11, C(1, 1), eps_C);

    w << -1e-5, 1e-5;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    EXPECT_NEAR(0, sigma[0], eps_sigma);
    EXPECT_NEAR(0, sigma[1], eps_sigma);
    EXPECT_NEAR(0, C(0, 0), eps_C);
    EXPECT_NEAR(0, C(0, 1), eps_C);
    EXPECT_NEAR(0, C(1, 0), eps_C);
    EXPECT_NEAR(0, C(1, 1), eps_C);
}

TEST(MaterialLib_Fracture, MohrCoulomb2D_elastic)
{
    ParameterLib::ConstantParameter<double> const kn("", 1e11);
    ParameterLib::ConstantParameter<double> const ks("", 1e9);
    ParameterLib::ConstantParameter<double> const phi("", 15);
    ParameterLib::ConstantParameter<double> const psi("", 5);
    ParameterLib::ConstantParameter<double> const c("", 3e16);
    MohrCoulomb::MohrCoulomb<2>::MaterialProperties const mp{kn, ks, phi, psi,
                                                             c};

    double const aperture0 = 1;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb::MohrCoulomb<2> fractureModel{nonlinear_solver_parameters,
                                              penalty_aperture_cutoff,
                                              tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<2>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector2d const w_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma0 = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma_prev = sigma0;
    Eigen::Vector2d w(-1e-5, -1e-5);

    // Result vectors, not initialized.
    Eigen::Vector2d sigma;
    Eigen::Matrix2d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    EXPECT_NEAR(-1e4, sigma[0], eps_sigma);
    EXPECT_NEAR(-1e6, sigma[1], eps_sigma);
    EXPECT_NEAR(1e9, C(0, 0), eps_C);
    EXPECT_NEAR(0, C(0, 1), eps_C);
    EXPECT_NEAR(0, C(1, 0), eps_C);
    EXPECT_NEAR(1e11, C(1, 1), eps_C);

    w << -1e-5, 1e-5;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    EXPECT_NEAR(0, sigma[0], eps_sigma);
    EXPECT_NEAR(0, sigma[1], eps_sigma);
    EXPECT_NEAR(0, C(0, 0), eps_C);
    EXPECT_NEAR(0, C(0, 1), eps_C);
    EXPECT_NEAR(0, C(1, 0), eps_C);
    EXPECT_NEAR(0, C(1, 1), eps_C);
}

TEST(MaterialLib_Fracture, MohrCoulomb2D_negative_t)
{
    ParameterLib::ConstantParameter<double> const kn("", 50e9);
    ParameterLib::ConstantParameter<double> const ks("", 20e9);
    ParameterLib::ConstantParameter<double> const phi("", 15);
    ParameterLib::ConstantParameter<double> const psi("", 5);
    ParameterLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb::MohrCoulomb<2>::MaterialProperties const mp{kn, ks, phi, psi,
                                                             c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb::MohrCoulomb<2> fractureModel{nonlinear_solver_parameters,
                                              penalty_aperture_cutoff,
                                              tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<2>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector2d const w_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma0(-3.46e6, -2e6);
    Eigen::Vector2d const sigma_prev = sigma0;
    Eigen::Vector2d const w(-1.08e-5, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector2d sigma;
    Eigen::Matrix2d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    EXPECT_NEAR(-3.50360e6, sigma[0], eps_sigma);
    EXPECT_NEAR(-2.16271e6, sigma[1], eps_sigma);
    EXPECT_NEAR(1.10723e+09, C(0,0), eps_C);
    EXPECT_NEAR(1.26558e+10, C(0,1), eps_C);
    EXPECT_NEAR(4.13226e+09, C(1,0), eps_C);
    EXPECT_NEAR(4.72319e+10, C(1,1), eps_C);
}


TEST(MaterialLib_Fracture, MohrCoulomb2D_positive_t)
{
    ParameterLib::ConstantParameter<double> const kn("", 50e9);
    ParameterLib::ConstantParameter<double> const ks("", 20e9);
    ParameterLib::ConstantParameter<double> const phi("", 15);
    ParameterLib::ConstantParameter<double> const psi("", 5);
    ParameterLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb::MohrCoulomb<2>::MaterialProperties const mp{kn, ks, phi, psi,
                                                             c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb::MohrCoulomb<2> fractureModel{nonlinear_solver_parameters,
                                              penalty_aperture_cutoff,
                                              tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<2>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector2d const w_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma0(0, -2e6);
    Eigen::Vector2d const sigma_prev = sigma0;
    Eigen::Vector2d const w(20.08e-5, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector2d sigma;
    Eigen::Matrix2d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    ASSERT_NEAR(446608.0, state->getShearYieldFunctionValue(), eps_sigma);
    EXPECT_NEAR(253972.0, sigma[0], eps_sigma);
    EXPECT_NEAR(-2.94784e+06, sigma[1], eps_sigma);
    EXPECT_NEAR(1.10723e+09, C(0, 0), eps_C);
    EXPECT_NEAR(-1.26558e+10, C(0, 1), eps_C);
    EXPECT_NEAR(-4.13226e+09, C(1, 0), eps_C);
    EXPECT_NEAR(4.72319e+10, C(1, 1), eps_C);
}

TEST(MaterialLib_Fracture, MohrCoulomb3D_negative_t1)
{
    ParameterLib::ConstantParameter<double> const kn("", 50e9);
    ParameterLib::ConstantParameter<double> const ks("", 20e9);
    ParameterLib::ConstantParameter<double> const phi("", 15);
    ParameterLib::ConstantParameter<double> const psi("", 5);
    ParameterLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb::MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi,
                                                             c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb::MohrCoulomb<3> fractureModel{nonlinear_solver_parameters,
                                              penalty_aperture_cutoff,
                                              tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma0(-3.46e6, 0, -2e6);
    Eigen::Vector3d const sigma_prev = sigma0;
    Eigen::Vector3d const w(-1.08e-5, 0, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    ASSERT_NEAR(1.06608E+5, state->getShearYieldFunctionValue(), eps_sigma);
    EXPECT_NEAR(-3.50360e6, sigma[0], eps_sigma);
    EXPECT_NEAR(0.0, sigma[1], eps_sigma);
    EXPECT_NEAR(-2.16271e6, sigma[2], eps_sigma);
    EXPECT_NEAR(1.10723e+09, C(0, 0), eps_C);
    EXPECT_NEAR(0.0, C(0, 1), eps_C);
    EXPECT_NEAR(1.26558e+10, C(0, 2), eps_C);
    EXPECT_NEAR(0.0, C(1, 0), eps_C);
    EXPECT_NEAR(20.e9, C(1, 1), eps_C);
    EXPECT_NEAR(0.0, C(1, 2), eps_C);
    EXPECT_NEAR(4.13226e+09, C(2, 0), eps_C);
    EXPECT_NEAR(0.0, C(2, 1), eps_C);
    EXPECT_NEAR(4.72319e+10, C(2, 2), eps_C);
}

TEST(MaterialLib_Fracture, MohrCoulomb3D_positive_t1)
{
    ParameterLib::ConstantParameter<double> const kn("", 50e9);
    ParameterLib::ConstantParameter<double> const ks("", 20e9);
    ParameterLib::ConstantParameter<double> const phi("", 15);
    ParameterLib::ConstantParameter<double> const psi("", 5);
    ParameterLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb::MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi,
                                                             c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb::MohrCoulomb<3> fractureModel{nonlinear_solver_parameters,
                                              penalty_aperture_cutoff,
                                              tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma0(0, 0, -2e6);
    Eigen::Vector3d const sigma_prev = sigma0;
    Eigen::Vector3d const w(20.08e-5, 0, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    ASSERT_NEAR(446608.0, state->getShearYieldFunctionValue(), eps_sigma);
    EXPECT_NEAR(253972.0, sigma[0], eps_sigma);
    EXPECT_NEAR(0.0, sigma[1], eps_sigma);
    EXPECT_NEAR(-2.94784e+06, sigma[2], eps_sigma);
    EXPECT_NEAR(1.10723e+09, C(0, 0), eps_C);
    EXPECT_NEAR(0.0, C(0, 1), eps_C);
    EXPECT_NEAR(-1.26558e+10, C(0, 2), eps_C);
    EXPECT_NEAR(0.0, C(1, 0), eps_C);
    EXPECT_NEAR(20.e9, C(1, 1), eps_C);
    EXPECT_NEAR(0.0, C(1, 2), eps_C);
    EXPECT_NEAR(-4.13226e+09, C(2, 0), eps_C);
    EXPECT_NEAR(0.0, C(2, 1), eps_C);
    EXPECT_NEAR(4.72319e+10, C(2, 2), eps_C);
}


TEST(MaterialLib_Fracture, MohrCoulomb3D_positive_t2)
{
    ParameterLib::ConstantParameter<double> const kn("", 50e9);
    ParameterLib::ConstantParameter<double> const ks("", 20e9);
    ParameterLib::ConstantParameter<double> const phi("", 15);
    ParameterLib::ConstantParameter<double> const psi("", 5);
    ParameterLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb::MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi,
                                                             c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb::MohrCoulomb<3> fractureModel{nonlinear_solver_parameters,
                                              penalty_aperture_cutoff,
                                              tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma0(0, 0, -2e6);
    Eigen::Vector3d const sigma_prev = sigma0;
    Eigen::Vector3d const w(0, 20.08e-5, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    ASSERT_NEAR(446608.0, state->getShearYieldFunctionValue(), eps_sigma);
    EXPECT_NEAR(0.0, sigma[0], eps_sigma);
    EXPECT_NEAR(253972.0, sigma[1], eps_sigma);
    EXPECT_NEAR(-2.94784e+06, sigma[2], eps_sigma);
    EXPECT_NEAR(20.e9, C(0, 0), eps_C);
    EXPECT_NEAR(0.0, C(0, 1), eps_C);
    EXPECT_NEAR(0.0, C(0, 2), eps_C);
    EXPECT_NEAR(0.0, C(1, 0), eps_C);
    EXPECT_NEAR(1.10723e+09, C(1, 1), eps_C);
    EXPECT_NEAR(-1.26558e+10, C(1, 2), eps_C);
    EXPECT_NEAR(0.0, C(2, 0), eps_C);
    EXPECT_NEAR(-4.13226e+09, C(2, 1), eps_C);
    EXPECT_NEAR(4.72319e+10, C(2, 2), eps_C);
}


TEST(MaterialLib_Fracture, MohrCoulomb3D_negative_t1t2)
{
    ParameterLib::ConstantParameter<double> const kn("", 50e9);
    ParameterLib::ConstantParameter<double> const ks("", 20e9);
    ParameterLib::ConstantParameter<double> const phi("", 15);
    ParameterLib::ConstantParameter<double> const psi("", 5);
    ParameterLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb::MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi,
                                                             c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb::MohrCoulomb<3> fractureModel{
        {1000, 1e-9}, penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma0(-3.46e6 / std::sqrt(2), -3.46e6 / std::sqrt(2),
                                 -2e6);
    Eigen::Vector3d const sigma_prev = sigma0;
    Eigen::Vector3d const w(-1.08e-5/std::sqrt(2), -1.08e-5/std::sqrt(2), -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    ASSERT_NEAR(1.06608E+5, state->getShearYieldFunctionValue(), eps_sigma);
    EXPECT_NEAR(-3.50360e6 / std::sqrt(2), sigma[0], eps_sigma);
    EXPECT_NEAR(-3.50360e6 / std::sqrt(2), sigma[1], eps_sigma);
    EXPECT_NEAR(-2.16271e6, sigma[2], eps_sigma);
    EXPECT_NEAR(1.05536e+10, C(0, 0), eps_C);
    EXPECT_NEAR(-9.44638e+09, C(0, 1), eps_C);
    EXPECT_NEAR(8.94897e+09, C(0, 2), eps_C);
    EXPECT_NEAR(-9.44638e+09, C(1, 0), eps_C);
    EXPECT_NEAR(1.05536e+10, C(1, 1), eps_C);
    EXPECT_NEAR(8.94897e+09, C(1, 2), eps_C);
    EXPECT_NEAR(2.92195e+09, C(2, 0), eps_C);
    EXPECT_NEAR(2.92195e+09, C(2, 1), eps_C);
    EXPECT_NEAR(4.72319e+10, C(2, 2), eps_C);
}


TEST(MaterialLib_Fracture, MohrCoulomb3D_negative_t1_positive_t2)
{
    ParameterLib::ConstantParameter<double> const kn("", 50e9);
    ParameterLib::ConstantParameter<double> const ks("", 20e9);
    ParameterLib::ConstantParameter<double> const phi("", 15);
    ParameterLib::ConstantParameter<double> const psi("", 5);
    ParameterLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb::MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi,
                                                             c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;

    MohrCoulomb::MohrCoulomb<3> fractureModel{
        {1000, 1e-9}, penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());
    state->pushBackState();

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma0(-3.46e6 / std::sqrt(2), 3.46e6 / std::sqrt(2),
                                 -2e6);
    Eigen::Vector3d const sigma_prev = sigma0;
    Eigen::Vector3d const w(-1.08e-5/std::sqrt(2), 1.08e-5/std::sqrt(2), -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ParameterLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    ASSERT_NEAR(1.06608E+5, state->getShearYieldFunctionValue(), eps_sigma);
    EXPECT_NEAR(-3.50360e6 / std::sqrt(2), sigma[0], eps_sigma);
    EXPECT_NEAR(3.50360e6 / std::sqrt(2), sigma[1], eps_sigma);
    EXPECT_NEAR(-2.16271e6, sigma[2], eps_sigma);
    EXPECT_NEAR(1.05536e+10, C(0, 0), eps_C);
    EXPECT_NEAR(9.44638e+09, C(0, 1), eps_C);
    EXPECT_NEAR(8.94897e+09, C(0, 2), eps_C);
    EXPECT_NEAR(9.44638e+09, C(1, 0), eps_C);
    EXPECT_NEAR(1.05536e+10, C(1, 1), eps_C);
    EXPECT_NEAR(-8.94897e+09, C(1, 2), eps_C);
    EXPECT_NEAR(2.92195e+09, C(2, 0), eps_C);
    EXPECT_NEAR(-2.92195e+09, C(2, 1), eps_C);
    EXPECT_NEAR(4.72319e+10, C(2, 2), eps_C);
}

