AddTest(
    NAME PhaseField_3D_Unconfined_Compression
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS cube_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_cube_1e0_pcs_0_ts_10000_t_1.000000.vtu cube_1e0_pcs_0_ts_10000_t_1.000000.vtu displacement displacement 1e-6 1e-6
    expected_cube_1e0_pcs_0_ts_10000_t_1.000000.vtu cube_1e0_pcs_0_ts_10000_t_1.000000.vtu phasefield phasefield 1e-6 1e-6
   )

AddTest(
    NAME PhaseField_2D_StaticHydraulicFracture
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS 2D_static.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_2D_StaticCrack_pcs_1_ts_1_t_1.000000.vtu 2D_StaticCrack_pcs_1_ts_1_t_1.000000.vtu displacement displacement 1e-15 1e-15
    expected_2D_StaticCrack_pcs_1_ts_1_t_1.000000.vtu 2D_StaticCrack_pcs_1_ts_1_t_1.000000.vtu phasefield phasefield 1e-15 1e-15
   )

AddTest(
    NAME PhaseField_3D_beam
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d_stag_1pcs.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 18
    DIFF_DATA
    expected_beam3_stag1pcsAT2_pcs_1_ts_10_t_1.000000.vtu beam3_stag1pcsAT2_pcs_1_ts_10_t_1.000000.vtu displacement displacement 1e-5 0
    expected_beam3_stag1pcsAT2_pcs_1_ts_10_t_1.000000.vtu beam3_stag1pcsAT2_pcs_1_ts_10_t_1.000000.vtu phasefield phasefield 1e-6 0
   )


AddTest(
    NAME LARGE_PhaseField_Staggered_square_line
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_line_h_400.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_square_line_h_400_pcs_0_ts_10000_t_1.000000.vtu square_line_h_400_pcs_0_ts_10000_t_1.000000.vtu displacement displacement 1e-16 0
    expected_square_line_h_400_pcs_0_ts_10000_t_1.000000.vtu square_line_h_400_pcs_0_ts_10000_t_1.000000.vtu phasefield phasefield 1e-16 0
   )

AddTest(
    NAME LARGE_PhaseField_Staggered_square_shear
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_shear_h_400.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_square_shear_h_400_pcs_0_ts_10000_t_1.000000.vtu square_shear_h_400_pcs_0_ts_10000_t_1.000000.vtu displacement displacement 1e-16 0
    expected_square_shear_h_400_pcs_0_ts_10000_t_1.000000.vtu square_shear_h_400_pcs_0_ts_10000_t_1.000000.vtu phasefield phasefield 1e-16 0
   )

AddTest(
    NAME LARGE_PhaseField_propagation
    PATH PhaseField
    EXECUTABLE ogs
    EXECUTABLE_ARGS 2D_propagating_petsc_snesVI.prj -- -snes_monitor -ksp_type cg  -pc_type bjacobi -ksp_rtol 1.e-8 -ksp_atol 1.e-10  -snes_atol 1.e-8 -snes_rtol 1.e-8 -snes_max_it 200 -snes_type vinewtonssls
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_2D_PropagatingCrack_AT1_ell0p03_pcs_1_ts_23_t_0_230000_0.vtu 2D_PropagatingCrack_AT1_ell0p03_pcs_1_ts_23_t_0_230000_0.vtu displacement displacement 1e-16 0
    expected_2D_PropagatingCrack_AT1_ell0p03_pcs_1_ts_23_t_0_230000_0.vtu 2D_PropagatingCrack_AT1_ell0p03_pcs_1_ts_23_t_0_230000_0.vtu phasefield phasefield 1e-16 0
   )
