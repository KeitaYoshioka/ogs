AddTest(
    NAME PhaseField_Staggered_3D_Beam
    PATH PhaseFieldStaggered
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam3d_stag.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    ABSTOL 1e-4 RELTOL 1e-4
    DIFF_DATA
    expected_beam3d_stagAT2_pcs_0_ts_10_t_1.000000.vtu beam3d_stagAT2_pcs_0_ts_10_t_1.000000.vtu displacement displacement
    expected_beam3d_stagAT2_pcs_0_ts_10_t_1.000000.vtu beam3d_stagAT2_pcs_0_ts_10_t_1.000000.vtu phasefield phasefield
   )
