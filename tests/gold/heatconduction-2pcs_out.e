CDF      
      len_name      	time_step          num_dim       	num_nodes          num_elem      
num_el_blk        num_el_in_blk1        num_nod_per_el1       num_el_in_blk2     
   num_nod_per_el2       num_elem_var      num_info   �   len_line   Q         api_version       A   version       A   floating_point_word_size            	file_size               maximum_name_length              int64_status             title          tests/heatconduction-2pcs_out.e          
time_whole                            O@   	eb_status                             x   eb_prop1               name      ID              �   coordx                            �   coordy                            �   coordz                            �   eb_names                   
_FillValue                        	�   
coor_names                     
_FillValue                        �   node_num_map                    �      �   connect1                  	elem_type         EDGE2         �         connect2         	         	elem_type         EDGE2         P      �   elem_num_map                    x      �   name_elem_var         
             
_FillValue                        p   vals_elem_var1eb1                          �      OH   vals_elem_var2eb1                          �      O�   vals_elem_var1eb2                          P      P�   vals_elem_var2eb2                          P      P�   elem_var_tab         
                    p   info_records                      <�      �                                                                                                                                                                                             ?�������?ə�����?�333334?ٙ�����?�      ?�333333?�ffffff?陙����?�������?�������                                                                                                                                                                                                                                                                        ?�������?�������?�333334?�������?�      ?�333333?�ffffff?ə�����?�������?�������?љ�����?�333333?�������?�fffffg?�     ?ٙ�����?�333335?�������?�fffffi?�     ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      0                                                                                                                                                                                                                                                               1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          	   
                                                                                                                   	   	   
   
                                                                                                                                                       	   
                                                            T_analytical                                                                                                                                                                                                                                                    T_cell                                                                                                                                                                                                                                                                      ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                    ./phoenix-opt -i tests/heatconduction-2pcs.i### Version Info ###                Framework Information:                                                           MOOSE Version:           git commit d5122356b7 on 2022-09-20                     LibMesh Version:                                                                 PETSc Version:           3.16.6                                                  SLEPc Version:           3.16.2                                                  Current Time:            Sat Oct 15 20:50:26 2022                                Executable Timestamp:    Sat Oct 15 18:48:52 2022                                                                                                                                                                                                  ### Input File ###                                                                                                                                                []                                                                                 inactive                       = (no_default)                                    initial_from_file_timestep     = LATEST                                          initial_from_file_var          = INVALID                                         allow_negative_qweights        = 1                                               custom_blocks                  = (no_default)                                    custom_orders                  = (no_default)                                    element_order                  = AUTO                                            order                          = AUTO                                            side_order                     = AUTO                                            type                           = GAUSS                                         []                                                                                                                                                                [Components]                                                                       inactive                       = (no_default)                                                                                                                     [./comp-1]                                                                         inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = TestComponent                                     control_tags                 = Components                                        enable                       = 1                                                 length                       = 0.5                                               n_elems                      = 20                                                orientation                  = '0 0 1'                                           position                     = '0 0 0'                                         [../]                                                                                                                                                             [./comp-2]                                                                         inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = TestComponent                                     control_tags                 = Components                                        enable                       = 1                                                 length                       = 1                                                 n_elems                      = 10                                                orientation                  = '1 0 0'                                           position                     = '0 0 1'                                         [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      auto_preconditioning           = 1                                               inactive                       = (no_default)                                    isObjectAction                 = 1                                               type                           = Transient                                       abort_on_solve_fail            = 0                                               accept_on_max_fixed_point_iteration = 0                                          accept_on_max_picard_iteration = 0                                               auto_advance                   = INVALID                                         automatic_scaling              = INVALID                                         check_aux                      = 0                                               compute_initial_residual_before_preset_bcs = 0                                   compute_scaling_once           = 1                                               contact_line_search_allowed_lambda_cuts = 2                                      contact_line_search_ltol       = INVALID                                         control_tags                   = (no_default)                                    custom_abs_tol                 = 1e-50                                           custom_pp                      = INVALID                                         custom_rel_tol                 = 1e-08                                           direct_pp_value                = 0                                               disable_fixed_point_residual_norm_check = 0                                      disable_picard_residual_norm_check = 0                                           dt                             = 1                                               dtmax                          = 1e+30                                           dtmin                          = 1e-13                                           enable                         = 1                                               end_time                       = 10                                              error_on_dtmin                 = 1                                               fixed_point_abs_tol            = 1e-50                                           fixed_point_algorithm          = picard                                          fixed_point_force_norms        = 0                                               fixed_point_max_its            = 1                                               fixed_point_min_its            = 1                                               fixed_point_rel_tol            = 1e-08                                           ignore_variables_for_autoscaling = INVALID                                       l_abs_tol                      = 1e-50                                           l_max_its                      = 10000                                           l_tol                          = 1                                               line_search                    = default                                         line_search_package            = petsc                                           max_xfem_update                = 4294967295                                      mffd_type                      = wp                                              n_max_nonlinear_pingpong       = 100                                             n_startup_steps                = 0                                               nl_abs_div_tol                 = 1e+50                                           nl_abs_step_tol                = 0                                               nl_abs_tol                     = 1e-50                                           nl_div_tol                     = 1e+10                                           nl_forced_its                  = 0                                               nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 0                                               nl_rel_tol                     = 1                                               normalize_solution_diff_norm_by_dt = 1                                           num_grids                      = 1                                               num_steps                      = 4294967295                                      off_diagonals_in_auto_scaling  = 0                                               outputs                        = INVALID                                         petsc_options                  = INVALID                                         petsc_options_iname            = INVALID                                         petsc_options_value            = INVALID                                         picard_abs_tol                 = 1e-50                                           picard_custom_pp               = INVALID                                         picard_force_norms             = 0                                               picard_max_its                 = 1                                               picard_rel_tol                 = 1e-08                                           relaxation_factor              = 1                                               relaxed_variables              = (no_default)                                    reset_dt                       = 0                                               resid_vs_jac_scaling_param     = 0                                               residual_and_jacobian_together = 0                                               restart_file_base              = (no_default)                                    scaling_group_variables        = INVALID                                         scheme                         = implicit-euler                                  skip_exception_check           = 0                                               snesmf_reuse_base              = 1                                               solve_type                     = INVALID                                         splitting                      = INVALID                                         ss_check_tol                   = 1e-08                                           ss_tmin                        = 0                                               start_time                     = 0                                               steady_state_detection         = 0                                               steady_state_start_time        = 0                                               steady_state_tolerance         = 1e-08                                           time_period_ends               = INVALID                                         time_period_starts             = INVALID                                         time_periods                   = INVALID                                         timestep_tolerance             = 1e-13                                           trans_ss_check                 = 0                                               transformed_postprocessors     = (no_default)                                    transformed_variables          = (no_default)                                    update_xfem_at_timestep_begin  = 0                                               use_multiapp_dt                = 0                                               verbose                        = 0                                             []                                                                                                                                                                [Outputs]                                                                          append_date                    = 0                                               append_date_format             = INVALID                                         checkpoint                     = 0                                               color                          = 1                                               console                        = 1                                               controls                       = 0                                               csv                            = 0                                               dofmap                         = 0                                               execute_on                     = 'INITIAL TIMESTEP_END'                          exodus                         = 1                                               file_base                      = INVALID                                         gmv                            = 0                                               gnuplot                        = 0                                               hide                           = INVALID                                         inactive                       = (no_default)                                    interval                       = 1                                               json                           = 0                                               nemesis                        = 0                                               output_if_base_contains        = INVALID                                         perf_graph                     = 0                                               perf_graph_live                = 1                                               perf_graph_live_mem_limit      = 100                                             perf_graph_live_time_limit     = 5                                               print_linear_converged_reason  = 1                                               print_linear_residuals         = 1                                               print_mesh_changed_info        = 0                                               print_nonlinear_converged_reason = 1                                             print_nonlinear_residuals      = 1                                               print_perf_log                 = 0                                               show                           = INVALID                                         solution_history               = 0                                               sync_times                     = (no_default)                                    tecplot                        = 0                                               vtk                            = 0                                               xda                            = 0                                               xdr                            = 0                                               xml                            = 0                                             []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ?�      ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?���:�?�3�ի�?��D���9?�V�]q�?�OB�'x�?�Ǎ����?���~Ό?��DXI�X?�hv��ߍ?�,��@�?�,��@�?�hv��ߏ?��DXI�[?���~ΐ?�Ǎ����?�OB�'x�?�V�]q�?��D���;?�3�ի�?���:�?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�QZ�1Bp?ڔ��7.�?�$���?�O��?��� �w�?��� �w�?�O��?�$���?ڔ��7.�?�QZ�1Bp@       ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?�R��j}?�얼	�?؆�}�i?�}��?����Ʒ�?�^>����?�R�Gn�k?�Z婸??�))��B?��+:'?��+:'?�))��B?�Z婸??�R�Gn�k?�^>����?����Ʒ�?�}��?؆�}�i?�얼	�?�R��j}?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�}��?�:Km�?�Q�آ�?�:^��?��Y�m�?��Y�m�?�:^��?�Q�آ�?�:Km�?�}��@      ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?� g���?��+h���?؊��F��?���f�?��&���?�a��<�K?�W�,�?��"�No?�-t?��o��H?��o��H?�-t?��"�No?�W�,�?�a��<�K?��&���?���f�?؊��F��?��+h���?� g���?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�,L4�T+?�E��w?��ؠm�4?치_r^?�ׅ蚗�?�ׅ蚑�?치_ro?��ؠm��?�E�c ?�,L4�l_@      ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?� z��u�?��Hp
�O?؊�uI� ?��Xr{?��:���i?�bHĦ�?�W4�J ?���uЍ?�.�脭?����g'?����g'?�.�脭?���uЍ?�W4�J ?�bHĦ�?��:���i?��Xr|?؊�uI� ?��Hp
�O?� z��u�?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�0.`l�?�K@6��?���./?����\?��jڼZ!?��jڼU�?���Ő?����?�K@6��?�0.`}D@      ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?� {�U�?��I'�>�?؊��%?����)~?��;3��?�bޑ��?�W57B�S?�� �w�?�.�w3�?���E���?���E���?�.�w3�?�� �w�?�W57B�S?�bޑ��?��;3��?����)?؊��%?��I'�>�?� {�U�?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�0`�>��?�K���Da?��u<&7?�g���?�����?=?�����??�g��/?��u<%�?�K���C\?�0`�>�!@      ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?� {��?��I,��/?؊��w�?����sp?��;6��s?�b�F��?�W5;k�?�� ��ff?�.�5�?���J�^�?���J�^�?�.�5�?�� ��ff?�W5;k�?�b�F��?��;6��s?����sp?؊��w�?��I,��/?� {��?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�0i ��?�K�����?��~S4!G?�sz�%�?����d�?����d�?�sz�%�?��~S4!C?�K�����?�0i ��@      ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?� {�RK?��I,�փ?؊���?�����?��;6ι]?�b�^U�?�W5;�]?�� ��d?�.�S)q?���J�-�?���J�-�?�.�S)q?�� ��d?�W5;�]?�b�^U�?��;6ι]?�����?؊���?��I,�ք?� {�RK?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�0ið��?�K�ق�?��/nd�?�t�n�>?���/��?���/��?�t�n�>?��/nd�?�K�ق�?�0ið��@       ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?� {��7?��I,��a?؊�䴄?�����Z?��;6�8�?�b�^�n?�W5;�?�� ��R?�.�S��?���J���?���J���?�.�S��?�� ��R?�W5;�?�b�^�n?��;6�8�?�����Z?؊�䴄?��I,��b?� {��7?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�0i���?�K��U?��C�V?�t�'��?������?������?�t�'��?��C�U?�K��T?�0i���@"      ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?� {��C?��I,���?؊��:?�����?��;6�;�?�b�^�?�W5;�#?�� ���?�.�S�?���J���?���J���?�.�S�?�� ���?�W5;�#?�b�^�?��;6�;�?�����?؊��;?��I,���?� {��B?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�0i�m��?�K��'2?��E�p?�t����?��W{C?��W{C?�t����?��E�p?�K��'3?�0i�m��@$      ?��29�I?�ብ���?�}⦮�c?�N��.�?��GF ��?�U>�=�?�H��A?�k�2�F?�����?��.&`�?��.&`�?�����?�k�2�F?�H��A?�U>�=�?��GF ��?�N��.�?�}⦮�d?�ብ���?��29�R?� {��V?��I,��?؊��R?�����?��;6�;�?�b�^�5?�W5;�=?�� ���?�.�S�?���J���?���J���?�.�S�?�� ���?�W5;�=?�b�^�5?��;6�;�?�����?؊��R?��I,��?� {��V?�g�Su?�.+D� ?栞f;�?� =,l?�$�/�\?�$�/�\?� =,k?栞f;�?�.+D��?�g�Sy?�0iו`z?�K��`��?��E�UY?�t�&?���1 ?���1 ?�t�&?��E�UY?�K��`��?�0iו`z