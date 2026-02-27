

integer, allocatable    ::  rough_U_bot_inds(:,:), rough_U_top_inds(:,:), &
                            rough_V_bot_inds(:,:), rough_V_top_inds(:,:), &
                            rough_W_bot_inds(:,:), rough_W_top_inds(:,:), &
                            rough_T_bot_inds(:,:), rough_T_top_inds(:,:), &
                             bulk_U_bot_inds(:,:),  bulk_U_top_inds(:,:), &
                             bulk_T_bot_inds(:,:),  bulk_T_top_inds(:,:), &
                              aux_rough_slab(:,:), &
                               all_ijk_U  (:,:,:), &
                               all_ijk_V  (:,:,:), &
                               all_ijk_W  (:,:,:), &
                               all_ijk_T  (:,:,:), &
                               all_ijk_div(:,:,:)
real(dp), allocatable   ::       all_c_U  (:,:), &
                                 all_c_V  (:,:), &
                                 all_c_W  (:,:), &
                                 all_c_T  (:,:), &
                                 all_c_div(:,:), &
                                 yp_full(:), buffer_ibm(:)
integer                 ::         n_ghosts_U  , n_fluids_U  , &
                                   n_ghosts_V  , n_fluids_V  , &
                                   n_ghosts_W  , n_fluids_W  , &
                                   n_ghosts_T  , n_fluids_T  , &
                                   n_ghosts_div, n_fluids_div

type(diezDecomp_ibm_type) :: obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div