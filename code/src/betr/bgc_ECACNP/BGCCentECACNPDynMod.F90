module BGCCentECACNPDynMod
#include "shr_assert.h"
  !
  ! !DESCRIPTION:
  ! subroutines for stoichiometric configuration of the century bgc
  ! !History, created by Jinyun Tang, Dec, 2014.
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use decompMod          , only : bounds_type
  use clm_varcon         , only : spval
  use clm_varpar         , only : ndecomp_pools
  use ColumnType         , only : col
  use clm_varctl         , only : spinup_state
  implicit none
  save
  private

  public :: calc_cascade_matrix

  logical, public :: ldebug_bgc =.false.

contains

  !-------------------------------------------------------------------------------
  subroutine calc_cascade_matrix(nstvars, nreactions, cn_ratios, cp_ratios, n2_n2o_ratio_denit, pct_sand, &
       centurybgc_vars, nitrogen_limit_flag, phosphos_limit_flag, cascade_matrix)
    !
    ! !DESCRIPTION:
    ! calculate cascade matrix for the decomposition model
    !
    ! !USES:
    use clm_varcon                , only : nitrif_n2o_loss_frac
    use BGCCentECACNPParMod       , only : CNDecompBgcParamsInst, NutrientCompetitionParamsInst
    use MathfuncMod               , only : safe_div
    use BGCCentECACNPDynLibMod    , only : centurybgc_type

    ! !ARGUMENTS:
    integer                       , intent(in)     :: nstvars
    integer                       , intent(in)     :: nreactions
    type(centurybgc_type)         , intent(inout)  :: centurybgc_vars
    real(r8)                      , intent(in)     :: cn_ratios(centurybgc_vars%nom_pools)
    real(r8)                      , intent(in)     :: cp_ratios(centurybgc_vars%nom_pools)
    real(r8)                      , intent(in)     :: n2_n2o_ratio_denit                             !ratio of n2 to n2o during denitrification
    real(r8)                      , intent(in)     :: pct_sand
    real(r8)                      , intent(out)    :: cascade_matrix(nstvars, nreactions)
    logical                       , intent(out) :: nitrogen_limit_flag(centurybgc_vars%nom_pools) !
    logical                       , intent(out) :: phosphos_limit_flag(centurybgc_vars%nom_pools) !
    ! !LOCAL VARIABLES:
    real(r8) :: ftxt, f1, f2
    integer :: k, reac


    associate(                                                             & !
         lit1      => centurybgc_vars%lit1                               , & !
         lit2      => centurybgc_vars%lit2                               , & !
         lit3      => centurybgc_vars%lit3                               , & !
         som1      => centurybgc_vars%som1                               , & !
         som2      => centurybgc_vars%som2                               , & !
         som3      => centurybgc_vars%som3                               , & !
         cwd       => centurybgc_vars%cwd                                , & !
         c_loc     => centurybgc_vars%c_loc                              , & !
         n_loc     => centurybgc_vars%n_loc                              , & !
         p_loc     => centurybgc_vars%p_loc                              , & !
         nelms     => centurybgc_vars%nelms                              , & !
         lid_at_rt => centurybgc_vars%lid_at_rt                          , & !
         lid_o2    => centurybgc_vars%lid_o2                             , & !
         lid_co2   => centurybgc_vars%lid_co2                            , & !
         lid_nh4   => centurybgc_vars%lid_nh4                            , & !
         lid_ch4   => centurybgc_vars%lid_ch4                            , & !
         lid_ar    => centurybgc_vars%lid_ar                             , & !
         lid_no3   => centurybgc_vars%lid_no3                            , & !
         lid_n2o   => centurybgc_vars%lid_n2o                            , & !
         lid_n2    => centurybgc_vars%lid_n2                             , & !
         lid_co2_hr=> centurybgc_vars%lid_co2_hr                         , & !
         lid_n2o_nit=> centurybgc_vars%lid_n2o_nit                       , & !
         lid_minn_nh4_immob => centurybgc_vars%lid_minn_nh4_immob        , & !
         lid_minn_no3_immob => centurybgc_vars%lid_minn_no3_immob        , & !
         lid_minn_nh4_plant => centurybgc_vars%lid_minn_nh4_plant        , & !
         lid_minn_no3_plant => centurybgc_vars%lid_minn_no3_plant        , & !
         lid_nh4_nit        => centurybgc_vars%lid_nh4_nit               , & !
         lid_n2_paere => centurybgc_vars%lid_n2_paere                    , & !
         lid_ch4_paere => centurybgc_vars%lid_ch4_paere                  , & !
         lid_n2o_paere => centurybgc_vars%lid_n2o_paere                  , & !
         lid_o2_paere => centurybgc_vars%lid_o2_paere                    , & !
         lid_ar_paere => centurybgc_vars%lid_ar_paere                    , & !
         lid_co2_paere => centurybgc_vars%lid_co2_paere                  , & !
         lid_minp_solution => centurybgc_vars%lid_minp_solution          , & !
         lid_minp_secondary => centurybgc_vars%lid_minp_secondary        , & !
         lid_minp_occlude =>  centurybgc_vars%lid_minp_occlude           , & !
         lid_plant_minp => centurybgc_vars%lid_plant_minp                , & !
         lid_minp_immob => centurybgc_vars%lid_minp_immob                , & !
         is_aerobic_reac=> centurybgc_vars%is_aerobic_reac               , & !
         is_1st_order => centurybgc_vars%is_1st_order                    , & !
         primvarid    => centurybgc_vars%primvarid                       , & !
         lit1_dek_reac=> centurybgc_vars%lit1_dek_reac                   , & !
         lit2_dek_reac=> centurybgc_vars%lit2_dek_reac                   , & !
         lit3_dek_reac=> centurybgc_vars%lit3_dek_reac                   , & !
         som1_dek_reac=> centurybgc_vars%som1_dek_reac                   , & !
         som2_dek_reac=> centurybgc_vars%som2_dek_reac                   , & !
         som3_dek_reac=> centurybgc_vars%som3_dek_reac                   , & !
         cwd_dek_reac=> centurybgc_vars%cwd_dek_reac                     , & !
         lid_at_rt_reac=> centurybgc_vars%lid_at_rt_reac                 , & !
         lid_no3_den  => centurybgc_vars%lid_no3_den                     , & !
         lid_plant_minn_nh4_up_reac=> centurybgc_vars%lid_plant_minn_nh4_up_reac , & !
         lid_plant_minn_no3_up_reac=> centurybgc_vars%lid_plant_minn_no3_up_reac , & !
         lid_plant_minn_nh4  => centurybgc_vars%lid_plant_minn_nh4       , &
         lid_plant_minn_no3  => centurybgc_vars%lid_plant_minn_no3       , &
         lid_minp_secondary_reac => centurybgc_vars%lid_minp_secondary_reac    , & !
         lid_minp_occlude_reac  => centurybgc_vars%lid_minp_occlude_reac       , & !
         lid_minp_solution_reac => centurybgc_vars%lid_minp_solution_reac      , & !
         lid_plant_minp_up_reac => centurybgc_vars%lid_plant_minp_up_reac, & !

         lid_nh4_nit_reac => centurybgc_vars%lid_nh4_nit_reac            , & !
         lid_no3_den_reac => centurybgc_vars%lid_no3_den_reac            , & !
         lid_n2_aere_reac => centurybgc_vars%lid_n2_aere_reac            , & !
         lid_ch4_aere_reac=> centurybgc_vars%lid_ch4_aere_reac           , & !
         lid_n2o_aere_reac=> centurybgc_vars%lid_n2o_aere_reac           , & !
         lid_o2_aere_reac => centurybgc_vars%lid_o2_aere_reac            , & !
         lid_ar_aere_reac => centurybgc_vars%lid_ar_aere_reac            , & !
         lid_co2_aere_reac=> centurybgc_vars%lid_co2_aere_reac           , & !

         cn_ratios_tgt    => CNDecompBgcParamsInst%cn_ratios_tgt         , & !
         cp_ratios_tgt    => CNDecompBgcParamsInst%cp_ratios_tgt         , & !
         rf_l1s1_bgc      => CNDecompBgcParamsInst%rf_l1s1_bgc           , &
         rf_l2s1_bgc      => CNDecompBgcParamsInst%rf_l2s1_bgc           , & !
         rf_l3s2_bgc      => CNDecompBgcParamsInst%rf_l3s2_bgc           , & !
         rf_s2s1_bgc      => CNDecompBgcParamsInst%rf_s2s1_bgc           , & !
         rf_s3s1_bgc      => CNDecompBgcParamsInst%rf_s3s1_bgc           , & !
         cwd_fcel_bgc     => CNDecompBgcParamsInst%cwd_fcel_bgc          , & !
         cwd_flig_bgc     => CNDecompBgcParamsInst%cwd_flig_bgc            & !

         )

    !initialize all entries to zero
    cascade_matrix = 0._r8
    !higher [nh4] makes lower [no3] competitiveness
    !note all reactions are in the form products - substrates = 0, therefore
    !mass balance is automatically ensured.
    !set up first order reactions
    !---------------------------------------------------------------------------------
    !reaction1, lit1 -> s1
    reac=lit1_dek_reac
    !lit1 + 0.55*o2 -> 0.45 som1 + 0.55co2 + (1/cn_ratios(lit1) - 0.45/cn_ratios_tgt(som1))min_n+ (1/cp_ratios(lit1)-0.45/cp_ratios_tgt(som1))min_p
    cascade_matrix((lit1-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((lit1-1)*nelms+n_loc   ,reac)  = -safe_div(1._r8,cn_ratios(lit1))
    cascade_matrix((lit1-1)*nelms+p_loc   ,reac)  = -safe_div(1._r8,cp_ratios(lit1))

    cascade_matrix(lid_o2                 ,reac)  = -rf_l1s1_bgc
    cascade_matrix((som1-1)*nelms+c_loc   ,reac)  = 1._r8-rf_l1s1_bgc
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)  = safe_div(1._r8-rf_l1s1_bgc,cn_ratios_tgt(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)  = safe_div(1._r8-rf_l1s1_bgc,cp_ratios_tgt(som1))

    cascade_matrix(lid_co2                ,reac)  = rf_l1s1_bgc
    cascade_matrix(lid_nh4                ,reac)  = safe_div(1._r8,cn_ratios(lit1)) - safe_div(1._r8-rf_l1s1_bgc,cn_ratios_tgt(som1))
    cascade_matrix(lid_minp_solution         ,reac)  = safe_div(1._r8,cp_ratios(lit1)) - safe_div(1._r8-rf_l1s1_bgc,cp_ratios_tgt(som1))

    cascade_matrix(lid_minn_nh4_immob     ,reac)  = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_co2_hr             ,reac)  = rf_l1s1_bgc
    cascade_matrix(lid_minp_immob         ,reac)  = -cascade_matrix(lid_minp_solution  ,reac)

    primvarid(reac)       = (lit1-1)*nelms+c_loc
    is_aerobic_reac(reac) = .true.
    nitrogen_limit_flag(reac) = (cascade_matrix(lid_nh4, reac) < 0._r8)
    phosphos_limit_flag(reac) = (cascade_matrix(lid_minp_solution,reac) < 0._r8)
    is_1st_order(reac) = .true.
    !---------------------------------------------------------------------------------
    !reaction 2, lit2 -> s1
    reac = lit2_dek_reac
    !lit2 + 0.5 o2  -> 0.5 som1 + 0.5 co2 + (1/cn_ratios(lit2)-0.5/cn_ratios_tgt(som1))min_n +(1/cp_ratios(lit2)-0.5/cp_ratios_tgt(som1))min_p
    cascade_matrix((lit2-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((lit2-1)*nelms+n_loc   ,reac)   = -safe_div(1._r8,cn_ratios(lit2))
    cascade_matrix((lit2-1)*nelms+p_loc   ,reac)   = -safe_div(1._r8,cp_ratios(lit2))

    cascade_matrix(lid_o2                 ,reac)   = -rf_l2s1_bgc
    cascade_matrix((som1-1)*nelms+c_loc   ,reac)   =  1._r8-rf_l2s1_bgc
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)   =  safe_div(1._r8-rf_l2s1_bgc,cn_ratios_tgt(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)   =  safe_div(1._r8-rf_l2s1_bgc,cp_ratios_tgt(som1))

    cascade_matrix(lid_co2                ,reac)   =  rf_l2s1_bgc
    cascade_matrix(lid_nh4                ,reac)   = safe_div(1._r8,cn_ratios(lit2)) - safe_div(1._r8-rf_l2s1_bgc,cn_ratios_tgt(som1))
    cascade_matrix(lid_minp_solution         ,reac)   = safe_div(1._r8,cn_ratios(lit2)) - safe_div(1._r8-rf_l2s1_bgc,cp_ratios_tgt(som1))

    cascade_matrix(lid_minn_nh4_immob     ,reac)   = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac)   = -cascade_matrix(lid_minp_solution  ,reac)
    cascade_matrix(lid_co2_hr             ,reac)   = rf_l2s1_bgc

    primvarid(reac) = (lit2-1)*nelms+c_loc
    is_aerobic_reac(reac) = .true.
    nitrogen_limit_flag(reac) = (cascade_matrix(lid_nh4, reac) < 0._r8)
    phosphos_limit_flag(reac) = (cascade_matrix(lid_minp_solution,reac) < 0._r8)
    is_1st_order(reac) = .true.
    !---------------------------------------------------------------------------------

    !reaction 3, lit3->s2
    reac = lit3_dek_reac
    !lit3 + 0.5 o2 -> 0.5 som2 + 0.5 co2 + (1/cn_ratios(lit3) - 0.5/cn_ratios(som2))min_n + (1/cp_ratios(lit3)-0.5_r8/cp_ratios(som2))minp
    cascade_matrix((lit3-1)*nelms+c_loc   ,reac) = -1._r8
    cascade_matrix((lit3-1)*nelms+n_loc   ,reac) = -safe_div(1._r8,cn_ratios(lit3))
    cascade_matrix((lit3-1)*nelms+p_loc   ,reac) = -safe_div(1._r8,cp_ratios(lit3))

    cascade_matrix(lid_o2                 ,reac) = -rf_l3s2_bgc
    cascade_matrix((som2-1)*nelms+c_loc   ,reac) =  1._r8-rf_l3s2_bgc
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) =  safe_div(1._r8-rf_l3s2_bgc,cn_ratios_tgt(som2))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) =  safe_div(1._r8-rf_l3s2_bgc,cp_ratios_tgt(som2))

    cascade_matrix(lid_co2                ,reac) = rf_l3s2_bgc
    cascade_matrix(lid_nh4                ,reac) = safe_div(1._r8,cn_ratios(lit3)) - safe_div(1._r8-rf_l3s2_bgc,cn_ratios_tgt(som2))
    cascade_matrix(lid_minp_solution         ,reac) = safe_div(1._r8,cp_ratios(lit3)) - safe_div(1._r8-rf_l3s2_bgc,cp_ratios_tgt(som2))

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_solution  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = rf_l3s2_bgc

    primvarid(reac) = (lit3-1)*nelms+c_loc
    is_aerobic_reac(reac) = .true.
    nitrogen_limit_flag(reac) = (cascade_matrix(lid_nh4, reac) < 0._r8)
    phosphos_limit_flag(reac) = (cascade_matrix(lid_minp_solution,reac) < 0._r8)
    is_1st_order(reac) = .true.

    !---------------------------------------------------------------------------------
    !double check those stoichiometry parameters
    !reaction 4, the partition into som2 and som3 is soil texture dependent
    !SOM1 -> f1*SOM2 + f2*SOm3 + (1-f1-f2)*CO2 + (1/cn_ratios(SOM1)-f1/cn_ratios(SOM2)-f2/cn_ratios(SOM3))*min_n
    ! +(1/cp_ratios(SOM1)-f1/cp_ratios(SOM2)-f2/cp_ratios(SOM3))*min_p
    reac = som1_dek_reac

    ftxt = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - pct_sand)
    f1 = 0.996*(1._r8-ftxt)
    f2 = 0.004*(1._r8-ftxt)
    ftxt = 1._r8-f1-f2
    cascade_matrix((som1-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)  = -safe_div(1._r8,cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)  = -safe_div(1._r8,cp_ratios(som1))

    cascade_matrix(lid_o2                 ,reac) = -ftxt
    cascade_matrix((som3-1)*nelms+c_loc   ,reac)  = f2
    cascade_matrix((som3-1)*nelms+n_loc   ,reac)  = safe_div(f2,cn_ratios_tgt(som3))
    cascade_matrix((som3-1)*nelms+p_loc   ,reac)  = safe_div(f2,cp_ratios_tgt(som3))

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) = safe_div(f1,cn_ratios_tgt(som2))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) = safe_div(f1,cp_ratios_tgt(som2))

    cascade_matrix(lid_co2                ,reac) = ftxt
    cascade_matrix(lid_nh4                ,reac) = safe_div(1._r8,cn_ratios(som1))-safe_div(f1,cn_ratios_tgt(som2))-safe_div(f2,cn_ratios_tgt(som3))
    cascade_matrix(lid_minp_solution         ,reac) = safe_div(1._r8,cp_ratios(som1))-safe_div(f1,cp_ratios_tgt(som2))-safe_div(f2,cp_ratios_tgt(som3))

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_solution  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = ftxt

    primvarid(reac) = (som1-1)*nelms+c_loc
    is_aerobic_reac(reac) = .true.
    nitrogen_limit_flag(reac) = (cascade_matrix(lid_nh4, reac) < 0._r8)
    phosphos_limit_flag(reac) = (cascade_matrix(lid_minp_solution,reac) < 0._r8)
    is_1st_order(reac) = .true.

    !---------------------------------------------------------------------------------
    !reaction 5, som2->som1, som3
    reac = som2_dek_reac
    !som2 + 0.55 o2 -> 0.42 som1 + 0.03som3 + 0.55co2 + (1/cn_ratios(som2)-0.42/cn_ratios(som1)-0.03/cn_ratios(som3)) + (1/cp_raitos(som2)-0.42/cp_ratios(som1)-0.03/cp_ratios(som3))
    cascade_matrix((som2-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((som2-1)*nelms+n_loc   ,reac)   = -safe_div(1._r8,cn_ratios(som2))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac)   = -safe_div(1._r8,cp_ratios(som2))

    cascade_matrix(lid_o2                 ,reac)   = -rf_s2s1_bgc
    cascade_matrix((som1-1)*nelms+c_loc   ,reac)   =  0.93_r8*(1._r8-rf_s2s1_bgc)
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)   =  0.93_r8*safe_div(1._r8-rf_s2s1_bgc,cn_ratios_tgt(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)   =  0.93_r8*safe_div(1._r8-rf_s2s1_bgc,cp_ratios_tgt(som1))

    cascade_matrix((som3-1)*nelms+c_loc   ,reac)   =  0.07_r8*(1._r8-rf_s2s1_bgc)
    cascade_matrix((som3-1)*nelms+n_loc   ,reac)   =  0.07_r8*safe_div(1._r8-rf_s2s1_bgc,cn_ratios_tgt(som3))
    cascade_matrix((som3-1)*nelms+p_loc   ,reac)   =  0.07_r8*safe_div(1._r8-rf_s2s1_bgc,cp_ratios_tgt(som3))

    cascade_matrix(lid_co2                ,reac)   =  rf_s2s1_bgc
    cascade_matrix(lid_nh4                ,reac)   =  safe_div(1._r8,cn_ratios(som2))-0.93_r8*safe_div(1._r8-rf_s2s1_bgc,cn_ratios_tgt(som1)) &
         -0.07_r8*safe_div(1._r8-rf_s2s1_bgc,cn_ratios_tgt(som3))
    cascade_matrix(lid_minp_solution         ,reac)   =  safe_div(1._r8,cp_ratios(som2))-0.93_r8*safe_div(1._r8-rf_s2s1_bgc,cp_ratios_tgt(som1)) &
              -0.07_r8*safe_div(1._r8-rf_s2s1_bgc,cp_ratios_tgt(som3))

    cascade_matrix(lid_minn_nh4_immob     ,reac)   = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac)   = -cascade_matrix(lid_minp_solution  ,reac)
    cascade_matrix(lid_co2_hr             ,reac)   = rf_s2s1_bgc

    primvarid(reac) = (som2-1)*nelms+c_loc
    is_aerobic_reac(reac) = .true.
    nitrogen_limit_flag(reac) = (cascade_matrix(lid_nh4, reac) < 0._r8)
    phosphos_limit_flag(reac) = (cascade_matrix(lid_minp_solution,reac) < 0._r8)
    is_1st_order(reac) = .true.

    !---------------------------------------------------------------------------------
    !reaction 6, s3-> s1
    reac = som3_dek_reac
    !som3 + 0.55 o2 -> 0.45*som1 + 0.55co2 + (1/cn_ratios(som3)-0.45/cn_ratios(som1)) + (1/cp_ratios(som3)-0.45/cp_ratios(som1))
    cascade_matrix((som3-1)*nelms+c_loc   ,reac) = -1._r8
    cascade_matrix((som3-1)*nelms+n_loc   ,reac) = -safe_div(1._r8,cn_ratios(som3))
    cascade_matrix((som3-1)*nelms+p_loc   ,reac) = -safe_div(1._r8,cp_ratios(som3))

    cascade_matrix(lid_o2                 ,reac) = -rf_s3s1_bgc
    cascade_matrix((som1-1)*nelms+c_loc   ,reac) = 1._r8-rf_s3s1_bgc
    cascade_matrix((som1-1)*nelms+n_loc   ,reac) = safe_div(1._r8-rf_s3s1_bgc,cn_ratios_tgt(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac) = safe_div(1._r8-rf_s3s1_bgc,cp_ratios_tgt(som1))

    cascade_matrix(lid_co2                ,reac) = rf_s3s1_bgc
    cascade_matrix(lid_nh4                ,reac) = safe_div(1._r8,cn_ratios(som3)) - safe_div(1._r8-rf_s3s1_bgc,cn_ratios_tgt(som1))
    cascade_matrix(lid_minp_solution         ,reac) = safe_div(1._r8,cp_ratios(som3)) - safe_div(1._r8-rf_s3s1_bgc,cp_ratios_tgt(som1))

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_solution  ,reac)

    cascade_matrix(lid_co2_hr             ,reac) = rf_s3s1_bgc

    primvarid(reac) = (som3-1)*nelms+c_loc
    is_aerobic_reac(reac) = .true.
    nitrogen_limit_flag(reac) = (cascade_matrix(lid_nh4, reac) < 0._r8)
    phosphos_limit_flag(reac) = (cascade_matrix(lid_minp_solution,reac) < 0._r8)
    is_1st_order(reac) = .true.

    !---------------------------------------------------------------------------------
    !reaction 7, the partition into lit1 and lit2 is nutrient dependent, respires co2?
    reac = cwd_dek_reac
    !cwd + o2 -> 0.76lit2 + 0.24*lit3 + (1/cn_ratios(cwd)-0.76/cn_ratios(lit2)-0.24/cn_ratios(lit3)) + (1/cp_ratios(cwd)-0.76/cp_ratios(lit2)-0.24/cp_ratios(lit3))
    cascade_matrix((cwd-1)*nelms+c_loc    ,reac) = -1._r8
    cascade_matrix((cwd-1)*nelms+n_loc    ,reac) = -safe_div(1._r8,cn_ratios(cwd))
    cascade_matrix((cwd-1)*nelms+p_loc    ,reac) = -safe_div(1._r8,cp_ratios(cwd))

    cascade_matrix((lit2-1)*nelms+c_loc   ,reac) = cwd_fcel_bgc
    cascade_matrix((lit2-1)*nelms+n_loc   ,reac) = safe_div(cwd_fcel_bgc,cn_ratios_tgt(lit2))
    cascade_matrix((lit2-1)*nelms+p_loc   ,reac) = safe_div(cwd_fcel_bgc,cp_ratios_tgt(lit2))

    cascade_matrix((lit3-1)*nelms+c_loc   ,reac) = cwd_flig_bgc
    cascade_matrix((lit3-1)*nelms+n_loc   ,reac) = safe_div(cwd_flig_bgc,cn_ratios_tgt(lit3))
    cascade_matrix((lit3-1)*nelms+p_loc   ,reac) = safe_div(cwd_flig_bgc,cp_ratios_tgt(lit3))

    cascade_matrix(lid_nh4                ,reac) = safe_div(1._r8,cn_ratios(cwd)) - safe_div(cwd_fcel_bgc,cn_ratios_tgt(lit2)) - &
         safe_div(cwd_flig_bgc,cn_ratios_tgt(lit3))
    cascade_matrix(lid_minp_solution         ,reac) = safe_div(1._r8,cp_ratios(cwd)) - safe_div(cwd_fcel_bgc,cp_ratios_tgt(lit2)) - &
         safe_div(cwd_flig_bgc,cp_ratios_tgt(lit3))

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_solution  ,reac)

    primvarid(reac) = (cwd-1)*nelms+c_loc
    is_aerobic_reac(reac) = .true.
    nitrogen_limit_flag(reac) = (cascade_matrix(lid_nh4, reac) < 0._r8)
    phosphos_limit_flag(reac) = (cascade_matrix(lid_minp_solution,reac) < 0._r8)
    is_1st_order(reac) = .true.

    !---------------------------------------------------------------------------------
    !reaction 8, nitrification
    reac = lid_nh4_nit_reac
    !NH4(+) + (2-f)O2 + (2-f)OH(-)-> (1-f)NO3(-) + (f/2)N2O + (3-f/2) H2O
    cascade_matrix(lid_nh4 ,reac) = -1._r8
    cascade_matrix(lid_o2  ,reac) = -(2._r8 - nitrif_n2o_loss_frac)
    cascade_matrix(lid_no3 ,reac) = 1._r8 - nitrif_n2o_loss_frac
    cascade_matrix(lid_n2o, reac) = 0.5_r8 * nitrif_n2o_loss_frac

    cascade_matrix(lid_nh4_nit,reac) = 1._r8
    cascade_matrix(lid_n2o_nit,reac) = nitrif_n2o_loss_frac
    primvarid(reac) = lid_nh4
    is_aerobic_reac(reac) = .true.
    is_1st_order(reac) = .true.

    !---------------------------------------------------------------------------------
    !reaction 9, denitrification
    reac = lid_no3_den_reac
    !NO3(-) -> 0.5*f N2O + 0.5* (1-f) N2, where f is a function determined from the century denitrification model
    cascade_matrix(lid_no3 ,reac)    = -1._r8
    cascade_matrix(lid_n2o ,reac)    = 0.5_r8 * 1._r8/(1._r8+n2_n2o_ratio_denit)
    cascade_matrix(lid_n2  ,reac)    = 0.5_r8 * n2_n2o_ratio_denit/(1._r8+n2_n2o_ratio_denit)
    cascade_matrix(lid_no3_den,reac) = 1._r8
    primvarid(reac)                  = lid_no3
    is_1st_order(reac) = .true.

    !---------------------------------------------------------------------------------
    !reaction 10, inorganic P non-equilibrium adsorption
    !P_solution -> p_secondary
    reac = lid_minp_solution_reac
    cascade_matrix(lid_minp_solution,  reac) = -1._r8
    cascade_matrix(lid_minp_secondary, reac) = 1._r8
    primvarid(reac) = lid_minp_solution

    !----------------------------------------------------------------------
    !reaction 11, inorganic P non-equilibrium desorption
    ! p_secondary -> P_solution
    reac = lid_minp_secondary_reac
    cascade_matrix(lid_minp_solution,  reac) = 1._r8
    cascade_matrix(lid_minp_secondary, reac) = -1._r8
    primvarid(reac) = lid_minp_secondary

    !----------------------------------------------------------------------
    !reaction 12, inorganic P occlusion
    !p_secondary -> p_occlude
    reac = lid_minp_occlude_reac
    cascade_matrix(lid_minp_secondary,  reac) = -1._r8
    cascade_matrix(lid_minp_occlude  ,  reac) = 1._r8
    primvarid(reac) = lid_minp_secondary

    !----------------------------------------------------------------------
    !reaction 13, plant mineral nitrogen nh4 uptake
    reac = lid_plant_minn_nh4_up_reac
    !  nh4  -> plant_nitrogen
    cascade_matrix(lid_nh4, reac)        = -1._r8
    cascade_matrix(lid_plant_minn_nh4, reac) = 1._r8
    primvarid(reac) = lid_nh4

    !----------------------------------------------------------------------
    !reaction 14, plant mineral nitrogen no3 uptake
    reac = lid_plant_minn_no3_up_reac
    !  no3  -> plant_nitrogen
    cascade_matrix(lid_no3, reac)        = -1._r8
    cascade_matrix(lid_plant_minn_no3, reac) = 1._r8
    primvarid(reac) = lid_no3

    !----------------------------------------------------------------------
    !reaction 15, plant mineral phosphorus uptake
    reac = lid_plant_minp_up_reac
    ! p_solution -> plant_p
    cascade_matrix(lid_minp_solution, reac) = -1._r8
    cascade_matrix(lid_plant_minp, reac) = 1._r8
    primvarid(reac) = lid_minp_solution

    reac = lid_at_rt_reac
    !ar + o2 -> co2
    cascade_matrix(lid_co2, reac) =  1._r8
    cascade_matrix(lid_o2,  reac) = -1._r8
    is_aerobic_reac(reac) = .true.

    !--------------------------------------------------------------------
    !arenchyma transport
    !second primary variables
    reac                               = lid_o2_aere_reac
    cascade_matrix(lid_o2, reac)       = -1._r8
    cascade_matrix(lid_o2_paere, reac) = 1._r8
    is_aerobic_reac(reac) = .true.

    if ( spinup_state /= 1 ) then
       reac                                = lid_ch4_aere_reac
       cascade_matrix(lid_ch4, reac)       = -1._r8
       cascade_matrix(lid_ch4_paere, reac) = 1._r8

       reac                                = lid_ar_aere_reac
       cascade_matrix(lid_ar, reac)        = -1._r8
       cascade_matrix(lid_ar_paere, reac)  = 1._r8

       reac                                = lid_co2_aere_reac
       cascade_matrix(lid_co2, reac)       = -1._r8
       cascade_matrix(lid_co2_paere, reac) = 1._r8

       reac                                = lid_n2o_aere_reac
       cascade_matrix(lid_n2o, reac)       = -1._r8
       cascade_matrix(lid_n2o_paere, reac) = 1._r8

       reac                                = lid_n2_aere_reac
       cascade_matrix(lid_n2, reac)        = -1._r8
       cascade_matrix(lid_n2_paere, reac)  = 1._r8
    endif

  end associate
end subroutine calc_cascade_matrix



end module BGCCentECACNPDynMod
