begin
  start_write_time = get_cpu_time()

  DateYMD = systemfunc("date '+%m%d%Y'")

;---Set some NetCDF options before we open the file
  setfileoption("nc","preFill",False)
  setfileoption("nc","defineMode",True)
;---Open a new NetCDF file to write to
  fout_name="sbgc.cdomcnp_pars."+DateYMD+".nc"
  system("rm -f " + fout_name)
  fout=addfile(fout_name, "c")
;---create dimensions
  allpfts          = ispan(1,1,1)
  allpfts!0        ="allpfts"
  allpfts@long_nmae="allpft index"
  allpfts@units    ="none"

  alllithos          = ispan(1,1,1)
  alllithos!0        ="alllithos"
  alllithos@long_nmae="alllitho index"
  alllithos@units    ="none"

  pfts          = ispan(1,25,1)
  pfts!0        ="pfts"
  pfts@long_nmae="pft index"
  pfts@units    ="none"

  lithos          = ispan(1,16,1)
  lithos!0        ="lithos"
  lithos@long_nmae="lithology index"
  lithos@units    ="none"

  surfsoil          = ispan(1,2,1)
  surfsoil!0        ="surfsoil"
  surfsoil@long_name="at surface or in soil index"
  surfsoil@units    = "none"

  dims=dimsizes(allpfts)
  dim_names=(/allpfts!0/)
  dim_unlimited=(/False/)
  filedimdef(fout,dim_names,dims,dim_unlimited)
  delete(dims)

  dims=dimsizes(surfsoil)
  dim_names=(/surfsoil!0/)
  dim_unlimited=(/False/)
  filedimdef(fout,dim_names,dims,dim_unlimited)
  delete(dims)

  dims=dimsizes(alllithos)
  dim_names=(/alllithos!0/)
  dim_unlimited=(/False/)
  filedimdef(fout,dim_names,dims,dim_unlimited)
  delete(dims)

  dims=dimsizes(lithos)
  dim_names=(/lithos!0/)
  dim_unlimited=(/False/)
  filedimdef(fout,dim_names,dims,dim_unlimited)
  delete(dims)

  dims=dimsizes(pfts)
  dim_names=(/pfts!0/)
  dim_unlimited=(/False/)
  filedimdef(fout,dim_names,dims,dim_unlimited)

;---define varaible
  filevardef(fout,"cwd_fcel","float",(/allpfts!0/))
  filevardef(fout,"cwd_flig","float",(/allpfts!0/))
  filevardef(fout,"fwd_fcel","float",(/allpfts!0/))
  filevardef(fout,"fwd_flig","float",(/allpfts!0/))
  filevardef(fout,"lwd_fcel","float",(/allpfts!0/))
  filevardef(fout,"lwd_flig","float",(/allpfts!0/))
  filevardef(fout,"surface_tension_water","float",(/allpfts!0/))
  filevardef(fout,"rij_kro_a","float",(/allpfts!0/))
  filevardef(fout,"rij_kro_alpha","float",(/allpfts!0/))
  filevardef(fout,"rij_kro_beta","float",(/allpfts!0/))
  filevardef(fout,"rij_kro_gamma","float",(/allpfts!0/))
  filevardef(fout,"rij_kro_delta","float",(/allpfts!0/))
  filevardef(fout,"minpsi_bgc","float",(/allpfts!0/))
  filevardef(fout,"k_nitr_max","float",(/allpfts!0/))
  filevardef(fout,"Q10","float",(/allpfts!0/))
  filevardef(fout,"froz_q10","float",(/allpfts!0/))
  filevardef(fout,"decomp_depth_efolding","float",(/allpfts!0/))
  filevardef(fout,"rf_l3s2_bgc","float",(/allpfts!0/))
  filevardef(fout,"rf_s2s1_bgc","float",(/allpfts!0/))
  filevardef(fout,"rf_s3s1_bgc","float",(/allpfts!0/))
  filevardef(fout,"k_decay_pom","float",(/allpfts!0/))
  filevardef(fout,"k_decay_humus","float",(/allpfts!0/))
  filevardef(fout,"k_decay_cwd","float",(/allpfts!0/))
  filevardef(fout,"k_decay_lwd","float",(/allpfts!0/))
  filevardef(fout,"k_decay_fwd","float",(/allpfts!0/))
  filevardef(fout,"k_m_o2_bgc","float",(/allpfts!0/))
  filevardef(fout,"organic_max","float",(/allpfts!0/))
  filevardef(fout,"nitrif_n2o_loss_frac","float",(/allpfts!0/))
  filevardef(fout,"rf_l1s1_bgc","float",(/allpfts!0,surfsoil!0/))
  filevardef(fout,"rf_l2s1_bgc","float",(/allpfts!0,surfsoil!0/))
  filevardef(fout,"rf_s1s2a_bgc","float",(/allpfts!0,surfsoil!0/))
  filevardef(fout,"rf_s1s2b_bgc","float",(/allpfts!0,surfsoil!0/))
  filevardef(fout,"k_decay_lmet","float",(/allpfts!0,surfsoil!0/))
  filevardef(fout,"k_decay_lcel","float",(/allpfts!0,surfsoil!0/))
  filevardef(fout,"k_decay_llig","float",(/allpfts!0,surfsoil!0/))
  filevardef(fout,"k_decay_mic","float",(/allpfts!0,surfsoil!0/))
  filevardef(fout,"vmax_minp_soluble_to_secondary","float",(/pfts!0/))
  filevardef(fout,"frac_p_sec_to_sol","float",(/pfts!0/))
  filevardef(fout,"minp_secondary_decay","float",(/pfts!0/))
  filevardef(fout,"k_decay_dom","float",(/allpfts!0/))
  filevardef(fout,"km_mic_som","float",(/allpfts!0/))
  filevardef(fout,"km_decomp_nh4","float",(/allpfts!0/))
  filevardef(fout,"km_decomp_no3","float",(/allpfts!0/))
  filevardef(fout,"km_decomp_p","float",(/allpfts!0/))
  filevardef(fout,"km_den","float",(/allpfts!0/))
  filevardef(fout,"km_nit","float",(/allpfts!0/))
  filevardef(fout,"vmax_decomp_n","float",(/allpfts!0/))
  filevardef(fout,"vmax_decomp_p","float",(/allpfts!0/))
  filevardef(fout,"vmax_nit","float",(/allpfts!0/))
  filevardef(fout,"vmax_den","float",(/allpfts!0/))
;  filevardef(fout,"rock_name","float",(/lithos!0/))
  filevardef(fout,"b_weath","float",(/lithos!0/))
  filevardef(fout,"f_shield","float",(/lithos!0/))
  filevardef(fout,"P_weip","float",(/lithos!0/))
  filevardef(fout,"E_weath","float",(/lithos!0/))
  filevardef(fout,"T_ref_weath","float",(/alllithos!0/))
;---define variable attributes

   fatt=True
   fatt@long_name="fraction of cellulose carbon in coarse woody roots"
   fatt@units = "none"
   filevarattdef(fout,"cwd_fcel", fatt)

   fatt@long_name="fraction of lignin carbon in coarse woody roots"
   fatt@units = "none"
   filevarattdef(fout,"cwd_flig", fatt)

   fatt@long_name="fraction of cellulose carbon in fine woody branches"
   fatt@units = "none"
   filevarattdef(fout,"fwd_fcel", fatt)

   fatt@long_name="fraction of lignin carbon in fine woody branches"
   fatt@units = "none"
   filevarattdef(fout,"fwd_flig", fatt)

   fatt@long_name="fraction of cellulose carbon in large woody branches"
   fatt@units = "none"
   filevarattdef(fout,"lwd_fcel", fatt)

   fatt@long_name="fraction of lignin carbon in large woody branches"
   fatt@units = "none"
   filevarattdef(fout,"lwd_flig", fatt)

   fatt@long_name="surface tension of water"
   fatt@units = "N m-1"
   filevarattdef(fout,"surface_tension_water", fatt)

   fatt@long_name="parameter a for arah-vinten anaerobic fraction"
   fatt@units = "none"
   filevarattdef(fout,"rij_kro_a", fatt)

   fatt@long_name="parameter alpha for arah-vinten anaerobic fraction"
   fatt@units = "none"
   filevarattdef(fout,"rij_kro_alpha", fatt)

   fatt@long_name="parameter beta for arah-vinten anaerobic fraction"
   fatt@units = "none"
   filevarattdef(fout,"rij_kro_beta", fatt)

   fatt@long_name="parameter gamma for arah-vinten anaerobic fraction"
   fatt@units = "none"
   filevarattdef(fout,"rij_kro_gamma", fatt)

   fatt@long_name="parameter delta for arah-vinten anaerobic fraction"
   fatt@units = "none"
   filevarattdef(fout,"rij_kro_delta", fatt)

   fatt@long_name="minimum suction pressure with active aerobic respiration"
   fatt@units = "M Pa"
   filevarattdef(fout,"minpsi_bgc", fatt)

   fatt@long_name="maximum nitrification rate"
   fatt@units = "1/s"
   filevarattdef(fout,"k_nitr_max", fatt)

   fatt@long_name="Increase in reaction rate per 10 K increase in temperature for unfrozen condition"
   fatt@units = "none"
   filevarattdef(fout,"Q10", fatt)

   fatt@long_name="Increase in reaction rate per 10 K increase in temperature for frozen condition"
   fatt@units = "none"
   filevarattdef(fout,"froz_q10", fatt)

   fatt@long_name="The relaxing depth to impose the decrease of respiration along depth"
   fatt@units = "1/m"
   filevarattdef(fout,"decomp_depth_efolding", fatt)

   fatt@long_name="respiraiton efficiency during l3->s2"
   fatt@units = "none"
   filevarattdef(fout,"rf_l3s2_bgc", fatt)

   fatt@long_name="respiration efficiency during s2->s1"
   fatt@units = "none"
   filevarattdef(fout,"rf_s2s1_bgc", fatt)

   fatt@long_name="respiration efficiency during s3->s1"
   fatt@units = "none"
   filevarattdef(fout,"rf_s3s1_bgc", fatt)

   fatt@long_name="intrinsic decay rate of pom"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_pom", fatt)

   fatt@long_name="intrinsic decay rate of humus"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_humus", fatt)

   fatt@long_name="intrinsic decay rate of coarse woody root"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_cwd", fatt)

   fatt@long_name="intrinsic decay rate of large woody debries"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_lwd", fatt)

   fatt@long_name="intrinsic decay rate of fine woody root"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_fwd", fatt)

   fatt@long_name="half saturation constant for MM-based O2 binding"
   fatt@units = "mol/m3"
   filevarattdef(fout,"k_m_o2_bgc", fatt)

   fatt@long_name="maximum density of organic matter"
   fatt@units = "kg m-3"
   filevarattdef(fout,"organic_max", fatt)

   fatt@long_name="respiration efficiency during l1->s1"
   fatt@units = "none"
   filevarattdef(fout,"rf_l1s1_bgc", fatt)

   fatt@long_name="respiration efficiency during l2->s1"
   fatt@units = "none"
   filevarattdef(fout,"rf_l2s1_bgc", fatt)

   fatt@long_name="regression parameter a for computing respiration efficiency during s1->s2"
   fatt@units = "none"
   filevarattdef(fout,"rf_s1s2a_bgc", fatt)

   fatt@long_name="half saturation constatnt for microbial biomass for modeling linear decomposition"
   fatt@units = "mol C m-3"
   filevarattdef(fout,"km_mic_som", fatt)


   fatt@long_name="regression parameter b for computing respiration efficiency during s1->s2"
   fatt@units = "none"
   filevarattdef(fout,"rf_s1s2b_bgc", fatt)

   fatt@long_name="intrinsic decay rate of pool dom"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_dom", fatt)

   fatt@long_name="intrinsic decay rate of pool lit1"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_lmet", fatt)

   fatt@long_name="intrinsic decay rate of pool lit2"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_lcel", fatt)

   fatt@long_name="intrinsic decay rate of pool lit3"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_llig", fatt)

   fatt@long_name="intrinsic decay rate of pool _mic"
   fatt@units = "1/year"
   filevarattdef(fout,"k_decay_mic", fatt)

   fatt@long_name="maximum conversion rate of soprobed inorganic P into secondary P"
   fatt@units = "1/s"
   filevarattdef(fout,"vmax_minp_soluble_to_secondary", fatt)

   fatt@long_name="fraction of secondary inorganic P into soluble P during secondary inorganic P mobilization"
   fatt@units = "none"
   filevarattdef(fout,"frac_p_sec_to_sol", fatt)

   fatt@long_name="secondary inorganic P mobilization rate"
   fatt@units = "1/s"
   filevarattdef(fout,"minp_secondary_decay", fatt)

   fatt@long_name="fraction of nitrification product as n2o"
   fatt@units = "none"
   filevarattdef(fout,"nitrif_n2o_loss_frac", fatt)

   fatt@long_name="NH4 affinity for nitrification"
   fatt@units = "g N/m3"
   filevarattdef(fout,"km_nit", fatt)

   fatt@long_name="NO3 affinity for denitrification"
   fatt@units = "g N/m3"
   filevarattdef(fout,"km_den", fatt)

   fatt@long_name="NH4 affinity for decomposition"
   fatt@units = "g N/m3"
   filevarattdef(fout,"km_decomp_nh4", fatt)

   fatt@long_name="NO3 affinity for decomposition"
   fatt@units = "g N/m3"
   filevarattdef(fout,"km_decomp_no3", fatt)

   fatt@long_name="P affinity for decomposition"
   fatt@units = "g P/m3"
   filevarattdef(fout,"km_decomp_p", fatt)

   fatt@long_name="max N uptake rate for decomposers"
   fatt@units = "1/hour"
   filevarattdef(fout,"vmax_decomp_n", fatt)

   fatt@long_name="max P uptake rate for decomposers"
   fatt@units = "1/hour"
   filevarattdef(fout,"vmax_decomp_p", fatt)

   fatt@long_name="specific max nitrification rate"
   fatt@units = "1/day"
   filevarattdef(fout,"vmax_nit", fatt)

   fatt@long_name="specific max denitrification rate"
   fatt@units = "1/day"
   filevarattdef(fout,"vmax_den", fatt)

;   fatt@long_name="lithological name of the rocks"
;   fatt@units = "none"
;   filevarattdef(fout,"rock_name", fatt)

   fatt@long_name="weathering b parameter"
   fatt@units = "ton km-2 mm-1"
   filevarattdef(fout,"b_weath", fatt)

   fatt@long_name="weathering shileding factor"
   fatt@units = "none"
   filevarattdef(fout,"f_shield", fatt)

   fatt@long_name="P-content in % of the wathering flux"
   fatt@units = "none"
   filevarattdef(fout,"P_weip", fatt)

   fatt@long_name="activation energy of weathering"
   fatt@units = "kJ/mol"
   filevarattdef(fout,"E_weath", fatt)

   fatt@long_name="Reference temperature for rock weathering"
   fatt@units = "K"
   filevarattdef(fout,"T_ref_weath", fatt)

;   if (False) then
   fout->cwd_fcel=(/0.76/)
   fout->cwd_flig=(/0.26/)

   fout->fwd_fcel=(/0.78/)
   fout->fwd_flig=(/0.22/)

   fout->lwd_fcel=(/0.76/)
   fout->lwd_flig=(/0.26/)

   fout->surface_tension_water=(/0.073/)
   fout->Q10=(/2.0/)
   fout->froz_q10=(/10.0/)
   fout->decomp_depth_efolding=(/0.7/)
   fout->rf_l1s1_bgc=(/(/0.55/),(/0.55/)/)
   fout->rf_l2s1_bgc=(/(/0.45/),(/0.55/)/)
   fout->rf_l3s2_bgc=(/0.3/)
   fout->rf_s2s1_bgc=(/0.55/)
   fout->rf_s3s1_bgc=(/0.55/)
   fout->rf_s1s2a_bgc=(/(/0.60/),(/0.17/)/)
   fout->rf_s1s2b_bgc=(/(/0.00/),(/0.68/)/)
   fout->k_decay_lmet=(/(/14.8/),(/18.5/)/)
   fout->k_decay_lcel=(/(/3.9/),(/4.9/)/)
   fout->k_decay_llig=(/(/3.9/),(/4.9/)/)
   fout->k_decay_mic=(/(/6.7/),(/7.3/)/)
   fout->k_decay_pom=(/0.2/)
   fout->k_decay_humus=(/0.0045/)
   fout->k_decay_cwd=(/0.6/)
   fout->k_decay_fwd=(/1.5/)
   fout->k_decay_lwd=(/0.5/)
   fout->k_decay_dom=(/8.0/)
   fout->km_mic_som=(/3.0/)
   fout->k_m_o2_bgc = (/0.22/)
   fout->k_nitr_max = (/1.1574074e-06/)
   fout->minpsi_bgc =(/-10.0/)
   aa=1.5*10^(-7.48)
   fout->nitrif_n2o_loss_frac=(/1.e-4/)
   fout->organic_max = (/160.0/)
   ;the Arah and Vinten value is 1.8e-4
   fout->rij_kro_a  = (/aa/)
   ;the Arah and Vinten value is 0.77
   fout->rij_kro_alpha=(/1.26/)
   ;the Arah and Vinten value is 0.37
   fout->rij_kro_beta =(/0.6/)
   ;the Arah and Vinten value is 0.32
   fout->rij_kro_gamma = (/0.6/)
   ;the Arah and Vinten value is 0.45
   fout->rij_kro_delta = (/0.85/)

   fout->vmax_minp_soluble_to_secondary=(/1.e-7,1.e-7,1.e-7,1.e-7,1.e-7,1.e-7,\
                                          1.e-7,1.e-7,1.e-7,1.e-7,1.e-7,1.e-7,\
                                          1.e-7,1.e-7,1.e-7,1.e-7,1.e-7,1.e-7,\
                                          1.e-7,1.e-7,1.e-7,1.e-7,1.e-7,1.e-7,\
                                          1.e-7/)
   fout->frac_p_sec_to_sol = (/0.95,0.95,0.95,0.95,0.95,0.95,\
                               0.95,0.95,0.95,0.95,0.95,0.95,\
                               0.95,0.95,0.95,0.95,0.95,0.95,\
                               0.95,0.95,0.95,0.95,0.95,0.95,\
                               0.95/)

   fout->minp_secondary_decay = (/2.e-8,2.e-8,2.e-8,2.e-8,2.e-8,2.e-8,\
                                  2.e-8,2.e-8,2.e-8,2.e-8,2.e-8,2.e-8,\
                                  2.e-8,2.e-8,2.e-8,2.e-8,2.e-8,2.e-8,\
                                  2.e-8,2.e-8,2.e-8,2.e-8,2.e-8,2.e-8,\
                                  2.e-8/)
   fout->km_nit = (/0.76/)
   fout->km_den = (/0.11/)
   fout->km_decomp_nh4 = (/0.18/)
   fout->km_decomp_no3 = (/0.41/)
   fout->km_decomp_p   = (/0.2/)
   fout->vmax_decomp_n = (/5.e-3/)
   fout->vmax_decomp_p = (/1.e-3/)
   fout->vmax_nit = (/0.67/)
   fout->vmax_den = (/1.8/)
   fout->b_weath =(/0.021333,0.04054,0.027799,0.04054,0.057721,0.151243, \
                    0.020762,0.032036,0.02691,0.04054,0.      ,0.076876, \
                    0.02691,0.151243,0.,0./)
   fout->f_shield=(/0.5,0.7,0.4,0.6,0.5,0.8,0.5,0.5,0.4,0.6,0.,0.8,0.8,1.0,0.,0./)
   fout->E_weath =(/60.0,50.0,60.0,50.,60,60,60,60,60,50,0,46,60,60,0,0/)
   fout->P_weip =(/0.052,0.135,0.052,0.227,0.052,0.048,0.031,0.079,0.07,0.135,0,0.227,0.227,0.026,0,0/)
   fout->T_ref_weath=(/284.15/)
;   end if
   delete(fout)
   end_write_time = get_cpu_time()

end
