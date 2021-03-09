#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
 #include <math.h>
 #include <admodel.h>
  #include <time.h>
 ofstream CheckFile;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
 ofstream R_out;
 dvariable gammln(dvariable& xx)
 {
   RETURN_ARRAYS_INCREMENT();	
        dvariable x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
               0.1208650973866179e-2,-0.5395239384953e-5}; 
        int j;
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) 
          {
            ser += cof[j]/(y+1.);
          }
         dvariable value_=-tmp+log(2.5066282746310005*ser/x);
    RETURN_ARRAYS_DECREMENT();         
         return(value_); 
 }
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <2016sc.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_post = new ofstream("eval.csv");
 call_no = 0;
 spmo=0.;																							// spmo=deviation in fraction of year from time of fishery to mating 
  styr.allocate("styr");
  endyr.allocate("endyr");
  nirec.allocate("nirec");
  nlenm.allocate("nlenm");
  nobs_fish.allocate("nobs_fish");
  yrs_fish.allocate(1,nobs_fish,"yrs_fish");
  nsamples_fish.allocate(1,2,1,nobs_fish,"nsamples_fish");
  nobs_fish_discf.allocate("nobs_fish_discf");
  yrs_fish_discf.allocate(1,nobs_fish_discf,"yrs_fish_discf");
  nsamples_fish_discf.allocate(1,nobs_fish_discf,"nsamples_fish_discf");
  nobs_fish_discm.allocate("nobs_fish_discm");
  yrs_fish_discm.allocate(1,nobs_fish_discm,"yrs_fish_discm");
  nsamples_fish_discm.allocate(1,2,1,nobs_fish_discm,"nsamples_fish_discm");
  nobs_trawl.allocate("nobs_trawl");
  yrs_trawl.allocate(1,nobs_trawl,"yrs_trawl");
  nsamples_trawl.allocate(1,2,1,nobs_trawl,"nsamples_trawl");
  nobs_srv1.allocate("nobs_srv1");
  yrs_srv1.allocate(1,nobs_srv1,"yrs_srv1");
  nobs_srv1_length.allocate("nobs_srv1_length");
  yrs_srv1_length.allocate(1,nobs_srv1_length,"yrs_srv1_length");
  nsamples_srv1_length.allocate(1,2,1,2,1,2,1,nobs_srv1_length,"nsamples_srv1_length");
  yrs_srv2.allocate("yrs_srv2");
  nsamples_srv2_length.allocate(1,2,1,2,1,2,"nsamples_srv2_length");
  obs_p_srv2_lend.allocate(1,2,1,2,1,2,1,nlenm,"obs_p_srv2_lend");
  obs_srv2.allocate(1,2,"obs_srv2");
  obs_srv2_cv.allocate(1,2,1,2,"obs_srv2_cv");
  yrs_srv10.allocate("yrs_srv10");
  nsamples_srv10_length.allocate(1,2,1,2,1,2,"nsamples_srv10_length");
  obs_p_srv10_lend.allocate(1,2,1,2,1,2,1,nlenm,"obs_p_srv10_lend");
  obs_srv10.allocate(1,2,"obs_srv10");
  obs_srv10_cv.allocate(1,2,1,2,"obs_srv10_cv");
  obs_p_srv1_lend.allocate(1,2,1,2,1,2,1,nobs_srv1_length,1,nlenm,"obs_p_srv1_lend");
  obs_p_fish_retd.allocate(1,2,1,nobs_fish,1,nlenm,"obs_p_fish_retd");
  obs_p_fish_discfd.allocate(1,nobs_fish_discf,1,nlenm,"obs_p_fish_discfd");
  obs_p_fish_discmd.allocate(1,2,1,nobs_fish_discm,1,nlenm,"obs_p_fish_discmd");
  obs_p_trawld.allocate(1,2,1,nobs_trawl,1,nlenm,"obs_p_trawld");
  catch_numbers.allocate(styr,endyr,"catch_numbers");
 catch_numbers /= 1000.;                      												 //retained catch number of crab 
  catch_ret.allocate(styr,endyr,"catch_ret");
 catch_ret /= 1000.;
  catch_odisc.allocate(1,2,styr,endyr,"catch_odisc");
 catch_odisc /= 1000.;
  catch_trawl.allocate(styr,endyr,"catch_trawl");
 catch_trawl /= 1000.;
  obs_srv1.allocate(1,nobs_srv1,"obs_srv1");
 obs_srv1 /= 1000.;
  cv_srv1o.allocate(1,2,1,nobs_srv1,"cv_srv1o");
  maturity_logistic.allocate(1,nlenm,"maturity_logistic");
  maturity_average.allocate(1,2,1,nlenm,"maturity_average");
  maturity_old_average.allocate(1,2,1,nlenm,"maturity_old_average");
  maturity.allocate(1,2,styr,endyr,1,nlenm,"maturity");
  maturity_old.allocate(1,2,styr,endyr,1,nlenm,"maturity_old");
  cv_mean_length_obs.allocate(1,2,1,2,"cv_mean_length_obs");
  length_bins.allocate(1,nlenm,"length_bins");
  catch_midptIn.allocate(styr,endyr,"catch_midptIn");
  meanlength.allocate(1,13,"meanlength");
  catch_fracmature.allocate(1,22,"catch_fracmature");
  cpue.allocate(styr,endyr,"cpue");
  catch_ghl.allocate(styr,endyr,"catch_ghl");
  nobs_growf.allocate("nobs_growf");
  femalegrowdatx.allocate(1,nobs_growf,"femalegrowdatx");
  femalegrowdaty.allocate(1,nobs_growf,"femalegrowdaty");
  nobs_growm.allocate("nobs_growm");
  malegrowdatx.allocate(1,nobs_growm,"malegrowdatx");
  malegrowdaty.allocate(1,nobs_growm,"malegrowdaty");
cout<<" "<<endl;  	 
cout<<obs_srv2_cv<<endl;    
cout<<" "<<endl;  	 
cout<<catch_ret<<endl;  	
cout<<" "<<endl;  	 
cout<<catch_odisc<<endl;  	  
cout<<catch_trawl<<endl;  	 
cout<<" "<<endl;  	 
cout<<catch_fracmature<<endl;  	 
cout<<" "<<endl;  	 
cout<<malegrowdaty<<endl;  
cout<<" "<<endl;  	 
cout<<nobs_growm<<endl;  		 
 ad_comm::change_datafile_name("2016sc.ctl");
  styr_fut.allocate("styr_fut");
  endyr_fut.allocate("endyr_fut");
  nsellen.allocate("nsellen");
  nsellen_srv1.allocate("nsellen_srv1");
  p_const.allocate("p_const");
  q1.allocate("q1");
  M_in.allocate(1,2,"M_in");
  M_matn_in.allocate(1,2,"M_matn_in");
  M_mato_in.allocate(1,2,"M_mato_in");
  m_disc.allocate("m_disc");
  m_trawl.allocate("m_trawl");
  median_rec.allocate(1,2,"median_rec");
  median_rec_yrs.allocate("median_rec_yrs");
  var_rec_obs.allocate("var_rec_obs");
  sd_var_rec.allocate("sd_var_rec");
  sep_rec_dev.allocate("sep_rec_dev");
  sel_som.allocate(1,5,"sel_som");
  sel_avg_Nyrs.allocate("sel_avg_Nyrs");
  alpha_wt_imm_f.allocate("alpha_wt_imm_f");
  beta_wt_imm_f.allocate("beta_wt_imm_f");
  alpha_wt_mat_f.allocate("alpha_wt_mat_f");
  beta_wt_mat_f.allocate("beta_wt_mat_f");
  alpha_wt_m.allocate("alpha_wt_m");
  beta_wt_m.allocate("beta_wt_m");
  linff_obs.allocate("linff_obs");
  sd_linff.allocate("sd_linff");
  linfm_obs.allocate("linfm_obs");
  sd_linfm.allocate("sd_linfm");
  growthkf_obs.allocate("growthkf_obs");
  sd_growthkf.allocate("sd_growthkf");
  growthkm_obs.allocate("growthkm_obs");
  sd_growthkm.allocate("sd_growthkm");
  af_obs.allocate("af_obs");
  sd_af.allocate("sd_af");
  am_obs.allocate("am_obs");
  sd_am.allocate("sd_am");
  bf_obs.allocate("bf_obs");
  sd_bf.allocate("sd_bf");
  bm_obs.allocate("bm_obs");
  sd_bm.allocate("sd_bm");
  a1_obs.allocate("a1_obs");
  sd_a1.allocate("sd_a1");
  b1_obs.allocate("b1_obs");
  sd_b1.allocate("sd_b1");
  sd_meetpt.allocate("sd_meetpt");
  var_last_obs.allocate("var_last_obs");
  sd_var_last.allocate("sd_var_last");
  mate_ratio.allocate("mate_ratio");
  fraction_new_error.allocate("fraction_new_error");
  nages.allocate("nages");
  matest_n.allocate("matest_n");
  matestm_n.allocate("matestm_n");
  wt_like.allocate(1,8,"wt_like");
  like_wght.allocate(1,7,"like_wght");
  like_wght_mbio.allocate("like_wght_mbio");
  like_wght_rec.allocate("like_wght_rec");
  like_wght_recf.allocate("like_wght_recf");
  like_wght_sexr.allocate("like_wght_sexr");
  like_wght_sel50.allocate("like_wght_sel50");
  like_wght_fph1.allocate("like_wght_fph1");
  like_wght_fph2.allocate("like_wght_fph2");
  like_wght_fdev.allocate("like_wght_fdev");
  wght_total_catch.allocate("wght_total_catch");
  wght_female_potcatch.allocate("wght_female_potcatch");
  cpue_cv.allocate("cpue_cv");
  wt_lmlike.allocate("wt_lmlike");
  old_shell_constraint.allocate("old_shell_constraint");
  disclen_mult.allocate(1,2,"disclen_mult");
  selsmo_wght.allocate("selsmo_wght");
  femSel_wght.allocate("femSel_wght");
  init_yr_len_smooth_f.allocate("init_yr_len_smooth_f");
  init_yr_len_smooth_m.allocate("init_yr_len_smooth_m");
  natm_mult_wght.allocate("natm_mult_wght");
  natm_mult_var.allocate("natm_mult_var");
  natm_Immult_wght.allocate("natm_Immult_wght");
  natm_Immult_var.allocate("natm_Immult_var");
  smooth_mat_wght.allocate("smooth_mat_wght");
  mat_est_wght.allocate("mat_est_wght");
  mat_est_vs_obs_sd.allocate("mat_est_vs_obs_sd");
  smooth_mat_wght_f.allocate("smooth_mat_wght_f");
  mat_est_vs_obs_sd_f.allocate("mat_est_vs_obs_sd_f");
  growth_data_wght_m.allocate("growth_data_wght_m");
  growth_data_wght_f.allocate("growth_data_wght_f");
  extra_wght_ind_m.allocate("extra_wght_ind_m");
  smooth_disc_catch.allocate("smooth_disc_catch");
  disc_catch_wght_fem.allocate("disc_catch_wght_fem");
  fmort_phase.allocate("fmort_phase");
  rec_phase.allocate("rec_phase");
  growth_phase.allocate("growth_phase");
  growth_phase2.allocate("growth_phase2");
  maturity_phase.allocate("maturity_phase");
  natM_phase.allocate("natM_phase");
  phase_moltingp.allocate("phase_moltingp");
  phase_fishsel.allocate("phase_fishsel");
  survsel_phase.allocate("survsel_phase");
  survsel1_phase.allocate("survsel1_phase");
  phase_fut.allocate("phase_fut");
  phase_logistic_sel.allocate("phase_logistic_sel");
  phase_selcoffs.allocate("phase_selcoffs");
  growth_switch.allocate("growth_switch");
  somertonsel.allocate("somertonsel");
  monot_sel.allocate("monot_sel");
  monot_sel_srv1.allocate("monot_sel_srv1");
  maturity_switch.allocate("maturity_switch");
  f_penalties.allocate("f_penalties");
  retro_years.allocate("retro_years");
 Nproj = 101;
   cout<<"to local calcs"<<endl;
   styr_rec=styr-nirec; 																									  //year to start estimating recruits to get initial age comp
   if(nsellen>nlenm) nsellen=nlenm; 																				 //make sure nselages not greater than nages
   if(nsellen_srv1>nlenm) nsellen_srv1=nlenm; 															 //same as above for survey
   obs_srv1=obs_srv1*1000000;              																		 //survey numbers read in are millions of crab
   obs_srv2=obs_srv2*1000000;              																		  //survey numbers read in are millions of crab
   catch_ret=catch_ret/2204.6; 
   cout<<"end of local calcs"<<endl;
}

void model_parameters::initializationfunction(void)
{
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  af.allocate(-100,10.0,growth_phase,"af");
  am.allocate(-50.0,10.0,growth_phase,"am");
  bf.allocate(1.0,10.0,growth_phase,"bf");
  bm.allocate(1.0,5.0,growth_phase,"bm");
  bf1.allocate(1.0,2.0,growth_phase2,"bf1");
  deltaf.allocate(5.0,50.0,growth_phase2,"deltaf");
  st_gr.allocate(0.5,0.5,-2,"st_gr");
  growth_beta.allocate(1,2,0.749,0.751,-2,"growth_beta");
  matest50f.allocate(30,60,-4,"matest50f");
  matestslpf.allocate(0.05,0.5,-4,"matestslpf");
  matest50m.allocate(60,110,-4,"matest50m");
  matestslpm.allocate(0.02,0.5,-4,"matestslpm");
  mateste.allocate(1,matestm_n,-6.0,-0.000000001,maturity_phase,"mateste");
  matestfe.allocate(1,matest_n,-6.0,-0.000000001,maturity_phase,"matestfe");
  moltp_af.allocate(0.01,3.0,-6,"moltp_af");
  moltp_bf.allocate(20,200,-6,"moltp_bf");
  moltp_am.allocate(0.04,3.0,-5,"moltp_am");
  moltp_bm.allocate(130.0,300.0,-5,"moltp_bm");
  moltp_ammat.allocate(.0025,3.0,phase_moltingp,"moltp_ammat");
  moltp_bmmat.allocate(1,120,phase_moltingp,"moltp_bmmat");
  mean_log_rec_f.allocate("mean_log_rec_f");
  mean_log_rec_m.allocate(sep_rec_dev,"mean_log_rec_m");
  rec_devf.allocate(styr,endyr,-15,15,rec_phase,"rec_devf");
  rec_devm.allocate(styr,endyr,-15,15,sep_rec_dev*rec_phase,"rec_devm");
  alpha1_rec.allocate(11.49,11.51,-8,"alpha1_rec");
  beta_rec.allocate(3.99,4.01,-8,"beta_rec");
  mnatlen_styr.allocate(1,2,1,nlenm,-3,15,1,"mnatlen_styr");
  fnatlen_styr.allocate(1,2,1,12,-10,15,1,"fnatlen_styr");
  log_avg_fmort.allocate(1,"log_avg_fmort");
  fmort_dev.allocate(styr,endyr-1,-5,5,fmort_phase,"fmort_dev");
  log_avg_fmortdf.allocate(-8.0,-0.0001,fmort_phase,"log_avg_fmortdf");
  fmortdf_dev.allocate(styr,endyr-1,-15,15,fmort_phase,"fmortdf_dev");
  log_avg_fmortt.allocate(-8.0,-0.0001,fmort_phase,"log_avg_fmortt");
  fmortt_dev_era1.allocate(styr,1991,-15,15,fmort_phase,"fmortt_dev_era1");
  fmortt_dev_era2.allocate(1992,endyr-1,-15,15,fmort_phase,"fmortt_dev_era2");
  discard_mult.allocate(0.999,1.01,-fmort_phase,"discard_mult");
  Fem_F_prop_constant.allocate(0.0000001,1,fmort_phase,"Fem_F_prop_constant");
  log_avg_sel50_mn.allocate(4,5.0,phase_logistic_sel,"log_avg_sel50_mn");
  log_sel50_dev_mn.allocate(styr,endyr,-5,5,-phase_logistic_sel+2,"log_sel50_dev_mn");
  log_avg_sel50_mo.allocate(4,5.0,-phase_logistic_sel,"log_avg_sel50_mo");
  log_sel50_dev_mo.allocate(styr,endyr,-5,5,-phase_logistic_sel+2,"log_sel50_dev_mo");
  fish_slope_mn.allocate(0.1,0.5,phase_logistic_sel,"fish_slope_mn");
  fish_slope_mo.allocate(0.1,0.8,-phase_logistic_sel,"fish_slope_mo");
  fish_fit_slope_mn.allocate(.05,.5,phase_logistic_sel,"fish_fit_slope_mn");
  fish_fit_sel50_mn.allocate(85.0,120.0,phase_logistic_sel,"fish_fit_sel50_mn");
  fish_fit_slope_mo.allocate(.05,.5,-phase_logistic_sel,"fish_fit_slope_mo");
  fish_fit_sel50_mo.allocate(90.0,300.0,-phase_logistic_sel,"fish_fit_sel50_mo");
  fish_slope_mo2.allocate(1.9,2.0,phase_fishsel,"fish_slope_mo2");
  fish_sel50_mo2.allocate(159.0,160.0,phase_fishsel,"fish_sel50_mo2");
  fish_slope_mn2.allocate(.01,2.0,phase_fishsel,"fish_slope_mn2");
  fish_sel50_mn2.allocate(100.0,160.0,phase_fishsel,"fish_sel50_mn2");
  fish_disc_slope_f.allocate(0.10,0.70,phase_logistic_sel,"fish_disc_slope_f");
  fish_disc_sel50_f.allocate(1,5,phase_logistic_sel,"fish_disc_sel50_f");
  log_dev_50f.allocate(1995,endyr-1,-5,5,-phase_logistic_sel,"log_dev_50f");
  fish_disc_slope_tf.allocate(.01,.3,phase_logistic_sel,"fish_disc_slope_tf");
  fish_disc_sel50_tf.allocate(30.,120.0,phase_logistic_sel,"fish_disc_sel50_tf");
  srv1_q.allocate(0.2,1.000,-survsel1_phase+1,"srv1_q");
  srv1_q_f.allocate(0.2,1.000,-survsel1_phase+1,"srv1_q_f");
  srv1_sel95.allocate(30.0,150.0,-survsel1_phase,"srv1_sel95");
  srv1_sel50.allocate(0.0,150.0,-survsel1_phase,"srv1_sel50");
  srv2_q.allocate(0.2,1.000,survsel1_phase+1,"srv2_q");
  srv2_q_f.allocate(0.2,1.000,survsel1_phase+1,"srv2_q_f");
  srv2_sel95.allocate(50.0,160.0,survsel1_phase,"srv2_sel95");
  srv2_sel50.allocate(0.0,80.0,survsel1_phase,"srv2_sel50");
  srv3_q.allocate(0.20,1.0000,survsel1_phase,"srv3_q");
  srv3_sel95.allocate(40.0,200.0,survsel1_phase,"srv3_sel95");
  srv3_sel50.allocate(25.0,90.0,survsel1_phase,"srv3_sel50");
  srv3_q_f.allocate(0.20,1.0000,survsel1_phase,"srv3_q_f");
  srv3_sel95_f.allocate(40.0,150.0,survsel1_phase+1,"srv3_sel95_f");
  srv3_sel50_f.allocate(0.0,90.0,survsel1_phase+1,"srv3_sel50_f");
  srvind_q.allocate(0.10,1.000,survsel1_phase+2,"srvind_q");
  srvind_q_f.allocate(0.01,1.000,survsel1_phase+2,"srvind_q_f");
  srvind_sel95_f.allocate(5.0,120.0,survsel1_phase,"srvind_sel95_f");
  srvind_sel50_f.allocate(1.0,110.0,survsel1_phase,"srvind_sel50_f");
  srv10ind_q.allocate(0.1,1.000,survsel1_phase+2,"srv10ind_q");
  srv10ind_q_f.allocate(0.01,1.00,survsel1_phase+2,"srv10ind_q_f");
  srv10ind_sel95_f.allocate(0.05,50.0,-survsel1_phase,"srv10ind_sel95_f");
  srv10ind_sel50_f.allocate(0.0,50.0,-survsel1_phase,"srv10ind_sel50_f");
  selsmo10ind.allocate(1,nlenm,-20,20,survsel1_phase,"selsmo10ind");
  selsmo09ind.allocate(1,nlenm,-20,20,survsel1_phase,"selsmo09ind");
  Mmult_imat.allocate(0.2000,2.000001,natM_phase,"Mmult_imat");
  Mmult.allocate(0.20000,2.0000001,natM_phase,"Mmult");
  Mmultf.allocate(0.20000,2.0000001,natM_phase,"Mmultf");
  cpueq.allocate(0.0000877,0.00877,5,"cpueq");
  proprecn.allocate(1.0,1.0,-2,"proprecn");
  wtf.allocate(1,2,1,nlenm,"wtf");
  #ifndef NO_AD_INITIALIZE
    wtf.initialize();
  #endif
  wtm.allocate(1,nlenm,"wtm");
  #ifndef NO_AD_INITIALIZE
    wtm.initialize();
  #endif
  fish_sel50_mn.allocate(styr,endyr+Nproj,"fish_sel50_mn");
  #ifndef NO_AD_INITIALIZE
    fish_sel50_mn.initialize();
  #endif
  fish_sel50_mo.allocate(styr,endyr+Nproj,"fish_sel50_mo");
  #ifndef NO_AD_INITIALIZE
    fish_sel50_mo.initialize();
  #endif
  log_sel_fish.allocate(1,nlenm,"log_sel_fish");
  #ifndef NO_AD_INITIALIZE
    log_sel_fish.initialize();
  #endif
  log_sel_fish_disc.allocate(1,2,1,nlenm,"log_sel_fish_disc");
  #ifndef NO_AD_INITIALIZE
    log_sel_fish_disc.initialize();
  #endif
  log_sel_srv1.allocate(1,2,1,nlenm,"log_sel_srv1");
  #ifndef NO_AD_INITIALIZE
    log_sel_srv1.initialize();
  #endif
  sel.allocate(1,2,styr,endyr+Nproj,1,nlenm,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  sel_discf.allocate(1995,endyr+Nproj,1,nlenm,"sel_discf");
  #ifndef NO_AD_INITIALIZE
    sel_discf.initialize();
  #endif
  sel_discf_e.allocate(1,nlenm,"sel_discf_e");
  #ifndef NO_AD_INITIALIZE
    sel_discf_e.initialize();
  #endif
  sel_fit.allocate(1,2,styr,endyr+Nproj,1,nlenm,"sel_fit");
  #ifndef NO_AD_INITIALIZE
    sel_fit.initialize();
  #endif
  sel_ret.allocate(1,2,styr,endyr+Nproj,1,nlenm,"sel_ret");
  #ifndef NO_AD_INITIALIZE
    sel_ret.initialize();
  #endif
  sel_trawl.allocate(1,2,1,nlenm,"sel_trawl");
  #ifndef NO_AD_INITIALIZE
    sel_trawl.initialize();
  #endif
  sel_srv1.allocate(1,2,1,nlenm,"sel_srv1");
  #ifndef NO_AD_INITIALIZE
    sel_srv1.initialize();
  #endif
  sel_srv2.allocate(1,2,1,nlenm,"sel_srv2");
  #ifndef NO_AD_INITIALIZE
    sel_srv2.initialize();
  #endif
  sel_srv3.allocate(1,2,1,nlenm,"sel_srv3");
  #ifndef NO_AD_INITIALIZE
    sel_srv3.initialize();
  #endif
  sel_srvind.allocate(1,2,1,nlenm,"sel_srvind");
  #ifndef NO_AD_INITIALIZE
    sel_srvind.initialize();
  #endif
  sel_srvnmfs.allocate(1,2,1,nlenm,"sel_srvnmfs");
  #ifndef NO_AD_INITIALIZE
    sel_srvnmfs.initialize();
  #endif
  sel_srv10ind.allocate(1,2,1,nlenm,"sel_srv10ind");
  #ifndef NO_AD_INITIALIZE
    sel_srv10ind.initialize();
  #endif
  sel_srv10nmfs.allocate(1,2,1,nlenm,"sel_srv10nmfs");
  #ifndef NO_AD_INITIALIZE
    sel_srv10nmfs.initialize();
  #endif
  avgsel_fish.allocate("avgsel_fish");
  #ifndef NO_AD_INITIALIZE
  avgsel_fish.initialize();
  #endif
  avgsel_fish_disc.allocate(1,2,"avgsel_fish_disc");
  #ifndef NO_AD_INITIALIZE
    avgsel_fish_disc.initialize();
  #endif
  avgsel_srv1.allocate(1,2,"avgsel_srv1");
  #ifndef NO_AD_INITIALIZE
    avgsel_srv1.initialize();
  #endif
  popn_lmale.allocate(styr,endyr+Nproj,"popn_lmale");
  #ifndef NO_AD_INITIALIZE
    popn_lmale.initialize();
  #endif
  popn_disc.allocate(1,2,styr,endyr+Nproj,"popn_disc");
  #ifndef NO_AD_INITIALIZE
    popn_disc.initialize();
  #endif
  totn_srv1.allocate(1,2,styr,endyr+Nproj,"totn_srv1");
  #ifndef NO_AD_INITIALIZE
    totn_srv1.initialize();
  #endif
  totn_srv2.allocate(1,2,1,2,"totn_srv2");
  #ifndef NO_AD_INITIALIZE
    totn_srv2.initialize();
  #endif
  totn_srv10.allocate(1,2,1,2,"totn_srv10");
  #ifndef NO_AD_INITIALIZE
    totn_srv10.initialize();
  #endif
  totn_trawl.allocate(1,2,styr,endyr+Nproj,"totn_trawl");
  #ifndef NO_AD_INITIALIZE
    totn_trawl.initialize();
  #endif
  explbiom.allocate(styr,endyr+Nproj,"explbiom");
  #ifndef NO_AD_INITIALIZE
    explbiom.initialize();
  #endif
  pred_bio.allocate(styr,endyr+Nproj,"pred_bio");
  #ifndef NO_AD_INITIALIZE
    pred_bio.initialize();
  #endif
  fspbio.allocate(styr,endyr+Nproj,"fspbio");
  #ifndef NO_AD_INITIALIZE
    fspbio.initialize();
  #endif
  mspbio.allocate(styr,endyr+Nproj,"mspbio");
  #ifndef NO_AD_INITIALIZE
    mspbio.initialize();
  #endif
  pred_srv1.allocate(1,2,styr,endyr+Nproj,1,nlenm,"pred_srv1");
  #ifndef NO_AD_INITIALIZE
    pred_srv1.initialize();
  #endif
  pred_p_fish.allocate(1,2,styr,endyr+Nproj,1,nlenm,"pred_p_fish");
  #ifndef NO_AD_INITIALIZE
    pred_p_fish.initialize();
  #endif
  pred_p_fish_fit.allocate(1,2,styr,endyr+Nproj,1,nlenm,"pred_p_fish_fit");
  #ifndef NO_AD_INITIALIZE
    pred_p_fish_fit.initialize();
  #endif
  pred_p_fish_discm.allocate(1,2,styr,endyr+Nproj,1,nlenm,"pred_p_fish_discm");
  #ifndef NO_AD_INITIALIZE
    pred_p_fish_discm.initialize();
  #endif
  pred_p_fish_discf.allocate(styr,endyr+Nproj,1,nlenm,"pred_p_fish_discf");
  #ifndef NO_AD_INITIALIZE
    pred_p_fish_discf.initialize();
  #endif
  pred_p_trawl.allocate(1,2,styr,endyr+Nproj,1,nlenm,"pred_p_trawl");
  #ifndef NO_AD_INITIALIZE
    pred_p_trawl.initialize();
  #endif
  pred_p_srv1_len.allocate(1,2,1,2,styr,endyr+Nproj,1,nlenm,"pred_p_srv1_len");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv1_len.initialize();
  #endif
  pred_p_srv1_len_new.allocate(1,2,1,2,styr,endyr+Nproj,1,nlenm,"pred_p_srv1_len_new");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv1_len_new.initialize();
  #endif
  pred_p_srv1_len_old.allocate(1,2,1,2,styr,endyr+Nproj,1,nlenm,"pred_p_srv1_len_old");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv1_len_old.initialize();
  #endif
  pred_p_srv2_len.allocate(1,2,1,2,1,nlenm,"pred_p_srv2_len");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv2_len.initialize();
  #endif
  pred_p_srv2_len_ind.allocate(1,2,1,nlenm,"pred_p_srv2_len_ind");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv2_len_ind.initialize();
  #endif
  pred_p_srv2_len_nmfs.allocate(1,2,1,nlenm,"pred_p_srv2_len_nmfs");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv2_len_nmfs.initialize();
  #endif
  pred_p_srv10_len_ind.allocate(1,2,1,nlenm,"pred_p_srv10_len_ind");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv10_len_ind.initialize();
  #endif
  pred_p_srv10_len_nmfs.allocate(1,2,1,nlenm,"pred_p_srv10_len_nmfs");
  #ifndef NO_AD_INITIALIZE
    pred_p_srv10_len_nmfs.initialize();
  #endif
  pred_catch.allocate(styr,endyr+Nproj,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  pred_catch_ret.allocate(styr,endyr+Nproj,"pred_catch_ret");
  #ifndef NO_AD_INITIALIZE
    pred_catch_ret.initialize();
  #endif
  pred_catch_disc.allocate(1,2,styr,endyr+Nproj,"pred_catch_disc");
  #ifndef NO_AD_INITIALIZE
    pred_catch_disc.initialize();
  #endif
  pred_catch_trawl.allocate(styr,endyr+Nproj,"pred_catch_trawl");
  #ifndef NO_AD_INITIALIZE
    pred_catch_trawl.initialize();
  #endif
  catch_male_ret.allocate(styr,endyr+Nproj,1,nlenm,"catch_male_ret");
  #ifndef NO_AD_INITIALIZE
    catch_male_ret.initialize();
  #endif
  catch_lmale.allocate(styr,endyr+Nproj,1,nlenm,"catch_lmale");
  #ifndef NO_AD_INITIALIZE
    catch_lmale.initialize();
  #endif
  catch_discp.allocate(1,2,styr,endyr+Nproj,1,nlenm,"catch_discp");
  #ifndef NO_AD_INITIALIZE
    catch_discp.initialize();
  #endif
  natlength.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength");
  #ifndef NO_AD_INITIALIZE
    natlength.initialize();
  #endif
  natlength_iold.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength_iold");
  #ifndef NO_AD_INITIALIZE
    natlength_iold.initialize();
  #endif
  natlength_inew.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength_inew");
  #ifndef NO_AD_INITIALIZE
    natlength_inew.initialize();
  #endif
  natlength_mold.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength_mold");
  #ifndef NO_AD_INITIALIZE
    natlength_mold.initialize();
  #endif
  natlength_mnew.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength_mnew");
  #ifndef NO_AD_INITIALIZE
    natlength_mnew.initialize();
  #endif
  natlength_old.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength_old");
  #ifndef NO_AD_INITIALIZE
    natlength_old.initialize();
  #endif
  natlength_new.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength_new");
  #ifndef NO_AD_INITIALIZE
    natlength_new.initialize();
  #endif
  natlength_i.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength_i");
  #ifndef NO_AD_INITIALIZE
    natlength_i.initialize();
  #endif
  natlength_mat.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natlength_mat");
  #ifndef NO_AD_INITIALIZE
    natlength_mat.initialize();
  #endif
  len_len.allocate(1,2,1,nlenm,1,nlenm,"len_len");
  #ifndef NO_AD_INITIALIZE
    len_len.initialize();
  #endif
  moltp.allocate(1,2,1,nlenm,"moltp");
  #ifndef NO_AD_INITIALIZE
    moltp.initialize();
  #endif
  moltp_mat.allocate(1,2,1,nlenm,"moltp_mat");
  #ifndef NO_AD_INITIALIZE
    moltp_mat.initialize();
  #endif
  Ftot.allocate(1,2,1,2,styr,endyr+Nproj,1,nlenm,"Ftot");
  #ifndef NO_AD_INITIALIZE
    Ftot.initialize();
  #endif
  F.allocate(1,2,styr,endyr+Nproj,1,nlenm,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  F_ret.allocate(1,2,styr,endyr+Nproj,1,nlenm,"F_ret");
  #ifndef NO_AD_INITIALIZE
    F_ret.initialize();
  #endif
  Fmat.allocate(1,2,styr,endyr+Nproj,1,nlenm,"Fmat");
  #ifndef NO_AD_INITIALIZE
    Fmat.initialize();
  #endif
  Fmat_ret.allocate(1,2,styr,endyr+Nproj,1,nlenm,"Fmat_ret");
  #ifndef NO_AD_INITIALIZE
    Fmat_ret.initialize();
  #endif
  Fimat.allocate(1,2,styr,endyr+Nproj,1,nlenm,"Fimat");
  #ifndef NO_AD_INITIALIZE
    Fimat.initialize();
  #endif
  Fimat_ret.allocate(1,2,styr,endyr+Nproj,1,nlenm,"Fimat_ret");
  #ifndef NO_AD_INITIALIZE
    Fimat_ret.initialize();
  #endif
  Fdiscm.allocate(1,2,styr,endyr+Nproj,1,nlenm,"Fdiscm");
  #ifndef NO_AD_INITIALIZE
    Fdiscm.initialize();
  #endif
  Fdiscf.allocate(styr,endyr+Nproj,1,nlenm,"Fdiscf");
  #ifndef NO_AD_INITIALIZE
    Fdiscf.initialize();
  #endif
  Fdisct.allocate(1,2,styr,endyr+Nproj,1,nlenm,"Fdisct");
  #ifndef NO_AD_INITIALIZE
    Fdisct.initialize();
  #endif
  S.allocate(1,2,1,2,styr,endyr+Nproj,1,nlenm,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  Simat.allocate(1,2,1,2,styr,endyr+Nproj,1,nlenm,"Simat");
  #ifndef NO_AD_INITIALIZE
    Simat.initialize();
  #endif
  Smat.allocate(1,2,1,2,styr,endyr+Nproj,1,nlenm,"Smat");
  #ifndef NO_AD_INITIALIZE
    Smat.initialize();
  #endif
  fmort.allocate(styr,endyr+Nproj,"fmort");
  #ifndef NO_AD_INITIALIZE
    fmort.initialize();
  #endif
  fmortret.allocate(styr,endyr+Nproj,"fmortret");
  #ifndef NO_AD_INITIALIZE
    fmortret.initialize();
  #endif
  fmortdf.allocate(styr,endyr+Nproj,"fmortdf");
  #ifndef NO_AD_INITIALIZE
    fmortdf.initialize();
  #endif
  fmortt.allocate(styr,endyr+Nproj,"fmortt");
  #ifndef NO_AD_INITIALIZE
    fmortt.initialize();
  #endif
  fmort_disc.allocate(styr,endyr+Nproj,"fmort_disc");
  #ifndef NO_AD_INITIALIZE
    fmort_disc.initialize();
  #endif
  rbar.allocate("rbar");
  #ifndef NO_AD_INITIALIZE
  rbar.initialize();
  #endif
  surv.allocate(1,2,"surv");
  #ifndef NO_AD_INITIALIZE
    surv.initialize();
  #endif
  offset.allocate(1,5,"offset");
  #ifndef NO_AD_INITIALIZE
    offset.initialize();
  #endif
  matest.allocate(1,nlenm,"matest");
  #ifndef NO_AD_INITIALIZE
    matest.initialize();
  #endif
  matestf.allocate(1,nlenm,"matestf");
  #ifndef NO_AD_INITIALIZE
    matestf.initialize();
  #endif
  rec_like.allocate("rec_like");
  #ifndef NO_AD_INITIALIZE
  rec_like.initialize();
  #endif
  catch_like1.allocate("catch_like1");
  #ifndef NO_AD_INITIALIZE
  catch_like1.initialize();
  #endif
  catch_like2.allocate("catch_like2");
  #ifndef NO_AD_INITIALIZE
  catch_like2.initialize();
  #endif
  catch_liket.allocate("catch_liket");
  #ifndef NO_AD_INITIALIZE
  catch_liket.initialize();
  #endif
  catch_likef.allocate("catch_likef");
  #ifndef NO_AD_INITIALIZE
  catch_likef.initialize();
  #endif
  growth_like.allocate("growth_like");
  #ifndef NO_AD_INITIALIZE
  growth_like.initialize();
  #endif
  len_likeyr.allocate(1,2,1,2,1,2,1,nobs_srv1_length,"len_likeyr");
  #ifndef NO_AD_INITIALIZE
    len_likeyr.initialize();
  #endif
  len_like.allocate(1,7,"len_like");
  #ifndef NO_AD_INITIALIZE
    len_like.initialize();
  #endif
  sel_like.allocate(1,4,"sel_like");
  #ifndef NO_AD_INITIALIZE
    sel_like.initialize();
  #endif
  sel_like_50m.allocate("sel_like_50m");
  #ifndef NO_AD_INITIALIZE
  sel_like_50m.initialize();
  #endif
  like_natm.allocate("like_natm");
  #ifndef NO_AD_INITIALIZE
  like_natm.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  surv_like.allocate("surv_like");
  #ifndef NO_AD_INITIALIZE
  surv_like.initialize();
  #endif
  surv2_like.allocate("surv2_like");
  #ifndef NO_AD_INITIALIZE
  surv2_like.initialize();
  #endif
  surv3_like.allocate("surv3_like");
  #ifndef NO_AD_INITIALIZE
  surv3_like.initialize();
  #endif
  surv10_like.allocate("surv10_like");
  #ifndef NO_AD_INITIALIZE
  surv10_like.initialize();
  #endif
  sexr_like.allocate("sexr_like");
  #ifndef NO_AD_INITIALIZE
  sexr_like.initialize();
  #endif
  like_initsmo.allocate("like_initsmo");
  #ifndef NO_AD_INITIALIZE
  like_initsmo.initialize();
  #endif
  like_q.allocate("like_q");
  #ifndef NO_AD_INITIALIZE
  like_q.initialize();
  #endif
  like_linff.allocate("like_linff");
  #ifndef NO_AD_INITIALIZE
  like_linff.initialize();
  #endif
  like_linfm.allocate("like_linfm");
  #ifndef NO_AD_INITIALIZE
  like_linfm.initialize();
  #endif
  like_growthkf.allocate("like_growthkf");
  #ifndef NO_AD_INITIALIZE
  like_growthkf.initialize();
  #endif
  like_growthkm.allocate("like_growthkm");
  #ifndef NO_AD_INITIALIZE
  like_growthkm.initialize();
  #endif
  like_af.allocate("like_af");
  #ifndef NO_AD_INITIALIZE
  like_af.initialize();
  #endif
  like_am.allocate("like_am");
  #ifndef NO_AD_INITIALIZE
  like_am.initialize();
  #endif
  like_bf.allocate("like_bf");
  #ifndef NO_AD_INITIALIZE
  like_bf.initialize();
  #endif
  like_bf1.allocate("like_bf1");
  #ifndef NO_AD_INITIALIZE
  like_bf1.initialize();
  #endif
  like_bm.allocate("like_bm");
  #ifndef NO_AD_INITIALIZE
  like_bm.initialize();
  #endif
  like_a1.allocate("like_a1");
  #ifndef NO_AD_INITIALIZE
  like_a1.initialize();
  #endif
  like_b1.allocate("like_b1");
  #ifndef NO_AD_INITIALIZE
  like_b1.initialize();
  #endif
  like_meetpt.allocate("like_meetpt");
  #ifndef NO_AD_INITIALIZE
  like_meetpt.initialize();
  #endif
  like_srvsel.allocate("like_srvsel");
  #ifndef NO_AD_INITIALIZE
  like_srvsel.initialize();
  #endif
  selsmo_like.allocate("selsmo_like");
  #ifndef NO_AD_INITIALIZE
  selsmo_like.initialize();
  #endif
  discf_like.allocate("discf_like");
  #ifndef NO_AD_INITIALIZE
  discf_like.initialize();
  #endif
  largemale_like.allocate("largemale_like");
  #ifndef NO_AD_INITIALIZE
  largemale_like.initialize();
  #endif
  like_mat.allocate("like_mat");
  #ifndef NO_AD_INITIALIZE
  like_mat.initialize();
  #endif
  like_initnum.allocate("like_initnum");
  #ifndef NO_AD_INITIALIZE
  like_initnum.initialize();
  #endif
  len_like10_ind.allocate("len_like10_ind");
  #ifndef NO_AD_INITIALIZE
  len_like10_ind.initialize();
  #endif
  len_like10_nmfs.allocate("len_like10_nmfs");
  #ifndef NO_AD_INITIALIZE
  len_like10_nmfs.initialize();
  #endif
  surv10nmfs_like.allocate("surv10nmfs_like");
  #ifndef NO_AD_INITIALIZE
  surv10nmfs_like.initialize();
  #endif
  fspbios.allocate(styr,endyr+Nproj,"fspbios");
  mspbios.allocate(styr,endyr+Nproj,"mspbios");
  legal_malesd.allocate(styr,endyr+Nproj,"legal_malesd");
  recf_sd.allocate(styr,endyr+Nproj-1,"recf_sd");
  recm_sd.allocate(styr,endyr+Nproj-1,"recm_sd");
  depletion.allocate("depletion");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  tmp.allocate("tmp");
  #ifndef NO_AD_INITIALIZE
  tmp.initialize();
  #endif
  tmpn.allocate(1,nlenm,"tmpn");
  #ifndef NO_AD_INITIALIZE
    tmpn.initialize();
  #endif
  pred_sexr.allocate(styr,endyr+Nproj,"pred_sexr");
  #ifndef NO_AD_INITIALIZE
    pred_sexr.initialize();
  #endif
  preds_sexr.allocate(styr,endyr+Nproj,"preds_sexr");
  #ifndef NO_AD_INITIALIZE
    preds_sexr.initialize();
  #endif
  predpop_sexr.allocate(styr,endyr+Nproj,"predpop_sexr");
  #ifndef NO_AD_INITIALIZE
    predpop_sexr.initialize();
  #endif
  maxsel_fish.allocate("maxsel_fish");
  #ifndef NO_AD_INITIALIZE
  maxsel_fish.initialize();
  #endif
  maxsel_srv1.allocate("maxsel_srv1");
  #ifndef NO_AD_INITIALIZE
  maxsel_srv1.initialize();
  #endif
  maxseld.allocate("maxseld");
  #ifndef NO_AD_INITIALIZE
  maxseld.initialize();
  #endif
  mean_length.allocate(1,2,1,nlenm,"mean_length");
  #ifndef NO_AD_INITIALIZE
    mean_length.initialize();
  #endif
  sd_mean_length.allocate(1,2,1,nlenm,"sd_mean_length");
  #ifndef NO_AD_INITIALIZE
    sd_mean_length.initialize();
  #endif
  sum_len.allocate(1,2,1,nlenm,"sum_len");
  #ifndef NO_AD_INITIALIZE
    sum_len.initialize();
  #endif
  tmp1.allocate("tmp1");
  #ifndef NO_AD_INITIALIZE
  tmp1.initialize();
  #endif
  tmp2.allocate("tmp2");
  #ifndef NO_AD_INITIALIZE
  tmp2.initialize();
  #endif
  tmp3.allocate("tmp3");
  #ifndef NO_AD_INITIALIZE
  tmp3.initialize();
  #endif
  sumfish.allocate(1,nobs_fish,"sumfish");
  #ifndef NO_AD_INITIALIZE
    sumfish.initialize();
  #endif
  sumfishret.allocate(1,nobs_fish,"sumfishret");
  #ifndef NO_AD_INITIALIZE
    sumfishret.initialize();
  #endif
  sumfishtot.allocate(1,nobs_fish_discm,"sumfishtot");
  #ifndef NO_AD_INITIALIZE
    sumfishtot.initialize();
  #endif
  sumfishdiscm.allocate(1,nobs_fish_discm,"sumfishdiscm");
  #ifndef NO_AD_INITIALIZE
    sumfishdiscm.initialize();
  #endif
  sumfishdiscf.allocate(1,nobs_fish_discf,"sumfishdiscf");
  #ifndef NO_AD_INITIALIZE
    sumfishdiscf.initialize();
  #endif
  sumtrawl.allocate(1,nobs_trawl,"sumtrawl");
  #ifndef NO_AD_INITIALIZE
    sumtrawl.initialize();
  #endif
  sumsrv.allocate(1,nobs_srv1_length,"sumsrv");
  #ifndef NO_AD_INITIALIZE
    sumsrv.initialize();
  #endif
  obs_p_srv2_len.allocate(1,2,1,2,1,2,1,nlenm,"obs_p_srv2_len");
  #ifndef NO_AD_INITIALIZE
    obs_p_srv2_len.initialize();
  #endif
  obs_p_srv10_len.allocate(1,2,1,2,1,2,1,nlenm,"obs_p_srv10_len");
  #ifndef NO_AD_INITIALIZE
    obs_p_srv10_len.initialize();
  #endif
  obs_p_srv1_len.allocate(1,2,1,2,1,2,1,nobs_srv1_length,1,nlenm,"obs_p_srv1_len");
  #ifndef NO_AD_INITIALIZE
    obs_p_srv1_len.initialize();
  #endif
  obs_p_srv1_len1.allocate(1,2,1,2,1,2,1,nobs_srv1_length,1,nlenm,"obs_p_srv1_len1");
  #ifndef NO_AD_INITIALIZE
    obs_p_srv1_len1.initialize();
  #endif
  obs_p_srv1_lenc.allocate(1,2,1,nobs_srv1_length,1,nlenm,"obs_p_srv1_lenc");
  #ifndef NO_AD_INITIALIZE
    obs_p_srv1_lenc.initialize();
  #endif
  obs_p_fish.allocate(1,2,1,2,1,nobs_fish,1,nlenm,"obs_p_fish");
  #ifndef NO_AD_INITIALIZE
    obs_p_fish.initialize();
  #endif
  obs_sexr.allocate(1,nobs_fish,"obs_sexr");
  #ifndef NO_AD_INITIALIZE
    obs_sexr.initialize();
  #endif
  obs_sexr_srv1_l.allocate(1,nobs_srv1_length,"obs_sexr_srv1_l");
  #ifndef NO_AD_INITIALIZE
    obs_sexr_srv1_l.initialize();
  #endif
  obs_p_fish_ret.allocate(1,2,1,nobs_fish,1,nlenm,"obs_p_fish_ret");
  #ifndef NO_AD_INITIALIZE
    obs_p_fish_ret.initialize();
  #endif
  obs_p_fish_tot.allocate(1,2,1,nobs_fish,1,nlenm,"obs_p_fish_tot");
  #ifndef NO_AD_INITIALIZE
    obs_p_fish_tot.initialize();
  #endif
  obs_p_fish_discm.allocate(1,2,1,nobs_fish_discm,1,nlenm,"obs_p_fish_discm");
  #ifndef NO_AD_INITIALIZE
    obs_p_fish_discm.initialize();
  #endif
  obs_p_fish_discf.allocate(1,nobs_fish_discf,1,nlenm,"obs_p_fish_discf");
  #ifndef NO_AD_INITIALIZE
    obs_p_fish_discf.initialize();
  #endif
  obs_p_trawl.allocate(1,2,1,nobs_trawl,1,nlenm,"obs_p_trawl");
  #ifndef NO_AD_INITIALIZE
    obs_p_trawl.initialize();
  #endif
  pf.allocate(1,nlenm,"pf");
  #ifndef NO_AD_INITIALIZE
    pf.initialize();
  #endif
  pb.allocate(1,nlenm,"pb");
  #ifndef NO_AD_INITIALIZE
    pb.initialize();
  #endif
  catch_tot.allocate(styr,endyr+Nproj,"catch_tot");
  #ifndef NO_AD_INITIALIZE
    catch_tot.initialize();
  #endif
  legal_males.allocate(styr,endyr+Nproj,"legal_males");
  #ifndef NO_AD_INITIALIZE
    legal_males.initialize();
  #endif
  legal_srv_males.allocate(styr,endyr+Nproj,"legal_srv_males");
  #ifndef NO_AD_INITIALIZE
    legal_srv_males.initialize();
  #endif
  popn.allocate(styr,endyr+Nproj,"popn");
  #ifndef NO_AD_INITIALIZE
    popn.initialize();
  #endif
  popn_fit.allocate(styr,endyr+Nproj,"popn_fit");
  #ifndef NO_AD_INITIALIZE
    popn_fit.initialize();
  #endif
  obs_srv1_num.allocate(1,2,styr,endyr+Nproj,1,nlenm,"obs_srv1_num");
  #ifndef NO_AD_INITIALIZE
    obs_srv1_num.initialize();
  #endif
  tmpo.allocate(1,2,styr,endyr+Nproj,"tmpo");
  #ifndef NO_AD_INITIALIZE
    tmpo.initialize();
  #endif
  tmpp.allocate(1,2,styr,endyr+Nproj,"tmpp");
  #ifndef NO_AD_INITIALIZE
    tmpp.initialize();
  #endif
  obs_srv1_bioms.allocate(1,2,styr,endyr+Nproj,"obs_srv1_bioms");
  #ifndef NO_AD_INITIALIZE
    obs_srv1_bioms.initialize();
  #endif
  obs_srv1_biom.allocate(styr,endyr+Nproj,"obs_srv1_biom");
  #ifndef NO_AD_INITIALIZE
    obs_srv1_biom.initialize();
  #endif
  obs_srv1_spbiom.allocate(1,2,styr,endyr+Nproj,"obs_srv1_spbiom");
  #ifndef NO_AD_INITIALIZE
    obs_srv1_spbiom.initialize();
  #endif
  obs_srv2_spbiom.allocate(1,2,1,2,"obs_srv2_spbiom");
  #ifndef NO_AD_INITIALIZE
    obs_srv2_spbiom.initialize();
  #endif
  obs_srv10_spbiom.allocate(1,2,1,2,"obs_srv10_spbiom");
  #ifndef NO_AD_INITIALIZE
    obs_srv10_spbiom.initialize();
  #endif
  obs_srv1_spnum.allocate(1,2,1,2,styr,endyr+Nproj,"obs_srv1_spnum");
  #ifndef NO_AD_INITIALIZE
    obs_srv1_spnum.initialize();
  #endif
  pred_srv1_biom.allocate(styr,endyr+Nproj,"pred_srv1_biom");
  #ifndef NO_AD_INITIALIZE
    pred_srv1_biom.initialize();
  #endif
  pred_srv1_bioms.allocate(1,2,styr,endyr+Nproj,"pred_srv1_bioms");
  #ifndef NO_AD_INITIALIZE
    pred_srv1_bioms.initialize();
  #endif
  obs_srv1t.allocate(styr,endyr+Nproj,"obs_srv1t");
  #ifndef NO_AD_INITIALIZE
    obs_srv1t.initialize();
  #endif
  cv_srv1.allocate(1,2,styr,endyr+Nproj,"cv_srv1");
  #ifndef NO_AD_INITIALIZE
    cv_srv1.initialize();
  #endif
  cv_srv1_nowt.allocate(1,2,styr,endyr+Nproj,"cv_srv1_nowt");
  #ifndef NO_AD_INITIALIZE
    cv_srv1_nowt.initialize();
  #endif
  obs_catchdm_biom.allocate(styr,endyr+Nproj,"obs_catchdm_biom");
  #ifndef NO_AD_INITIALIZE
    obs_catchdm_biom.initialize();
  #endif
  obs_catchdf_biom.allocate(styr,endyr+Nproj,"obs_catchdf_biom");
  #ifndef NO_AD_INITIALIZE
    obs_catchdf_biom.initialize();
  #endif
  obs_catcht_biom.allocate(styr,endyr+Nproj,"obs_catcht_biom");
  #ifndef NO_AD_INITIALIZE
    obs_catcht_biom.initialize();
  #endif
  obs_catchtot_biom.allocate(styr,endyr+Nproj,"obs_catchtot_biom");
  #ifndef NO_AD_INITIALIZE
    obs_catchtot_biom.initialize();
  #endif
  avgp.allocate(1,nlenm,"avgp");
  #ifndef NO_AD_INITIALIZE
    avgp.initialize();
  #endif
  avgpf.allocate(1,nlenm,"avgpf");
  #ifndef NO_AD_INITIALIZE
    avgpf.initialize();
  #endif
  avgpm.allocate(1,nlenm,"avgpm");
  #ifndef NO_AD_INITIALIZE
    avgpm.initialize();
  #endif
  biom_tmp.allocate(1,2,styr,endyr+Nproj,"biom_tmp");
  #ifndef NO_AD_INITIALIZE
    biom_tmp.initialize();
  #endif
  discardc.allocate(1,2,1,nlenm,"discardc");
  #ifndef NO_AD_INITIALIZE
    discardc.initialize();
  #endif
  rec_len.allocate(1,nlenm,"rec_len");
  #ifndef NO_AD_INITIALIZE
    rec_len.initialize();
  #endif
  avg_rec.allocate(1,2,"avg_rec");
  #ifndef NO_AD_INITIALIZE
    avg_rec.initialize();
  #endif
  obs_lmales.allocate(1,nobs_srv1_length,"obs_lmales");
  #ifndef NO_AD_INITIALIZE
    obs_lmales.initialize();
  #endif
  rec_dev.allocate(1,2,styr,endyr,"rec_dev");
  #ifndef NO_AD_INITIALIZE
    rec_dev.initialize();
  #endif
  avgsel.allocate(1,2,1,nlenm,"avgsel");
  #ifndef NO_AD_INITIALIZE
    avgsel.initialize();
  #endif
  avgsel_fit.allocate(1,2,1,nlenm,"avgsel_fit");
  #ifndef NO_AD_INITIALIZE
    avgsel_fit.initialize();
  #endif
  sel_avg_like.allocate("sel_avg_like");
  #ifndef NO_AD_INITIALIZE
  sel_avg_like.initialize();
  #endif
  effn_srv1.allocate(1,2,1,2,1,2,styr,endyr+Nproj,"effn_srv1");
  #ifndef NO_AD_INITIALIZE
    effn_srv1.initialize();
  #endif
  effn_fish_ret.allocate(1,2,styr,endyr+Nproj,"effn_fish_ret");
  #ifndef NO_AD_INITIALIZE
    effn_fish_ret.initialize();
  #endif
  effn_fish_tot.allocate(1,2,styr,endyr+Nproj,"effn_fish_tot");
  #ifndef NO_AD_INITIALIZE
    effn_fish_tot.initialize();
  #endif
  legal_males_bio.allocate(styr,endyr+Nproj,"legal_males_bio");
  #ifndef NO_AD_INITIALIZE
    legal_males_bio.initialize();
  #endif
  legal_srv_males_n.allocate(styr,endyr+Nproj,"legal_srv_males_n");
  #ifndef NO_AD_INITIALIZE
    legal_srv_males_n.initialize();
  #endif
  legal_srv_males_o.allocate(styr,endyr+Nproj,"legal_srv_males_o");
  #ifndef NO_AD_INITIALIZE
    legal_srv_males_o.initialize();
  #endif
  legal_srv_males_bio.allocate(styr,endyr+Nproj,"legal_srv_males_bio");
  #ifndef NO_AD_INITIALIZE
    legal_srv_males_bio.initialize();
  #endif
  obs_lmales_bio.allocate(1,nobs_srv1_length,"obs_lmales_bio");
  #ifndef NO_AD_INITIALIZE
    obs_lmales_bio.initialize();
  #endif
  natl_new_fishtime.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natl_new_fishtime");
  #ifndef NO_AD_INITIALIZE
    natl_new_fishtime.initialize();
  #endif
  natl_old_fishtime.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natl_old_fishtime");
  #ifndef NO_AD_INITIALIZE
    natl_old_fishtime.initialize();
  #endif
  natl_inew_fishtime.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natl_inew_fishtime");
  #ifndef NO_AD_INITIALIZE
    natl_inew_fishtime.initialize();
  #endif
  natl_iold_fishtime.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natl_iold_fishtime");
  #ifndef NO_AD_INITIALIZE
    natl_iold_fishtime.initialize();
  #endif
  natl_mnew_fishtime.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natl_mnew_fishtime");
  #ifndef NO_AD_INITIALIZE
    natl_mnew_fishtime.initialize();
  #endif
  natl_mold_fishtime.allocate(1,2,styr,endyr+Nproj,1,nlenm,"natl_mold_fishtime");
  #ifndef NO_AD_INITIALIZE
    natl_mold_fishtime.initialize();
  #endif
  popn_lmale_new.allocate(styr,endyr+Nproj,"popn_lmale_new");
  #ifndef NO_AD_INITIALIZE
    popn_lmale_new.initialize();
  #endif
  popn_lmale_old.allocate(styr,endyr+Nproj,"popn_lmale_old");
  #ifndef NO_AD_INITIALIZE
    popn_lmale_old.initialize();
  #endif
  catch_lmale_new.allocate(styr,endyr+Nproj,1,nlenm,"catch_lmale_new");
  #ifndef NO_AD_INITIALIZE
    catch_lmale_new.initialize();
  #endif
  catch_lmale_old.allocate(styr,endyr+Nproj,1,nlenm,"catch_lmale_old");
  #ifndef NO_AD_INITIALIZE
    catch_lmale_old.initialize();
  #endif
  catch_male_ret_new.allocate(styr,endyr+Nproj,1,nlenm,"catch_male_ret_new");
  #ifndef NO_AD_INITIALIZE
    catch_male_ret_new.initialize();
  #endif
  catch_male_ret_old.allocate(styr,endyr+Nproj,1,nlenm,"catch_male_ret_old");
  #ifndef NO_AD_INITIALIZE
    catch_male_ret_old.initialize();
  #endif
  popn_fit_new.allocate(styr,endyr+Nproj,"popn_fit_new");
  #ifndef NO_AD_INITIALIZE
    popn_fit_new.initialize();
  #endif
  popn_fit_old.allocate(styr,endyr+Nproj,"popn_fit_old");
  #ifndef NO_AD_INITIALIZE
    popn_fit_old.initialize();
  #endif
  popn_lmale_bio.allocate(styr,endyr+Nproj,"popn_lmale_bio");
  #ifndef NO_AD_INITIALIZE
    popn_lmale_bio.initialize();
  #endif
  like_mmat.allocate("like_mmat");
  #ifndef NO_AD_INITIALIZE
  like_mmat.initialize();
  #endif
  sumrecf.allocate("sumrecf");
  #ifndef NO_AD_INITIALIZE
  sumrecf.initialize();
  #endif
  sumrecm.allocate("sumrecm");
  #ifndef NO_AD_INITIALIZE
  sumrecm.initialize();
  #endif
  cv_rec.allocate(1,2,"cv_rec");
  #ifndef NO_AD_INITIALIZE
    cv_rec.initialize();
  #endif
  fspbio_srv1.allocate(styr,endyr+Nproj,"fspbio_srv1");
  #ifndef NO_AD_INITIALIZE
    fspbio_srv1.initialize();
  #endif
  mspbio_srv1.allocate(styr,endyr+Nproj,"mspbio_srv1");
  #ifndef NO_AD_INITIALIZE
    mspbio_srv1.initialize();
  #endif
  fspbio_srv2_ind.allocate("fspbio_srv2_ind");
  #ifndef NO_AD_INITIALIZE
  fspbio_srv2_ind.initialize();
  #endif
  mspbio_srv2_ind.allocate("mspbio_srv2_ind");
  #ifndef NO_AD_INITIALIZE
  mspbio_srv2_ind.initialize();
  #endif
  fspbio_srv2_nmfs.allocate("fspbio_srv2_nmfs");
  #ifndef NO_AD_INITIALIZE
  fspbio_srv2_nmfs.initialize();
  #endif
  mspbio_srv2_nmfs.allocate("mspbio_srv2_nmfs");
  #ifndef NO_AD_INITIALIZE
  mspbio_srv2_nmfs.initialize();
  #endif
  fspbio_srv10_ind.allocate("fspbio_srv10_ind");
  #ifndef NO_AD_INITIALIZE
  fspbio_srv10_ind.initialize();
  #endif
  mspbio_srv10_ind.allocate("mspbio_srv10_ind");
  #ifndef NO_AD_INITIALIZE
  mspbio_srv10_ind.initialize();
  #endif
  fspbio_srv10_nmfs.allocate("fspbio_srv10_nmfs");
  #ifndef NO_AD_INITIALIZE
  fspbio_srv10_nmfs.initialize();
  #endif
  mspbio_srv10_nmfs.allocate("mspbio_srv10_nmfs");
  #ifndef NO_AD_INITIALIZE
  mspbio_srv10_nmfs.initialize();
  #endif
  fspbio_srv1_num.allocate(1,2,styr,endyr+Nproj,"fspbio_srv1_num");
  #ifndef NO_AD_INITIALIZE
    fspbio_srv1_num.initialize();
  #endif
  mspbio_srv1_num.allocate(1,2,styr,endyr+Nproj,"mspbio_srv1_num");
  #ifndef NO_AD_INITIALIZE
    mspbio_srv1_num.initialize();
  #endif
  offset_srv.allocate(1,2,1,2,1,2,"offset_srv");
  #ifndef NO_AD_INITIALIZE
    offset_srv.initialize();
  #endif
  len_like_srv.allocate(1,2,1,2,1,2,"len_like_srv");
  #ifndef NO_AD_INITIALIZE
    len_like_srv.initialize();
  #endif
  growinc_90.allocate("growinc_90");
  #ifndef NO_AD_INITIALIZE
  growinc_90.initialize();
  #endif
  growinc_67.allocate("growinc_67");
  #ifndef NO_AD_INITIALIZE
  growinc_67.initialize();
  #endif
  hrate.allocate("hrate");
  #ifndef NO_AD_INITIALIZE
  hrate.initialize();
  #endif
  ghl.allocate("ghl");
  #ifndef NO_AD_INITIALIZE
  ghl.initialize();
  #endif
  ghl_number.allocate("ghl_number");
  #ifndef NO_AD_INITIALIZE
  ghl_number.initialize();
  #endif
  emspbio_matetime.allocate(styr,endyr+Nproj,"emspbio_matetime");
  #ifndef NO_AD_INITIALIZE
    emspbio_matetime.initialize();
  #endif
  efspbio_matetime.allocate(styr,endyr+Nproj,"efspbio_matetime");
  #ifndef NO_AD_INITIALIZE
    efspbio_matetime.initialize();
  #endif
  mspbio_matetime.allocate(styr,endyr+Nproj,"mspbio_matetime");
  #ifndef NO_AD_INITIALIZE
    mspbio_matetime.initialize();
  #endif
  mspbio_old_matetime.allocate(styr,endyr+Nproj,"mspbio_old_matetime");
  #ifndef NO_AD_INITIALIZE
    mspbio_old_matetime.initialize();
  #endif
  fspbio_new_matetime.allocate(styr,endyr+Nproj,"fspbio_new_matetime");
  #ifndef NO_AD_INITIALIZE
    fspbio_new_matetime.initialize();
  #endif
  efspbio_new_matetime.allocate(styr,endyr+Nproj,"efspbio_new_matetime");
  #ifndef NO_AD_INITIALIZE
    efspbio_new_matetime.initialize();
  #endif
  fspnum_new_matetime.allocate(styr,endyr+Nproj,"fspnum_new_matetime");
  #ifndef NO_AD_INITIALIZE
    fspnum_new_matetime.initialize();
  #endif
  fspbio_matetime.allocate(styr,endyr+Nproj,"fspbio_matetime");
  #ifndef NO_AD_INITIALIZE
    fspbio_matetime.initialize();
  #endif
  mspbio_fishtime.allocate(styr,endyr+Nproj,"mspbio_fishtime");
  #ifndef NO_AD_INITIALIZE
    mspbio_fishtime.initialize();
  #endif
  fspbio_fishtime.allocate(styr,endyr+Nproj,"fspbio_fishtime");
  #ifndef NO_AD_INITIALIZE
    fspbio_fishtime.initialize();
  #endif
  emspnum_old_matetime.allocate(styr,endyr+Nproj,"emspnum_old_matetime");
  #ifndef NO_AD_INITIALIZE
    emspnum_old_matetime.initialize();
  #endif
  mspnum_matetime.allocate(styr,endyr+Nproj,"mspnum_matetime");
  #ifndef NO_AD_INITIALIZE
    mspnum_matetime.initialize();
  #endif
  efspnum_matetime.allocate(styr,endyr+Nproj,"efspnum_matetime");
  #ifndef NO_AD_INITIALIZE
    efspnum_matetime.initialize();
  #endif
  pred_catch_gt101.allocate(styr,endyr+Nproj,"pred_catch_gt101");
  #ifndef NO_AD_INITIALIZE
    pred_catch_gt101.initialize();
  #endif
  pred_catch_no_gt101.allocate(styr,endyr+Nproj,"pred_catch_no_gt101");
  #ifndef NO_AD_INITIALIZE
    pred_catch_no_gt101.initialize();
  #endif
  num_males_gt101.allocate(styr,endyr+Nproj,"num_males_gt101");
  #ifndef NO_AD_INITIALIZE
    num_males_gt101.initialize();
  #endif
  bio_males_gt101.allocate(styr,endyr+Nproj,"bio_males_gt101");
  #ifndef NO_AD_INITIALIZE
    bio_males_gt101.initialize();
  #endif
  obs_tmp.allocate(styr,endyr+Nproj,"obs_tmp");
  #ifndef NO_AD_INITIALIZE
    obs_tmp.initialize();
  #endif
  catch_midpt.allocate(styr,endyr+Nproj,"catch_midpt");
  #ifndef NO_AD_INITIALIZE
    catch_midpt.initialize();
  #endif
  alpha.allocate(1,2,"alpha");
  #ifndef NO_AD_INITIALIZE
    alpha.initialize();
  #endif
  alpha_rec.allocate("alpha_rec");
  #ifndef NO_AD_INITIALIZE
  alpha_rec.initialize();
  #endif
  avg_beta.allocate("avg_beta");
  #ifndef NO_AD_INITIALIZE
  avg_beta.initialize();
  #endif
  devia.allocate("devia");
  #ifndef NO_AD_INITIALIZE
  devia.initialize();
  #endif
  recsum.allocate("recsum");
  #ifndef NO_AD_INITIALIZE
  recsum.initialize();
  #endif
  tmps.allocate(1,nlenm,"tmps");
  #ifndef NO_AD_INITIALIZE
    tmps.initialize();
  #endif
  tmpi.allocate("tmpi");
  #ifndef NO_AD_INITIALIZE
  tmpi.initialize();
  #endif
  catch_disc.allocate(1,2,styr,endyr+Nproj,"catch_disc");
  #ifndef NO_AD_INITIALIZE
    catch_disc.initialize();
  #endif
  old_mult.allocate("old_mult");
  #ifndef NO_AD_INITIALIZE
  old_mult.initialize();
  #endif
  maturity_est.allocate(1,2,1,nlenm,"maturity_est");
  #ifndef NO_AD_INITIALIZE
    maturity_est.initialize();
  #endif
  cpue_pred.allocate(styr,endyr+Nproj,"cpue_pred");
  #ifndef NO_AD_INITIALIZE
    cpue_pred.initialize();
  #endif
  cpue_like.allocate("cpue_like");
  #ifndef NO_AD_INITIALIZE
  cpue_like.initialize();
  #endif
  ftarget.allocate(styr,endyr+Nproj,"ftarget");
  #ifndef NO_AD_INITIALIZE
    ftarget.initialize();
  #endif
  pred_catch_target.allocate(styr,endyr+Nproj,"pred_catch_target");
  #ifndef NO_AD_INITIALIZE
    pred_catch_target.initialize();
  #endif
  ghl_like.allocate("ghl_like");
  #ifndef NO_AD_INITIALIZE
  ghl_like.initialize();
  #endif
  wt_lmlike2.allocate("wt_lmlike2");
  #ifndef NO_AD_INITIALIZE
  wt_lmlike2.initialize();
  #endif
  multi.allocate("multi");
  #ifndef NO_AD_INITIALIZE
  multi.initialize();
  #endif
  tmpmnew.allocate(1,nlenm,"tmpmnew");
  #ifndef NO_AD_INITIALIZE
    tmpmnew.initialize();
  #endif
  tmpinew.allocate(1,nlenm,"tmpinew");
  #ifndef NO_AD_INITIALIZE
    tmpinew.initialize();
  #endif
  tmpmold.allocate(1,nlenm,"tmpmold");
  #ifndef NO_AD_INITIALIZE
    tmpmold.initialize();
  #endif
  offset_srv2.allocate(1,2,"offset_srv2");
  #ifndef NO_AD_INITIALIZE
    offset_srv2.initialize();
  #endif
  offset_srv10.allocate(1,2,"offset_srv10");
  #ifndef NO_AD_INITIALIZE
    offset_srv10.initialize();
  #endif
  selsmo2_n.allocate(1,nlenm,"selsmo2_n");
  #ifndef NO_AD_INITIALIZE
    selsmo2_n.initialize();
  #endif
  selsmof2_n.allocate(1,nlenm,"selsmof2_n");
  #ifndef NO_AD_INITIALIZE
    selsmof2_n.initialize();
  #endif
  Fout.allocate(1,30,"Fout");
  #ifndef NO_AD_INITIALIZE
    Fout.initialize();
  #endif
  call_no.allocate("call_no");
  #ifndef NO_AD_INITIALIZE
  call_no.initialize();
  #endif
  M_matn.allocate(1,2,"M_matn");
  #ifndef NO_AD_INITIALIZE
    M_matn.initialize();
  #endif
  M_mato.allocate(1,2,"M_mato");
  #ifndef NO_AD_INITIALIZE
    M_mato.initialize();
  #endif
  M.allocate(1,2,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  tmpp1.allocate(1,nlenm,"tmpp1");
  #ifndef NO_AD_INITIALIZE
    tmpp1.initialize();
  #endif
  tmpp2.allocate(1,nlenm,"tmpp2");
  #ifndef NO_AD_INITIALIZE
    tmpp2.initialize();
  #endif
  tmpp3.allocate(1,nlenm,"tmpp3");
  #ifndef NO_AD_INITIALIZE
    tmpp3.initialize();
  #endif
  tmpp4.allocate(1,nlenm,"tmpp4");
  #ifndef NO_AD_INITIALIZE
    tmpp4.initialize();
  #endif
  len_like_immat.allocate("len_like_immat");
  #ifndef NO_AD_INITIALIZE
  len_like_immat.initialize();
  #endif
  len_like_mat.allocate("len_like_mat");
  #ifndef NO_AD_INITIALIZE
  len_like_mat.initialize();
  #endif
  len_like_ex.allocate("len_like_ex");
  #ifndef NO_AD_INITIALIZE
  len_like_ex.initialize();
  #endif
  sumobs.allocate("sumobs");
  #ifndef NO_AD_INITIALIZE
  sumobs.initialize();
  #endif
  sumobs2.allocate("sumobs2");
  #ifndef NO_AD_INITIALIZE
  sumobs2.initialize();
  #endif
  sumobs10.allocate("sumobs10");
  #ifndef NO_AD_INITIALIZE
  sumobs10.initialize();
  #endif
  sumobs102.allocate("sumobs102");
  #ifndef NO_AD_INITIALIZE
  sumobs102.initialize();
  #endif
  lp_q.allocate("lp_q");
  FutRec.allocate("FutRec");
  #ifndef NO_AD_INITIALIZE
  FutRec.initialize();
  #endif
  FutMort.allocate("FutMort");
  #ifndef NO_AD_INITIALIZE
  FutMort.initialize();
  #endif
  Bzero.allocate("Bzero");
  #ifndef NO_AD_INITIALIZE
  Bzero.initialize();
  #endif
  Target.allocate("Target");
  #ifndef NO_AD_INITIALIZE
  Target.initialize();
  #endif
  F35.allocate("F35");
  #ifndef NO_AD_INITIALIZE
  F35.initialize();
  #endif
  SBPRF35.allocate("SBPRF35");
  #ifndef NO_AD_INITIALIZE
  SBPRF35.initialize();
  #endif
  Ratio.allocate("Ratio");
  #ifndef NO_AD_INITIALIZE
  Ratio.initialize();
  #endif
  Btest.allocate("Btest");
  #ifndef NO_AD_INITIALIZE
  Btest.initialize();
  #endif
  OFL.allocate("OFL");
  #ifndef NO_AD_INITIALIZE
  OFL.initialize();
  #endif
  FOFL.allocate("FOFL");
  #ifndef NO_AD_INITIALIZE
  FOFL.initialize();
  #endif
  Bmsy.allocate("Bmsy");
  #ifndef NO_AD_INITIALIZE
  Bmsy.initialize();
  #endif
  mean_log_rec.allocate(1,2,"mean_log_rec");
  #ifndef NO_AD_INITIALIZE
    mean_log_rec.initialize();
  #endif
  Lbar.allocate(1,2,1,2,1,2,1,nobs_srv1_length,"Lbar");
  #ifndef NO_AD_INITIALIZE
    Lbar.initialize();
  #endif
  Lbar_hat_old.allocate(1,2,1,2,1,nobs_srv1_length,"Lbar_hat_old");
  #ifndef NO_AD_INITIALIZE
    Lbar_hat_old.initialize();
  #endif
  Lbar_hat_new.allocate(1,2,1,2,1,nobs_srv1_length,"Lbar_hat_new");
  #ifndef NO_AD_INITIALIZE
    Lbar_hat_new.initialize();
  #endif
  SE_Lbar_hat_old.allocate(1,2,1,2,1,nobs_srv1_length,"SE_Lbar_hat_old");
  #ifndef NO_AD_INITIALIZE
    SE_Lbar_hat_old.initialize();
  #endif
  SE_Lbar_hat_new.allocate(1,2,1,2,1,nobs_srv1_length,"SE_Lbar_hat_new");
  #ifndef NO_AD_INITIALIZE
    SE_Lbar_hat_new.initialize();
  #endif
  Francis_var_temp_new.allocate(1,2,1,2,1,nobs_srv1_length,"Francis_var_temp_new");
  #ifndef NO_AD_INITIALIZE
    Francis_var_temp_new.initialize();
  #endif
  Francis_var_temp_old.allocate(1,2,1,2,1,nobs_srv1_length,"Francis_var_temp_old");
  #ifndef NO_AD_INITIALIZE
    Francis_var_temp_old.initialize();
  #endif
  Francis_weight_m.allocate("Francis_weight_m");
  #ifndef NO_AD_INITIALIZE
  Francis_weight_m.initialize();
  #endif
  Francis_weight_f.allocate("Francis_weight_f");
  #ifndef NO_AD_INITIALIZE
  Francis_weight_f.initialize();
  #endif
  countFem.allocate("countFem");
  #ifndef NO_AD_INITIALIZE
  countFem.initialize();
  #endif
  countMal.allocate("countMal");
  #ifndef NO_AD_INITIALIZE
  countMal.initialize();
  #endif
  totalFem.allocate("totalFem");
  #ifndef NO_AD_INITIALIZE
  totalFem.initialize();
  #endif
  totalMal.allocate("totalMal");
  #ifndef NO_AD_INITIALIZE
  totalMal.initialize();
  #endif
  FemMeanVarTerm.allocate("FemMeanVarTerm");
  #ifndef NO_AD_INITIALIZE
  FemMeanVarTerm.initialize();
  #endif
  MaleMeanVarTerm.allocate("MaleMeanVarTerm");
  #ifndef NO_AD_INITIALIZE
  MaleMeanVarTerm.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
   int mat;
   int m;
   
  // Calculate weight at length
  for(i=1;i<=nlenm;i++)
  {
   wtm(i) = alpha_wt_m*pow(length_bins(i),beta_wt_m);
   wtf(1,i) = alpha_wt_imm_f*pow(length_bins(i),beta_wt_imm_f);
   wtf(2,i) = alpha_wt_mat_f*pow(length_bins(i),beta_wt_mat_f);
  }
    
  //change weight to tons
   wtm=wtm*.000001;
   wtf=wtf*.000001; 
   
   for(m=styr;m<=endyr;m++)
    catch_midpt(m) = catch_midptIn(m);
   for(m=endyr+1;m<=endyr+Nproj;m++)
    catch_midpt(m) = catch_midptIn(endyr);
   // mnatlen_styr(1) = log(10*(obs_p_srv1_lend(1,1,2,1)+obs_p_srv1_lend(2,1,2,1)+1e-02));
   // mnatlen_styr(2) = log(10*(obs_p_srv1_lend(1,2,2,1)+obs_p_srv1_lend(2,2,2,1)+1e-02));
   
   // for(j=1;j<=12;j++){
   // fnatlen_styr(1,j) = log(10*(obs_p_srv1_lend(1,1,1,1,j)+obs_p_srv1_lend(2,1,1,1,j)+1e-02));
   // fnatlen_styr(2,j) = log(10*(obs_p_srv1_lend(1,2,1,1,j)+obs_p_srv1_lend(2,2,1,1,j)+1e-02));
   // }
   cout<<1<<endl;
   if(maturity_switch>0){
     maturity_average(2) = maturity_logistic;
      for(i=styr;i<=endyr;i++){
         maturity(2,i)= maturity_logistic;
        }
     }
  sumfishdiscm.initialize();
  sumfishret.initialize();
  sumfishtot.initialize();
  sumfishdiscf.initialize();
  sumtrawl.initialize();
  //get total catch - sum up catches mult by assumed mortality
  for(i=styr; i<=endyr; i++)
  {
    catch_trawl(i)=catch_trawl(i)*m_trawl;
    catch_disc(1,i)=catch_odisc(1,i)*m_disc;
    catch_disc(2,i)=catch_odisc(2,i)*m_disc;      
    catch_tot(i)=catch_numbers(i)+catch_disc(2,i);
  }
  //sum each  fishery data
  for(i=1; i<=nobs_fish_discm; i++)
  { 
     sumfishdiscm(i)+=sum(obs_p_fish_discmd(1,i));
     sumfishdiscm(i)+=sum(obs_p_fish_discmd(2,i));
  } 
  
  for(i=1; i<=nobs_fish; i++){
    for(j=1; j<=2; j++){
      sumfishret(i)+=sum(obs_p_fish_retd(j,i));
    }
  } 
  
  for(i=1; i<=nobs_fish_discf; i++){
    sumfishdiscf(i)+=sum(obs_p_fish_discfd(i));
  }
  
  for(i=1; i<=nobs_trawl; i++){
    for(k=1;k<=2;k++){
      sumtrawl(i)+=sum(obs_p_trawld(k,i));
    }
  }
   cout<<2 <<endl;
  for(i=1; i<=nobs_srv1_length;i++)
  {
    obs_p_srv1_lenc(1,i)=obs_p_srv1_lend(1,1,1,i)+obs_p_srv1_lend(1,2,1,i)+obs_p_srv1_lend(2,1,1,i)+obs_p_srv1_lend(2,2,1,i);
    obs_p_srv1_lenc(2,i)=obs_p_srv1_lend(1,1,2,i)+obs_p_srv1_lend(1,2,2,i)+obs_p_srv1_lend(2,1,2,i)+obs_p_srv1_lend(2,2,2,i);
    obs_sexr_srv1_l(i)=sum(obs_p_srv1_lenc(1,i))/(sum(obs_p_srv1_lenc(1,i))+sum(obs_p_srv1_lenc(2,i)));
  }
  
  offset.initialize();
  offset_srv2.initialize(); 
   for(k=1; k<=2; k++)
    for (i=1; i <= nobs_fish_discm; i++)
    for (j=1; j<=nlenm; j++)
    {
      if(k<2){   //new shell
        obs_p_fish_discm(k,i,j)=((obs_p_fish_discmd(k,i,j)*fraction_new_error)/sumfishdiscm(i))*(catch_disc(2,yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)));
         }
        else{ 
          obs_p_fish_discm(k,i,j)=((obs_p_fish_discmd(k,i,j)+obs_p_fish_discmd(1,i,j)*(1.-fraction_new_error))/sumfishdiscm(i))*(catch_disc(2,yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)));
           }
    }
  for(i=1;i<=nobs_fish;i++)
   {
    obs_p_fish_ret(1,i) = obs_p_fish_retd(1,i)*fraction_new_error;
    obs_p_fish_ret(2,i) = obs_p_fish_retd(2,i)+(1.-fraction_new_error)*obs_p_fish_retd(1,i);
   }
    //make observations proportions by year      
         //fishery offset
   for (i=1; i <= nobs_fish; i++)
    {
        
      for (j=1; j<=nlenm; j++)
       {
        for(k=1; k<=2; k++)
           obs_p_fish_ret(k,i,j)=(obs_p_fish_ret(k,i,j)/sumfishret(i));
               offset(1)-=nsamples_fish(1,i)*(obs_p_fish_ret(1,i,j)+obs_p_fish_ret(2,i,j))*log(obs_p_fish_ret(1,i,j)+obs_p_fish_ret(2,i,j)+p_const);
        }
	}
   cout<<3<<endl;
   for(k=1; k<=2; k++)
    {
       for (i=1; i <= nobs_fish_discm; i++)
       {
         //make observations proportions by year      
         //fishery offset
           for (j=1; j<=nlenm; j++)
           {
               obs_p_fish_tot(k,i,j)=((obs_p_fish_ret(k,i+yrs_fish_discm(1)-styr,j)*catch_numbers(yrs_fish_discm(i)))/catch_tot(yrs_fish_discm(i)))+obs_p_fish_discm(k,i,j);
             //old and new shell together
             if (k<2)
                offset(2)-=nsamples_fish(k,i)*(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j) )*log(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j)+p_const);
           }
       }
    }
     cout<<3<<endl;
  //make observations proportions by year      
        //fishery offset
  for (i=1; i <= nobs_fish_discf; i++)
   for (j=1; j<=nlenm; j++)
    {
        obs_p_fish_discf(i,j)=((obs_p_fish_discfd(i,j))/sumfishdiscf(i));
        offset(3)-=nsamples_fish_discf(i)*obs_p_fish_discf(i,j)*log(obs_p_fish_discf(i,j)+p_const);
	}        
   cout<<3<<endl;
  for(k=1; k<=2; k++)
   for (i=1; i <= nobs_trawl; i++)
    {
      //make observations proportions by year      
      //fishery offset
      for (j=1; j<=nlenm; j++)
      {
        obs_p_trawl(k,i,j)=((obs_p_trawld(k,i,j))/sumtrawl(i));
        offset(5)-=nsamples_trawl(k,i)*obs_p_trawl(k,i,j)*log(obs_p_trawl(k,i,j)+p_const);
      }
      obs_catcht_biom(yrs_trawl(i))=(obs_p_trawl(1,i)*catch_trawl(yrs_trawl(i)))*wtf(2)+(obs_p_trawl(2,i)*catch_trawl(yrs_trawl(i)))*wtm;
    }
   cout<<4<<endl;
  sumsrv.initialize();
 for(ll=1; ll<=nobs_srv1_length; ll++)
  for(mat=1;mat<=2;mat++)
   for(k=1; k<=2; k++)
    for(j=1; j<=2; j++)
    sumsrv(ll)+=sum(obs_p_srv1_lend(mat,k,j,ll));
  for(mat=1; mat<=2; mat++) //maturity
   for(l=1; l<=2; l++) //shell condition
    for(k=1; k<=2;k++) //sex
    for (i=1; i <= nobs_srv1_length; i++)
    for (j=1; j<=nlenm; j++)
    {
 //only do new/old shell correction for mature crab
       if(mat<2){
                obs_p_srv1_len1(mat,l,k,i,j)=(obs_p_srv1_lend(mat,l,k,i,j))/sumsrv(i);
                }
        else{
        if(l<2){
                obs_p_srv1_len1(mat,l,k,i,j)=(obs_p_srv1_lend(mat,l,k,i,j)*fraction_new_error)/sumsrv(i);
         }
        else{
               obs_p_srv1_len1(mat,l,k,i,j)=(obs_p_srv1_lend(mat,l,k,i,j)+obs_p_srv1_lend(mat,1,k,i,j)*(1.-fraction_new_error))/sumsrv(i);
            }
         }
    }
   if(maturity_switch>0)
   {
    for(i=1; i <= nobs_srv1_length; i++)
    {
     //oldshell
        tmps = (obs_p_srv1_len1(1,2,2,i)+obs_p_srv1_len1(2,2,2,i));
       obs_p_srv1_len1(2,2,2,i) = elem_prod(maturity_old_average(2),tmps);
       obs_p_srv1_len1(1,2,2,i) = elem_prod(1.0-maturity_old_average(2),tmps);
    }
   }
   avgpf=0;
   avgpm=0;
   for(i=1; i <= nobs_trawl; i++)
   {
    avgpf+=obs_p_trawl(1,i);
    avgpm+=obs_p_trawl(2,i);
   }
   avgpf=avgpf/nobs_trawl;
   avgpm=avgpm/nobs_trawl;
   for(i=styr;i<=(yrs_trawl(1)-1);i++)
    obs_catcht_biom(i)=avgpf*catch_trawl(i)*wtf(2)+avgpm*catch_trawl(i)*wtm;
  obs_catchdm_biom.initialize();
  obs_catchdf_biom.initialize();
  avgp.initialize();
  for(i=1;i<=nobs_fish_discm;i++)
  {
   //obs_p_fish_discm are proportional to total catch (retained+discard)
    obs_catchdm_biom(yrs_fish_discm(i)) = catch_tot(yrs_fish_discm(i)) * (obs_p_fish_discm(1,i)+obs_p_fish_discm(2,i))*wtm;
    avgp += obs_p_fish_discm(1,i)+obs_p_fish_discm(2,i);
  }
  avgp=avgp/nobs_fish_discm;
  avgpf = mean(obs_p_fish_discf);
  // Historical approximation to observed catch biomass
  for(i=styr;i<=(yrs_fish_discm(1)-1);i++)
  {
   // cout<<obs_catchdm_biom<<endl;
    obs_catchdm_biom(i)=avgp*catch_tot(i)*wtm;
  } 
  obs_catchdm_biom(endyr)=avgp*catch_tot(endyr)*wtm;
  for(i=styr;i<=endyr;i++)
  {
    obs_catchtot_biom(i)=obs_catchdm_biom(i)+catch_ret(i);
    obs_catchdf_biom(i)=(avgpf*catch_disc(1,i))*wtf(2);
  }
  
  // Compute the moulting probabilities
  get_moltingp();
  // estimate growth function
  get_growth();
  // Set maturity
  get_maturity();
  what=1;
  // mean_log_rec(1) = mean_log_rec_f;
  // mean_log_rec(2) = mean_log_rec_f;
  
  // if(active(mean_log_rec_m))
	// mean_log_rec(2) = mean_log_rec_m; 
 
  
   cout<<"end prelim calcs"<<endl;
}

void model_parameters::userfunction(void)
{
  f =0.0;
  ofstream& post= *pad_post;
    M(1)=M_in(1)*Mmult_imat;    			 //natural mortality immature females then males
    M(2)=M_in(2)*Mmult_imat;    			 //natural mortality immature females then males
    M_matn(2)=M_matn_in(2)*Mmult;  //natural mortality mature new shell female/male
    M_mato(2)=M_mato_in(2)*Mmult;  //natural mortality mature old shell female/male
    M_matn(1)=M_matn_in(1)*Mmultf;  //natural mortality mature new shell female/male
    M_mato(1)=M_mato_in(1)*Mmultf;  //natural mortality mature old shell female/male
  lp_q=srv3_q;
  get_maturity();
  // cout<<"maturity"<<endl;
  if(what<2)
   {
    get_obs_survey();
    what=what+1;
   }
   get_moltingp();
  // cout<<"molting"<<endl;
   if(active(am) || active(af))
    get_growth();
   // cout<<" end of growth "<<endl;
   get_selectivity();
   // cout<<" end of sel "<<endl;
   get_mortality();
   // cout<<" end of mort "<<endl;
   get_numbers_at_len();
   // cout<<" end of numbers at len "<<endl;
   get_catch_at_len();
   // cout<<" end of catch at len "<<endl;
   evaluate_the_objective_function();
   if(mceval_phase())
   {
	   Find_F35();
	   Find_OFL();
	   post<<f<<" "<<Bmsy<<" "<< F35 << " " << FOFL << " " << OFL << endl;
	   post<<fspbio_srv1<<endl;
	   post<<mspbio_srv1<<endl;
	   post<<mspbio_matetime<<endl;
	   post<<fspbio_srv2_ind<<endl;
	   post<<mspbio_srv2_ind<<endl;
	   post<<fspbio_srv2_nmfs<<endl;
	   post<<mspbio_srv2_nmfs<<endl;
	   post<<fspbio_srv10_ind<<endl;
	   post<<mspbio_srv10_ind<<endl;
	   post<<fspbio_srv10_nmfs<<endl;
	   post<<mspbio_srv10_nmfs<<endl;
	   post<<pred_catch<<endl;
	   post<<pred_catch_disc<<endl;
	   post<<pred_catch_ret<<endl;
	   post<<pred_catch_trawl<<endl;
      //Don't need these now because they can be taken from the .psv file, but I'm scared to get rid of them yet
	   // post<<af<< " "<<am<<" "<<bf<<" "<<bm<<" "<<b1<< " " <<bf1 << " " << deltam << " " << deltaf<< " " <<st_gr<< " "<<growth_beta<<" ";
	   // post<<matest50f<< " "<<matestslpf<<" "<<matest50m<<" "<<matestslpm<<" "<<mateste<< " " <<matestfe << " " ;
	   // post<<moltp_af<< " "<<moltp_bf<<" "<<moltp_am<<" "<<moltp_bm<<" "<<moltp_ammat<< " " <<moltp_bmmat << " " ;  
 	   // post<<mean_log_rec_m<<" "<<mean_log_rec_f<<" "<<rec_devf<<" "<<rec_devf<<" "<<alpha1_rec<< " " <<beta_rec << " " ;  
	   // post<<mnatlen_styr<<" "<< fnatlen_styr << " ";
	   // post<<log_avg_fmort<< " "<<fmort_dev<<" "<<log_avg_fmortdf<<" "<<fmortdf_dev<<" "<<log_avg_fmortt<<" " <<fmortt_dev<< " " <<discard_mult << " ";
	   // post<<log_avg_sel50_mn<< " "<<log_sel50_dev_mn<<" "<<log_avg_sel50_mo<<" "<<log_sel50_dev_mo<<" "<<fish_slope_mn<<" " <<fish_slope_mo << " ";
 	   // post<<fish_fit_slope_mn<< " "<<fish_fit_sel50_mn<<" "<<fish_fit_slope_mo<<" "<<fish_fit_sel50_mo << " " ;  	
 	   // post<<fish_slope_mo2<< " "<<fish_sel50_mo2<<" "<<fish_slope_mn2<<" "<<fish_sel50_mn2 << " " ;  	
 	   // post<<fish_disc_slope_f<< " "<<fish_disc_sel50_f<<" "<<log_dev_50f<<" "<<fish_disc_slope_tf << " " <<fish_disc_sel50_tf<<" " ;
	   // post<<srv1_q<< " "<<srv1_sel95<<" "<<srv1_sel50<<" "<<srv2_q<<" "<<srv2_sel95<<" " <<srv2_sel50 << " "<<srv3_q<<" "<<srv3_sel95<<" " <<srv3_sel50 << " ";   
	   // post<<srv3_sel95_f<<" "<< srv3_sel50_f << " "<<" ";   
 	   // post<<srvind_q<< " "<<srvind_q_f<<" "<<srvind_sel95<<" "<<srvind_sel50<<" "<<srvind_sel95_f<< " " <<srvind_sel50_f << " " ;  
 	   // post<<srvnmfs_sel95<< " "<<srvnmfs_sel50<<" "<<srvnmfs_sel95_f<<" "<<srvnmfs_sel50_f << " " ;  	
 	   // post<<srv10ind_q<< " "<<srv10ind_q_f<<" "<<srv10ind_sel95<<" "<<srv10ind_sel50<<" "<<srv10ind_sel95_f<< " " <<srv10ind_sel50_f << " " ;  
 	   // post<<srv10nmfs_sel95<< " "<<srv10nmfs_sel50<<" "<<srv10nmfs_sel95_f<<" "<<srv10nmfs_sel50_f << " " ;  	
 	   // post<<selsmo10ind<< " "<<selsmo09ind<<" "<<Mmult_imat<<" "<<Mmult<<" "<<Mmultf << " " ;  
	   // post<<cpueq<<" "<< proprecn << " "<<endl;
   }
}

void model_parameters::get_maturity(void)
{
  ofstream& post= *pad_post;
  if(active(matestfe))
   {
     for(j=1;j<=nlenm;j++)
      {
	      if(j<(matest_n+1)){
	      matestf(j)=matestfe(j);
          }
          else{
	          matestf(j)=0.0;
            }
           maturity_est(1,j) = mfexp(matestf(j));
       }
   }
   else
    {    
	  maturity_est(1) = maturity_average(1);
    }
  if(active(mateste))
   {
     for(j=1;j<=nlenm;j++)
      {
	      if(j<(matestm_n+1)){
	      matest(j)=mateste(j);
          }
          else{
	          matest(j)=0.0;
            }
          maturity_est(2,j) = mfexp(matest(j));
      }
    }
    else{    
       maturity_est(2) = maturity_average(2);
    }
}

void model_parameters::get_catch_tot(void)
{
  ofstream& post= *pad_post;
  for(i=styr;i<=endyr;i++)
  {
    catch_disc(1,i)=catch_odisc(1,i)*discard_mult*m_disc;
    catch_disc(2,i)=catch_odisc(2,i)*discard_mult*m_disc;      
    catch_tot(i)=catch_numbers(i)+catch_disc(2,i);
  }
   for(k=1; k<=2; k++)
    for (i=1; i <= nobs_fish_discm; i++)
    for (j=1; j<=nlenm; j++)
    {
      if(k<2)
	   {
        obs_p_fish_discm(k,i,j)=((obs_p_fish_discmd(k,i,j))/sumfishdiscm(i))*(catch_disc(2,yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)));
        }
        else
		{
          obs_p_fish_discm(k,i,j)=((obs_p_fish_discmd(k,i,j)*old_mult)/sumfishdiscm(i))*(catch_disc(2,yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)));
        }
    }
   offset(2)=0.0;
  //make observations proportions by year; fishery offset
   for(k=1; k<=2; k++)
    for (i=1; i <= nobs_fish_discm; i++)
    for (j=1; j<=nlenm; j++)
    {
        obs_p_fish_tot(k,i,j)=(obs_p_fish_retd(k,i+yrs_fish_discm(1)-styr,j)/sumfishret(i+yrs_fish_discm(1)-styr))*(catch_numbers(yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)))+obs_p_fish_discm(k,i,j);
   //old and new shell together
       if (k<2){
         offset(2)-=nsamples_fish(k,i)*(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j) )*log(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j)+p_const);
        }
    }
  obs_catchdm_biom.initialize();
  obs_catchdf_biom.initialize();
  avgp.initialize();
  for(i=1;i<=nobs_fish_discm;i++)
  {
    obs_catchdm_biom(yrs_fish_discm(i)) = catch_tot(yrs_fish_discm(i)) * (obs_p_fish_discm(1,i)+obs_p_fish_discm(2,i))*wtm;
    avgp += obs_p_fish_discm(1,i)+obs_p_fish_discm(2,i);
  }
  avgp=avgp/nobs_fish_discm;
  avgpf = (obs_p_fish_discf(1)+obs_p_fish_discf(2)+obs_p_fish_discf(3)+obs_p_fish_discf(4)+obs_p_fish_discf(5))/5.;
  for(i=styr;i<=(yrs_fish_discm(1)-1);i++)
  {
    obs_catchdm_biom(i)=avgp*catch_tot(i)*wtm;
  } 
  obs_catchdm_biom(endyr)=avgp*catch_tot(endyr)*wtm;
  // cout <<obs_catchdm_biom<<endl;exit(1);
  for(i=styr;i<=endyr;i++)
  {
    obs_catchtot_biom(i)=obs_catchdm_biom(i)+catch_ret(i);
    obs_catchdf_biom(i)=(avgpf*catch_disc(1,i))*wtf(2);
  }
}

void model_parameters::get_obs_survey(void)
{
  ofstream& post= *pad_post;
   //don't do the move with fraction_new_error -adjustment made in prelimn calcs for obs_p_srv1_len1
  //move new shells to old shells for mature only due to error in shell condition age
         for(k=1;k<=2;k++){
          obs_p_srv1_len(2,1,k)=obs_p_srv1_len1(2,1,k);
          obs_p_srv1_len(2,2,k)=obs_p_srv1_len1(2,2,k);
         }
          obs_p_srv1_len(1)=obs_p_srv1_len1(1);
 int mat;
  offset(4)=0.0;
  offset_srv2.initialize();
  for(k=1; k<=2;k++) //sex
   for (i=1; i <= nobs_srv1_length; i++)
    for (j=1; j<=nlenm; j++)
    {
          offset(4)-=nsamples_srv1_length(2,2,k,i)*(obs_p_srv1_len(1,1,k,i,j)+
                       obs_p_srv1_len(1,2,k,i,j)+
                       obs_p_srv1_len(2,1,k,i,j)+obs_p_srv1_len(2,2,k,i,j))*
                            log(obs_p_srv1_len(1,1,k,i,j)+obs_p_srv1_len(1,2,k,i,j)+
                                   obs_p_srv1_len(2,1,k,i,j)+obs_p_srv1_len(2,2,k,i,j)+p_const);
    }
    sumobs=sum(obs_p_srv2_lend(1,1,2)(5,nlenm)+obs_p_srv2_lend(1,2,2)(5,nlenm)+obs_p_srv2_lend(1,1,1)(5,nlenm)+obs_p_srv2_lend(1,2,1)(5,nlenm));
    sumobs2=sum(obs_p_srv2_lend(2,1,2)(5,nlenm)+obs_p_srv2_lend(2,2,2)(5,nlenm)+obs_p_srv2_lend(2,1,1)(5,nlenm)+obs_p_srv2_lend(2,2,1)(5,nlenm));
  for(il2=1;il2<=2;il2++) 
    for(mat=1;mat<=2;mat++)
    for(j=1;j<=nlenm;j++)
    {
      if(j>4){
        obs_p_srv2_len(1,il2,mat,j)=obs_p_srv2_lend(1,il2,mat,j)/sumobs;
        obs_p_srv2_len(2,il2,mat,j)=obs_p_srv2_lend(2,il2,mat,j)/sumobs2;
        }
        else{                    
        obs_p_srv2_len(1,il2,mat,j)=0.0;
        obs_p_srv2_len(2,il2,mat,j)=0.0;
        }
   }
  //offset
 for(k=1; k<=2; k++)
  for(il2=1; il2<=2 ;il2++)
   for(mat=1;mat<=2;mat++)
    for(j=1;j<=nlenm;j++)
    {
      offset_srv2(k) -= nsamples_srv2_length(k,il2,mat)* obs_p_srv2_len(k,il2,mat,j)*log(obs_p_srv2_len(k,il2,mat,j)+p_const);
    }
    sumobs10=sum(obs_p_srv10_lend(1,1,2)(1,nlenm)+obs_p_srv10_lend(1,2,2)(1,nlenm)+obs_p_srv10_lend(1,1,1)(1,nlenm)+obs_p_srv10_lend(1,2,1)(1,nlenm));
    sumobs102=sum(obs_p_srv10_lend(2,1,2)(1,nlenm)+obs_p_srv10_lend(2,2,2)(1,nlenm)+obs_p_srv10_lend(2,1,1)(1,nlenm)+obs_p_srv10_lend(2,2,1)(1,nlenm));
 for(il2=1;il2<=2;il2++) 
   for(mat=1;mat<=2;mat++)
    for(j=1;j<=nlenm;j++)
    {
       obs_p_srv10_len(1,il2,mat,j)=obs_p_srv10_lend(1,il2,mat,j)/sumobs10;
       obs_p_srv10_len(2,il2,mat,j)=obs_p_srv10_lend(2,il2,mat,j)/sumobs102;
    }
  //offset
   offset_srv10.initialize();
 for(k=1; k<=2; k++)
  for(il2=1; il2<=2 ;il2++)
   for(mat=1;mat<=2;mat++)
	for(j=1;j<=nlenm;j++)
	 offset_srv10(k) -= nsamples_srv10_length(k,il2,mat)* obs_p_srv10_len(k,il2,mat,j)*log(obs_p_srv10_len(k,il2,mat,j)+p_const);
  obs_srv1_num.initialize();
  obs_srv1_biom.initialize();
  obs_srv1_bioms.initialize();
  obs_srv1_spbiom.initialize();
  obs_srv1_spnum.initialize();
  obs_srv2_spbiom.initialize();
  obs_srv10_spbiom.initialize();
  for(i=1;i<=nobs_srv1;i++)
    obs_srv1t(yrs_srv1(i))=obs_srv1(i);
  for(mat=1;mat<=2;mat++)  //maturity status
   for(l=1;l<=2;l++)  //shell condition
    for(k=1;k<=2;k++)  //sex
    for (i=1; i <= nobs_srv1_length; i++)
    {
         obs_srv1_num(k,yrs_srv1_length(i)) += obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i));
        if(k<2){
          obs_srv1_bioms(k,yrs_srv1_length(i))+=obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtf(mat);
         obs_srv1_biom(yrs_srv1_length(i))+= obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtf(mat);
         }
        else{
         obs_srv1_bioms(k,yrs_srv1_length(i))+=obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtm;
         obs_srv1_biom(yrs_srv1_length(i))+= obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtm;
        }        
       if(mat>1)
       {
        if(k<2){
        obs_srv1_spbiom(k,yrs_srv1_length(i)) += obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtf(mat);
        }
        else{
        obs_srv1_spbiom(k,yrs_srv1_length(i)) += obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtm;
        }
        obs_srv1_spnum(l,k,yrs_srv1_length(i)) += sum(obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i)));
       }
    }
   obs_srv2_spbiom(1,1) = (obs_p_srv2_len(1,1,2)*sumobs*wtf(2))*1000.;
   obs_srv2_spbiom(2,1) = (obs_p_srv2_len(2,1,2)*sumobs2*wtf(2))*1000.;
   obs_srv2_spbiom(1,2) = (obs_p_srv2_len(1,2,2)*sumobs*wtm)*1000.;
   obs_srv2_spbiom(2,2) = (obs_p_srv2_len(2,2,2)*sumobs2*wtm)*1000.;
   obs_srv10_spbiom(1,1) = (obs_p_srv10_len(1,1,2)*sumobs10*wtf(2))*1000.;
   obs_srv10_spbiom(2,1) = (obs_p_srv10_len(2,1,2)*sumobs102*wtf(2))*1000.;
   obs_srv10_spbiom(1,2) = (obs_p_srv10_len(1,2,2)*sumobs10*wtm)*1000.;
   obs_srv10_spbiom(2,2) = (obs_p_srv10_len(2,2,2)*sumobs102*wtm)*1000.;
}

dvector model_parameters::bootstrap(dvector& data,int& sample_size, long int seed)
{
  ofstream& post= *pad_post;
 //data contains the data from which the bootstrap sample is to be taken
 //sample_size is the size of the desired bootstrap sample
 //seed is a seed for the random number generator
 //use like this
 //  dvector boot_sample=bootstrap(data,100,1231);  takes a sample with replacement of size
 // 100 from data using seed 1231 and puts it into the vector boot_sample.
 dvector tmp(1,sample_size);  //declare vector tmp of length sample_size
 tmp.fill_randu(seed);       // gives uniform random numbers from 0 to 1 size sample_size
 tmp=(tmp*data.size())+1.0;      // makes tmp random numbers from 0 to length of data+1
 ivector iselect(tmp);       // takes the integer part of tmp and puts it into a vector iselect
 tmp=data(iselect);          // takes the bootstrap sample from data puts into tmp
 return(tmp);                // return the bootstrap sample in tmp
}

void model_parameters::get_growth(void)
{
  ofstream& post= *pad_post;
  int ilen,il2,sex;
  double devia ;
  dvariable lensum;
  len_len.initialize();
  rec_len.initialize();
  for(ilen=1;ilen<=nlenm;ilen++)
  { 
  if(growth_switch==1)
   {
     mean_length(1,ilen) = ((af+bf*length_bins(ilen))*(1-cumd_norm((length_bins(ilen)-deltaf)/st_gr))+((af+(bf-bf1)*deltaf)+bf1*length_bins(ilen))*(cumd_norm((length_bins(ilen)-deltaf)/st_gr)));
     mean_length(2,ilen) = am+bm*length_bins(ilen);
 //     mean_length(2,ilen) = ((am+bm*length_bins(ilen))*(1-cumd_norm((length_bins(ilen)-deltam)/st_gr))+((am+(bm-b1)*deltam)+b1*length_bins(ilen))*(cumd_norm((length_bins(ilen)-deltam)/st_gr)));
  }
  }
   for (sex=1;sex<=2;sex++)
  {
    for(ilen=1;ilen<=nlenm;ilen++)
    {
      // subract the 2.5 from the midpoint of the length bin to get the lower bound
      alpha(sex) = (mean_length(sex,ilen)-(length_bins(ilen)-2.5))/growth_beta(sex);
      lensum = 0;
      for(il2=ilen;il2<=ilen+min(8,nlenm-ilen);il2++)
      {
          devia = length_bins(il2)+2.5-length_bins(ilen);
           len_len(sex,ilen,il2) = pow(devia,(alpha(sex)-1.))*mfexp(-devia/growth_beta(sex));
          lensum += len_len(sex,ilen,il2);
          // cout << devia<<" "<<alpha<<" "<<growth_beta(sex)<<endl;
   //standardize so each row sums to 1.0
      }  
      len_len(sex,ilen) /= lensum;
    }
  }
 // Fraction recruiting
  recsum=0.0;
  alpha_rec=alpha1_rec/beta_rec;
 //do only first 6 bins
  for(ilen=1;ilen<=6;ilen++)
  {
    devia = length_bins(ilen)+2.5-length_bins(1);
    rec_len(ilen) = pow(devia,alpha_rec-1.)*mfexp(-devia/beta_rec);
    recsum += rec_len(ilen);
  }
  //standardize so each row sums to 1.0
  for(ilen=1;ilen<=nlenm;ilen++)
   rec_len(ilen) = rec_len(ilen)/recsum;
}

void model_parameters::get_moltingp(void)
{
  ofstream& post= *pad_post;
  for(j=1;j<=nlenm;j++)
  {
 //  these next two lines are logistic molting females then males
     moltp(1,j)=1-(1./(1.+mfexp(-1.*moltp_af*(length_bins(j)-moltp_bf))));
      moltp(2,j)=1-(1./(1.+mfexp(-1.*moltp_am*(length_bins(j)-moltp_bm))));
     //set molting prob for mature females at 0.0
      moltp_mat(1,j)=0.0;
     if(phase_moltingp>0)
       {
       moltp_mat(2,j)=1-(1./(1.+mfexp(-1.*moltp_ammat*(length_bins(j)-moltp_bmmat))));
        }
        else
        {
 //mature male molting probability equal zero - terminal molt--I ithought we all agreed that there is a terminal molt?
          moltp_mat(2,j)=0.0;
        }
    }
}

void model_parameters::get_selectivity(void)
{
  ofstream& post= *pad_post;
  dvariable tempSeln, tempSelo;
  int x;
 for(iy=styr;iy<=endyr+Nproj;iy++)
   {
    if(iy<1978)
    {
     fish_sel50_mn(iy)=mfexp(log_avg_sel50_mn);
     if(active(log_avg_sel50_mo))
      fish_sel50_mo(iy)=mfexp(log_avg_sel50_mo);
    }
    if(iy>=1978 & iy<endyr)
    {      
     fish_sel50_mn(iy)=mfexp(log_avg_sel50_mn+log_sel50_dev_mn(iy));
     if(active(log_avg_sel50_mo))
       fish_sel50_mo(iy)=mfexp(log_avg_sel50_mo+log_sel50_dev_mo(iy));
    }
    if(iy>=endyr)  
    {  
     tempSeln = 0;
	 tempSelo = 0;
     for(x=(endyr-sel_avg_Nyrs);x<=(endyr-1);x++)
	 {
      tempSeln += log_sel50_dev_mn(x);
	  tempSelo += log_sel50_dev_mo(x);
	 }
	 tempSeln = tempSeln/sel_avg_Nyrs;
	 tempSelo = tempSelo/sel_avg_Nyrs;
     fish_sel50_mn(iy)=mfexp(log_avg_sel50_mn+tempSeln);
     if(active(log_avg_sel50_mo))
       fish_sel50_mo(iy)=mfexp(log_avg_sel50_mo+tempSelo);
    }	
  //logistic selectivity curve
 for (j=1;j<=nlenm;j++)
  { 
    sel(1,iy,j)=1./(1.+mfexp(-1.*fish_slope_mn*(length_bins(j)-fish_sel50_mn(iy))));
   //set new and old sel same
    sel(2,iy,j)= sel(1,iy,j);
	//for dome shaped add this part
    if(phase_fishsel>0)
     {
      tmp2=1./(1.+mfexp(fish_slope_mn2*(length_bins(j)-fish_sel50_mn2)));
      tmp3=1./(1.+mfexp(fish_slope_mo2*(length_bins(j)-fish_sel50_mo2)));
      sel(1,iy,j)=sel(1,iy,j)*tmp2;
      sel(2,iy,j)=sel(1,iy,j);
      }
      sel_ret(1,iy,j)=1./(1.+mfexp(-1.*fish_fit_slope_mn*(length_bins(j)-fish_fit_sel50_mn)));
      sel_fit(1,iy,j)=sel_ret(1,iy,j)*sel(1,iy,j);
      sel_ret(2,iy,j)=sel_ret(1,iy,j);
      sel_fit(2,iy,j)=sel_ret(2,iy,j)*sel(2,iy,j);
   }
  }  //end of year loop
 //female discards ascending logistic curve 
  for (j=1;j<=nlenm;j++)
  {
   for(iy=yrs_fish_discf(1);iy<(yrs_fish_discf(nobs_fish_discf)+.5);iy++)
    sel_discf(iy,j)=1./(1.+mfexp(-1.*fish_disc_slope_f*(length_bins(j)-mfexp(fish_disc_sel50_f+log_dev_50f(iy)))));
   for(iy=yrs_fish_discf(nobs_fish_discf);iy<endyr+Nproj;iy++)
     sel_discf(iy,j)=1./(1.+mfexp(-1.*fish_disc_slope_f*(length_bins(j)-mfexp(fish_disc_sel50_f+log_dev_50f(endyr-1)))));
    sel_trawl(1,j)=1./(1.+mfexp(-1.*fish_disc_slope_tf*(length_bins(j)-fish_disc_sel50_tf)));
    sel_trawl(2,j)=sel_trawl(1,j); 
 if(j<=nsellen_srv1)
   {
    if(somertonsel<0)
    { 
       sel_srv3(1,j)=1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv3_sel50_f)/(srv3_sel95_f-srv3_sel50_f)));
       sel_srv3(2,j)=sel_som(1)/(1.+mfexp(-1.*(sel_som(2)+(sel_som(3)*length_bins(j)))));           
    }
    else
    { 
    sel_srv3(2,j)=srv3_q*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv3_sel50)/(srv3_sel95-srv3_sel50)));
    sel_srv3(1,j)=srv3_q_f*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv3_sel50_f)/(srv3_sel95_f-srv3_sel50_f)));
    }
   sel_srvind(2,j)= srvind_q/(1+mfexp(-selsmo09ind(j)));
   sel_srvind(1,j)= srvind_q_f* 1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srvind_sel50_f)/(srvind_sel95_f-srvind_sel50_f)));
   sel_srvnmfs(1,j)= sel_srvind(1,j)*sel_srv3(1,j);
   sel_srvnmfs(2,j)= sel_srvind(2,j)*sel_srv3(2,j);
   // surv sel 2010 study area
   sel_srv10ind(2,j)= srv10ind_q*1/(1+mfexp(-selsmo10ind(j)));
   sel_srv10ind(1,j)= srv10ind_q_f* 1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv10ind_sel50_f)/(srv10ind_sel95_f-srv10ind_sel50_f)));
   sel_srv10nmfs(1,j)=sel_srv10ind(1,j)*sel_srv3(1,j);
   sel_srv10nmfs(2,j)=sel_srv10ind(2,j)*sel_srv3(2,j); 
   // this sets time periods 1 and 2 survey selectivities to somerton otto as well
  if(survsel1_phase<0){
    sel_srv1(2,j)=sel_srv3(2,j);
    sel_srv2(2,j)=sel_srv3(2,j);
    sel_srv1(1,j)=sel_srv3(1,j);
    sel_srv2(1,j)=sel_srv3(1,j);
    }
    else
    { 
    sel_srv1(2,j)=srv1_q*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv1_sel50)/(srv1_sel95-srv1_sel50)));
    sel_srv2(2,j)=srv2_q*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv2_sel50)/(srv2_sel95-srv2_sel50)));
    sel_srv1(1,j)=srv1_q_f*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv1_sel50)/(srv1_sel95-srv1_sel50)));
    sel_srv2(1,j)=srv2_q_f*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv2_sel50)/(srv2_sel95-srv2_sel50)));
    }
    }
    else
   {
    sel_srv1(1,j)=sel_srv1(1,j-1);
    sel_srv1(2,j)=sel_srv1(2,j-1);
    }
  }
  for(iy=styr;iy<=endyr+Nproj;iy++)
    {
     maxsel_fish=max(sel(1,iy));
     if(maxsel_fish<max(sel(2,iy)))
	   maxsel_fish=max(sel(2,iy));
	 sel(1,iy)=sel(1,iy)/maxsel_fish;
     sel(2,iy)=sel(2,iy)/maxsel_fish;
     sel_fit(1,iy)=elem_prod(sel_ret(1,iy),sel(1,iy));
     sel_fit(2,iy)=elem_prod(sel_ret(2,iy),sel(2,iy));
    }
}

void model_parameters::get_mortality(void)
{
  ofstream& post= *pad_post;
 int x;
 fmort(styr,endyr-1) = mfexp(log_avg_fmort+fmort_dev);
 fmort(endyr)=mfexp(log_avg_fmort);  // last year not fit to
 fmortdf(styr,endyr-1)=mfexp(log_avg_fmortdf+fmortdf_dev);
 fmortdf(endyr)=mfexp(log_avg_fmortdf);
 for(x=styr;x<=1991;x++)
  fmortt(x)=mfexp(log_avg_fmortt+fmortt_dev_era1(x));
 for(x=1992;x<=endyr-1;x++)
  fmortt(x)=mfexp(log_avg_fmortt+fmortt_dev_era2(x));
  fmortt(endyr)=mfexp(log_avg_fmortt);
  for (i=styr;i<=endyr;i++)
 {
   if(i>=yrs_fish_discf(1) && i<endyr)
   {
    Fdiscf(i)=sel_discf(i)*fmortdf(i);
   }
   else
   {
    sel_discf_e=(1./(1.+mfexp(-1.*fish_disc_slope_f*(length_bins-mfexp(fish_disc_sel50_f)))));
    Fdiscf(i)=sel_discf_e*fmortdf(i);
    }
   Fdisct(1,i)=sel_trawl(1)*fmortt(i);   
   Fdisct(2,i)=sel_trawl(2)*fmortt(i);   
   for(k=1;k<=2;k++) //over new (k=1) and old (k=2) shell...
   { 
     F(k,i) = sel(k,i)*fmort(i);       
     F_ret(k,i)=sel_fit(k,i)*fmort(i);
     Fmat(k,i) = sel(k,i)*fmort(i);       
     Fmat_ret(k,i)=sel_fit(k,i)*fmort(i);
     Fimat(k,i) = sel(k,i)*fmort(i);      // Fishing mort on immature males new or old shell  
     Fimat_ret(k,i)= sel_fit(k,i)*fmort(i);
     Ftot(1,k,i)=Fdiscf(i) + Fdisct(1,i);
     S(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Smat(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Simat(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Smat(2,k,i)=mfexp(-1.0*(Fmat(k,i)+Fdisct(2,i)));
     Simat(2,k,i)=mfexp(-1.0*(Fimat(k,i)+Fdisct(2,i)));
   } 
 }
}

void model_parameters::get_numbers_at_len(void)
{
  ofstream& post= *pad_post;
  int itmp;
  natlength_new.initialize();
  natlength_old.initialize();
  natlength_inew.initialize();
  natlength_iold.initialize();
  natlength_mnew.initialize();
  natlength_mold.initialize();
  natlength_i.initialize();
  natlength_mat.initialize();
  natlength.initialize();
  efspbio_matetime.initialize();
  emspbio_matetime.initialize();
  mspbio_old_matetime.initialize();
  fspbio_new_matetime.initialize();
  efspbio_new_matetime.initialize();
  fspnum_new_matetime.initialize();
  efspnum_matetime.initialize();
  emspnum_old_matetime.initialize();
  mspnum_matetime.initialize();
  bio_males_gt101.initialize();
  num_males_gt101.initialize();
   popn.initialize();
   explbiom.initialize();
   pred_bio.initialize();
   pred_srv1.initialize();
   pred_srv1_biom.initialize();
   pred_srv1_bioms.initialize();
   fspbio.initialize();
   mspbio.initialize(); 
    //this is dumb but necesssary?
  rec_dev(1)					=rec_devf;
  rec_dev(2)					=rec_devf;
  if(active(rec_devm))
	rec_dev(2) = rec_devm;
  mean_log_rec(1)	 			=mean_log_rec_f;
  mean_log_rec(2)	 			=mean_log_rec_f;
  if(active(mean_log_rec_m))
	mean_log_rec(2) = mean_log_rec_m;
  dvar_matrix mat_not_average(1,2,1,nlenm);
   int jk;
  for(k=1;k<=2;k++)
   for(jk=1;jk<=nlenm;jk++)
    mat_not_average(k,jk)	=1.0-maturity_average(k,jk);
  for(j=1;j<=12;j++)
   {
      natlength_mnew(1,styr,j)	 	= maturity_est(1,j)*mfexp(fnatlen_styr(1,j));
      natlength_inew(1,styr,j) 		= (1.0-maturity_est(1,j))*mfexp(fnatlen_styr(1,j));
      natlength_mold(1,styr,j) 		= mfexp(fnatlen_styr(2,j));
    }
   for(j=13;j<=nlenm;j++)
	{
      natlength_mnew(1,styr,j)	= 0.0;
      natlength_inew(1,styr,j) 	= 0.0;
      natlength_mold(1,styr,j) 	= 0.0;
    }
   natlength_iold(1,styr) 		= 0.0;
  //males
   natlength_inew(2,styr) 		= mfexp(mnatlen_styr(1));
   natlength_mnew(2,styr)		= mfexp(mnatlen_styr(2));
   natlength_mold(2,styr) 		= 0.0;
   natlength_iold(2,styr) 		= 0.0;
  for(k=1;k<=2;k++)  //k is sex 
   {
	 natlength_new(k,styr)		=natlength_inew(k,styr)+natlength_mnew(k,styr);
      natlength_old(k,styr)		=natlength_iold(k,styr)+natlength_mold(k,styr);
      natlength_mat(k,styr)		=natlength_mnew(k,styr)+natlength_mold(k,styr);
      natlength_i(k,styr)		=natlength_inew(k,styr)+natlength_iold(k,styr);
      natlength(k,styr)			=natlength_inew(k,styr)+natlength_iold(k,styr) + natlength_mnew(k,styr)+natlength_mold(k,styr);
    }
 for (ipass=styr;ipass<=endyr;ipass++)
    get_num_at_len_yr();
}

void model_parameters::get_num_at_len_yr(void)
{
  ofstream& post= *pad_post;
 i = ipass;
 // if(i < endyr)
 // {
 for(k=1;k<=2;k++)
  {
      // Numbers advancing to new shell...
      dvar_vector tmp 			= elem_prod(moltp(k)*mfexp(-(1-catch_midpt(i))*M(k)),elem_prod(Simat(k,1,i),mfexp(-catch_midpt(i)*M(k))*natlength_inew(k,i)));
      //is this is the same as:
	  // dvar_vector tmp 				= elem_prod(moltp(k),elem_prod(Simat(k,1,i),mfexp(*M(k))*natlength_inew(k,i)));
      natlength_new(k,i+1)		=  tmp * len_len(k);
      dvar_vector tmpo 			= elem_prod(moltp(k)*mfexp(-(1-catch_midpt(i))*M(k)), elem_prod(Simat(k,2,i),mfexp(-catch_midpt(i)*M(k))*natlength_iold(k,i)));
      natlength_new(k,i+1)		+=  tmpo * len_len(k);
      natlength_iold(k,i+1) 	= mfexp(-(1-catch_midpt(i))*M(k))*(elem_prod(Simat(k,1,i),mfexp(-catch_midpt(i)*M(k))*natlength_inew(k,i)) +  elem_prod(Simat(k,2,i),mfexp(-catch_midpt(i)*M(k))*natlength_iold(k,i))) - tmp-tmpo;
      dvar_vector tmpm 			= elem_prod(moltp_mat(k)*mfexp(-(1-catch_midpt(i))*M_matn(k)), elem_prod(Smat(k,1,i),mfexp(-catch_midpt(i)*M_matn(k))*natlength_mnew(k,i)));
      natlength_mnew(k,i+1)		=   tmpm * len_len(k);
      dvar_vector tmpmo			= elem_prod(moltp_mat(k)*mfexp(-(1-catch_midpt(i))*M_mato(k)), elem_prod(Smat(k,2,i),mfexp(-catch_midpt(i)*M_mato(k))*natlength_mold(k,i)));
      natlength_mnew(k,i+1) 	+=   tmpmo * len_len(k);
      natlength_mold(k,i+1)		= (mfexp(-(1-catch_midpt(i))*M_matn(k)) * elem_prod(Smat(k,1,i),mfexp(-catch_midpt(i)*M_matn(k))*natlength_mnew(k,i))) +  (mfexp(-(1-catch_midpt(i))*M_mato(k)) * elem_prod(Smat(k,2,i),mfexp(-catch_midpt(i)*M_mato(k))*natlength_mold(k,i))) - tmpm-tmpmo;
 //this is for estimating the fractino of new shell that move to old shell to fit
      natlength_mnew(k,i+1) 	+= elem_prod(maturity_est(k),natlength_new(k,i+1));
      natlength_inew(k,i+1) 	= elem_prod(1.0-maturity_est(k),natlength_new(k,i+1));
 if(i<=endyr)
      natlength_inew(k,i+1) 	+= mfexp(mean_log_rec(k)+rec_dev(k,i))*rec_len*proprecn;
 else
      natlength_inew(k,i+1) 	+= FutRec*rec_len*proprecn;
      natlength_mat(k,i+1)    	 = natlength_mnew(k,i+1) + natlength_mold(k,i+1);
      natlength_i(k,i+1)      	 = natlength_inew(k,i+1) + natlength_iold(k,i+1);
      natlength(k,i+1)        	 = natlength_mat(k,i+1)  + natlength_i(k,i+1);
      natlength_old(k,i+1)    	 = natlength_mold(k,i+1) + natlength_iold(k,i+1);
      natlength_new(k,i+1)    	 = natlength_inew(k,i+1) + natlength_mnew(k,i+1);
      natl_inew_fishtime(k,i)	 = mfexp(-catch_midpt(i)*M(k))*natlength_inew(k,i);
      natl_iold_fishtime(k,i) 	= mfexp(-catch_midpt(i)*M(k))*natlength_iold(k,i);
      natl_mnew_fishtime(k,i)	=  mfexp(-catch_midpt(i)*M_matn(k))*natlength_mnew(k,i);
      natl_mold_fishtime(k,i) 	=  mfexp(-catch_midpt(i)*M_mato(k))*natlength_mold(k,i);
      natl_new_fishtime(k,i)	= natl_inew_fishtime(k,i)+natl_mnew_fishtime(k,i); 
      natl_old_fishtime(k,i) 	= natl_iold_fishtime(k,i)+natl_mold_fishtime(k,i);
      mspbio_matetime(i) 		= (elem_prod(Smat(2,1,i)*mfexp(-(spmo)*M_matn(2)),mfexp(-catch_midpt(i)*M_matn(2))*natlength_mnew(2,i))+elem_prod(Smat(2,2,i)*mfexp(-(spmo)*M_mato(2)),mfexp(-catch_midpt(i)*M_mato(2))*natlength_mold(2,i)))*wtm;
      fspbio_matetime(i) 		= (elem_prod(Smat(1,1,i)*mfexp(-(spmo)*M_matn(1)),mfexp(-catch_midpt(i)*M_matn(1))*natlength_mnew(1,i))+elem_prod(Smat(1,2,i)*mfexp(-(spmo)*M_mato(1)),mfexp(-catch_midpt(i)*M_mato(1))*natlength_mold(1,i)))*wtf(2);
      mspbio_fishtime(i) 		= (natl_mnew_fishtime(2,i)+natl_mold_fishtime(2,i))*wtm;
      fspbio_fishtime(i) 		= (natl_mnew_fishtime(1,i)+natl_mold_fishtime(1,i))*wtf(2);
    if(k>1)
	{
     for(j=1;j<=nlenm;j++)
       {
        emspnum_old_matetime(i) 		+= Smat(2,2,i,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(i)*M_mato(2))*natlength_mold(2,i,j);
        mspnum_matetime(i) 				+= Smat(2,1,i,j)*mfexp(-(spmo)*M_matn(2))*mfexp(-catch_midpt(i)*M_matn(2))*natlength_mnew(2,i,j) + Smat(2,2,i,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(i)*M_mato(2))*natlength_mold(2,i,j);
        mspbio_old_matetime(i) 			+= (Smat(2,2,i,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(i)*M_mato(2))*natlength_mold(2,i,j))*wtm(j);
        fspnum_new_matetime(i) 			+= (Smat(1,1,i,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(i)*M_matn(1))*natlength_mnew(1,i,j));
        fspbio_new_matetime(i)			+= (Smat(1,1,i,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(i)*M_matn(1))*natlength_mnew(1,i,j))*wtf(2,j);
        efspnum_matetime(i) 			+= (Smat(1,1,i,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(i)*M_matn(1))*natlength_mnew(1,i,j)+Smat(1,2,i,j)*mfexp(-(spmo)*M_mato(1))*mfexp(-catch_midpt(i)*M_mato(1))*natlength_mold(1,i,j));
        if(j>15){
          num_males_gt101(i)			+= natl_inew_fishtime(k,i,j) + natl_iold_fishtime(k,i,j) + natl_mnew_fishtime(k,i,j) + natl_mold_fishtime(k,i,j);
          bio_males_gt101(i)			+= (natl_inew_fishtime(k,i,j) + natl_iold_fishtime(k,i,j) + natl_mnew_fishtime(k,i,j) + natl_mold_fishtime(k,i,j))*wtm(j);
        if(j<17){
              num_males_gt101(i)		=num_males_gt101(i)*0.5;
              bio_males_gt101(i)		=bio_males_gt101(i)*0.5;
            }
          }
       }
      efspbio_matetime(i)				=fspbio_matetime(i);
      emspbio_matetime(i)				=mspbio_old_matetime(i);
      if(emspnum_old_matetime(i)<(efspnum_matetime(i)/mate_ratio))
        efspbio_matetime(i)				=fspbio_matetime(i)*((emspnum_old_matetime(i)*mate_ratio)/efspnum_matetime(i));
      efspbio_new_matetime(i)			=fspbio_new_matetime(i);
      if(emspnum_old_matetime(i)<fspnum_new_matetime(i)/mate_ratio)
        efspbio_new_matetime(i)			=fspbio_new_matetime(i)*((emspnum_old_matetime(i)*mate_ratio)/fspnum_new_matetime(i));
       popn_lmale_new(i) 				= natl_inew_fishtime(2,i)*sel(1,i)+natl_mnew_fishtime(2,i)*sel(1,i);
       popn_lmale_old(i) 				=  natl_iold_fishtime(2,i)*sel(2,i)+natl_mold_fishtime(2,i)*sel(2,i);
       popn_lmale_bio(i)   			 	= elem_prod(natl_new_fishtime(2,i),sel(1,i))*wtm+elem_prod(natl_old_fishtime(2,i),sel(2,i))*wtm;
       popn_lmale(i)       				= popn_lmale_new(i)+popn_lmale_old(i);
       popn_fit_new(i)  				= natl_inew_fishtime(2,i)*sel_fit(1,i)+natl_mnew_fishtime(2,i)*sel_fit(1,i);
       popn_fit_old(i)  				= natl_iold_fishtime(2,i)*sel_fit(2,i)+natl_mold_fishtime(2,i)*sel_fit(2,i);
       popn_fit(i)          			= popn_fit_new(i) + popn_fit_old(i);
      } //end of if(k>1) loop
     if(k<2)
         {
         if(i>=yrs_fish_discf(1) && i<endyr)
		  {
		    popn_disc(1,i)       = (natl_new_fishtime(1,i)+natl_old_fishtime(1,i))*sel_discf(i);
          }
         else
		   {
             popn_disc(1,i)       = (natl_new_fishtime(1,i)+natl_old_fishtime(1,i))*sel_discf_e;
            }
         }
    }
  //isn't this all repeated code?
 if(i ==endyr)
  {
  natl_inew_fishtime(1,endyr) = (mfexp(-catch_midpt(endyr)*M(1))*natlength_inew(1,endyr));
  natl_iold_fishtime(1,endyr) = (mfexp(-catch_midpt(endyr)*M(1))*natlength_iold(1,endyr));
  natl_inew_fishtime(2,endyr) = (mfexp(-catch_midpt(endyr)*M(2))*natlength_inew(2,endyr));
  natl_iold_fishtime(2,endyr) = (mfexp(-catch_midpt(endyr)*M(2))*natlength_iold(2,endyr));
  natl_mnew_fishtime(1,endyr) = (mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr));
  natl_mold_fishtime(1,endyr) = (mfexp(-catch_midpt(endyr)*M_mato(1))*natlength_mold(1,endyr));
  natl_mnew_fishtime(2,endyr) = (mfexp(-catch_midpt(endyr)*M_matn(2))*natlength_mnew(2,endyr));
  natl_mold_fishtime(2,endyr) = (mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr));
  natl_new_fishtime(1,endyr)  = natl_inew_fishtime(1,endyr)+natl_mnew_fishtime(1,endyr);
  natl_old_fishtime(1,endyr)  = natl_iold_fishtime(1,endyr)+natl_mold_fishtime(1,endyr);
  natl_new_fishtime(2,endyr)  = natl_inew_fishtime(2,endyr)+natl_mnew_fishtime(2,endyr);
  natl_old_fishtime(2,endyr)  = natl_iold_fishtime(2,endyr)+natl_mold_fishtime(2,endyr);
  mspbio_matetime(endyr) = (elem_prod(Smat(2,1,endyr)*mfexp(-(spmo)*M_matn(2)),mfexp(-catch_midpt(endyr)*M_matn(2))*natlength_mnew(2,endyr))+elem_prod(Smat(2,2,endyr)*mfexp(-(spmo)*M_mato(2)),mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr)))*wtm;
  fspbio_matetime(endyr) = (elem_prod(Smat(1,1,endyr)*mfexp(-(spmo)*M_matn(1)),mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr))+elem_prod(Smat(1,2,endyr)*mfexp(-(spmo)*M_mato(1)),mfexp(-catch_midpt(endyr)*M_mato(1))*natlength_mold(1,endyr)))*wtf(2);
  mspbio_fishtime(endyr) = (natl_mnew_fishtime(2,endyr)+natl_mold_fishtime(2,endyr))*wtm;
  fspbio_fishtime(endyr) = (natl_mnew_fishtime(1,endyr)+natl_mold_fishtime(1,endyr))*wtf(2);
  for(j=1;j<=nlenm;j++)
   {
      emspnum_old_matetime(endyr) 	+= Smat(2,2,endyr,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr,j);
      mspnum_matetime(endyr) 		+= Smat(2,1,endyr,j)*mfexp(-(spmo)*M_matn(2))*mfexp(-catch_midpt(endyr)*M_matn(2))*natlength_mnew(2,endyr,j) + Smat(2,2,endyr,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr,j);
        mspbio_old_matetime(endyr) 	+= (Smat(2,2,endyr,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr,j))*wtm(j);
        fspnum_new_matetime(endyr) 	+= (Smat(1,1,endyr,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr,j));
        fspbio_new_matetime(endyr)	+= (Smat(1,1,endyr,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr,j))*wtf(2,j);
        efspnum_matetime(endyr) 	+= (Smat(1,1,endyr,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr,j)+Smat(1,2,endyr,j)*mfexp(-(spmo)*M_mato(1))*mfexp(-catch_midpt(endyr)*M_mato(1))*natlength_mold(1,endyr,j));
        if(j>15){
          num_males_gt101(endyr)	+= natl_inew_fishtime(2,endyr,j) + natl_iold_fishtime(2,endyr,j) + natl_mnew_fishtime(2,endyr,j) + natl_mold_fishtime(2,endyr,j);
          bio_males_gt101(endyr)	+= (natl_inew_fishtime(2,endyr,j) + natl_iold_fishtime(2,endyr,j) + natl_mnew_fishtime(2,endyr,j) + natl_mold_fishtime(2,endyr,j))*wtm(j);
        if(j<17){
              num_males_gt101(endyr)=num_males_gt101(endyr)*0.5;
              bio_males_gt101(endyr)=bio_males_gt101(endyr)*0.5;
            }
         }
    }
    efspbio_matetime(endyr)			=fspbio_matetime(endyr);
    emspbio_matetime(endyr)			=mspbio_old_matetime(endyr);
   if(emspnum_old_matetime(endyr)	>=efspnum_matetime(endyr)/mate_ratio)
    efspbio_matetime(endyr)			=fspbio_matetime(endyr)*((emspnum_old_matetime(endyr)*mate_ratio)/efspnum_matetime(endyr));
   efspbio_new_matetime(endyr)		=fspbio_new_matetime(endyr);
   if(emspnum_old_matetime(endyr)<fspnum_new_matetime(endyr)/mate_ratio)
    efspbio_new_matetime(endyr)		=fspbio_new_matetime(endyr)*((emspnum_old_matetime(endyr)*mate_ratio)/fspnum_new_matetime(endyr));
  popn_lmale_new(endyr)  			= natl_inew_fishtime(2,endyr)*sel(1,endyr)+natl_mnew_fishtime(2,endyr)*sel(1,endyr);
  popn_lmale_old(endyr) 			= natl_iold_fishtime(2,endyr)*sel(2,endyr)+natl_mold_fishtime(2,endyr)*sel(2,endyr);
  popn_lmale_bio(endyr)   			= elem_prod(natl_new_fishtime(2,endyr),sel(1,endyr))*wtm+elem_prod(natl_old_fishtime(2,endyr),sel(2,endyr))*wtm;
  popn_lmale(endyr)    				= popn_lmale_new(endyr) + popn_lmale_old(endyr);
  popn_fit_new(endyr)  				= natl_inew_fishtime(2,endyr)*sel_fit(1,endyr)+natl_mnew_fishtime(2,endyr)*sel_fit(1,endyr);
  popn_fit_old(endyr)  				= natl_iold_fishtime(2,endyr)*sel_fit(2,endyr)+natl_mold_fishtime(2,endyr)*sel_fit(2,endyr);
  popn_fit(endyr)    				= popn_fit_new(endyr)+ popn_fit_old(endyr);
  popn_disc(1,endyr) 				= (natl_new_fishtime(1,endyr)+natl_old_fishtime(1,endyr))*sel_discf(yrs_fish_discf(nobs_fish_discf));
  }
  //predicted survey values 
  for(k=1;k<=2;k++)
      {
       if(i<1982){
        totn_srv1(k,i)=(natlength(k,i)*sel_srv1(k));
        }
       if(i>1981 && i<1988){
        totn_srv1(k,i)=(natlength(k,i)*sel_srv2(k));
        }
       if(i>1987){
        totn_srv1(k,i)=(natlength(k,i)*sel_srv3(k));
        if(i==2009){
	        		totn_srv2(1,k)=(natlength(k,i)*sel_srvind(k));
       				totn_srv2(2,k)=(natlength(k,i)*sel_srvnmfs(k));
           		   }
           if(i==2010){
	        		totn_srv10(1,k)=(natlength(k,i)*sel_srv10ind(k));
       				totn_srv10(2,k)=(natlength(k,i)*sel_srv10nmfs(k));
           		   }
            }
        totn_trawl(k,i)= (natl_new_fishtime(k,i)+natl_old_fishtime(k,i))*sel_trawl(k);
   }
    predpop_sexr(i)=sum(natlength(1,i))/(sum(natlength(1,i))+sum(natlength(2,i)));
    fspbio(i)=natlength_mat(1,i)*wtf(2);
    mspbio(i)=natlength_mat(2,i)*wtm;
    if(i<1982)
	{
     fspbio_srv1(i) = q1*natlength_mat(1,i)*elem_prod(wtf(2),sel_srv1(1));
     mspbio_srv1(i) = q1*natlength_mat(2,i)*elem_prod(wtm,sel_srv1(2));
     fspbio_srv1_num(1,i) = q1*natlength_mnew(1,i)*sel_srv1(1);
     mspbio_srv1_num(1,i) = q1*natlength_mnew(2,i)*sel_srv1(2);
     fspbio_srv1_num(2,i) = q1*natlength_mold(1,i)*sel_srv1(1);
     mspbio_srv1_num(2,i) = q1*natlength_mold(2,i)*sel_srv1(2);
    }
    if(i>1981 && i<1988)
	{
     fspbio_srv1(i) = q1*natlength_mat(1,i)*elem_prod(wtf(2),sel_srv2(1));
     mspbio_srv1(i) = q1*natlength_mat(2,i)*elem_prod(wtm,sel_srv2(2));
     fspbio_srv1_num(1,i) = q1*natlength_mnew(1,i)*sel_srv2(1);
     mspbio_srv1_num(1,i) = q1*natlength_mnew(2,i)*sel_srv2(2);
     fspbio_srv1_num(2,i) = q1*natlength_mold(1,i)*sel_srv2(1);
     mspbio_srv1_num(2,i) = q1*natlength_mold(2,i)*sel_srv2(2);
    }
    if(i>1987)
	{
     fspbio_srv1(i) = q1*natlength_mat(1,i)*elem_prod(wtf(2),sel_srv3(1));
     mspbio_srv1(i) = q1*natlength_mat(2,i)*elem_prod(wtm,sel_srv3(2));
     fspbio_srv1_num(1,i) = q1*natlength_mnew(1,i)*sel_srv3(1);
     mspbio_srv1_num(1,i) = q1*natlength_mnew(2,i)*sel_srv3(2);
     fspbio_srv1_num(2,i) = q1*natlength_mold(1,i)*sel_srv3(1);
     mspbio_srv1_num(2,i) = q1*natlength_mold(2,i)*sel_srv3(2);
      if(i==2009)
        {
	     fspbio_srv2_ind = natlength_mat(1,i)*elem_prod(wtf(2),sel_srvind(1));
   		 mspbio_srv2_ind = natlength_mat(2,i)*elem_prod(wtm,sel_srvind(2));
		 fspbio_srv2_nmfs = natlength_mat(1,i)*elem_prod(wtf(2),sel_srvnmfs(1));
   		 mspbio_srv2_nmfs = natlength_mat(2,i)*elem_prod(wtm,sel_srvnmfs(2));
         pred_p_srv2_len_ind(1)=elem_prod(sel_srvind(1),natlength(1,i))/(totn_srv2(1,1)+totn_srv2(1,2));
         pred_p_srv2_len_nmfs(1)=elem_prod(sel_srvnmfs(1),natlength(1,i))/(totn_srv2(2,1)+totn_srv2(2,2));
         pred_p_srv2_len_ind(2)=elem_prod(sel_srvind(2),natlength(2,i))/(totn_srv2(1,1)+totn_srv2(1,2));
         pred_p_srv2_len_nmfs(2)=elem_prod(sel_srvnmfs(2),natlength(2,i))/(totn_srv2(2,1)+totn_srv2(2,2));
        }
      if(i==2010)
       {
	     fspbio_srv10_ind = natlength_mat(1,i)*elem_prod(wtf(2),sel_srv10ind(1));
   		 mspbio_srv10_ind = natlength_mat(2,i)*elem_prod(wtm,sel_srv10ind(2));
		 fspbio_srv10_nmfs = natlength_mat(1,i)*elem_prod(wtf(2),sel_srv10nmfs(1));
   		 mspbio_srv10_nmfs = natlength_mat(2,i)*elem_prod(wtm,sel_srv10nmfs(2));
         pred_p_srv10_len_ind(1)=elem_prod(sel_srv10ind(1),natlength(1,i))/(totn_srv10(1,1)+totn_srv10(1,2));
         pred_p_srv10_len_nmfs(1)=elem_prod(sel_srv10nmfs(1),natlength(1,i))/(totn_srv10(2,1)+totn_srv10(2,2));
         pred_p_srv10_len_ind(2)=elem_prod(sel_srv10ind(2),natlength(2,i))/(totn_srv10(1,1)+totn_srv10(1,2));
         pred_p_srv10_len_nmfs(2)=elem_prod(sel_srv10nmfs(2),natlength(2,i))/(totn_srv10(2,1)+totn_srv10(2,2));
	   }
    }
 for(k=1;k<=2;k++)
  {
    if(i<1982)
	{
      pred_srv1(k,i) = q1*elem_prod(natlength(k,i),sel_srv1(k));
     if(k<2)
	 {
	   pred_srv1_bioms(k,i) = q1*((natlength_inew(k,i)*elem_prod(sel_srv1(k),wtf(1)))+((natlength_mnew(k,i)+natlength_mold(k,i))*elem_prod(sel_srv1(k),wtf(2))));
      }
     else
	 {
      pred_srv1_bioms(k,i) = q1*(natlength(k,i)*elem_prod(sel_srv1(k),wtm));
      }
      pred_p_srv1_len_new(1,k,i)=elem_prod(sel_srv1(k),natlength_inew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(1,k,i)=elem_prod(sel_srv1(k),natlength_iold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_new(2,k,i)=elem_prod(sel_srv1(k),natlength_mnew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(2,k,i)=elem_prod(sel_srv1(k),natlength_mold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
    }
  if(i>1981 && i<1988)
  {
      pred_srv1(k,i) = q1*elem_prod(natlength(k,i),sel_srv2(k));
     if(k<2)
	 {
      pred_srv1_bioms(k,i) = q1*((natlength_inew(k,i)*elem_prod(sel_srv2(k),wtf(1)))+((natlength_mnew(k,i)+natlength_mold(k,i))*elem_prod(sel_srv2(k),wtf(2))));
      }
     else
	 {
      pred_srv1_bioms(k,i) = q1*(natlength(k,i)*elem_prod(sel_srv2(k),wtm));
      }
      pred_p_srv1_len_new(1,k,i)=elem_prod(sel_srv2(k),natlength_inew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(1,k,i)=elem_prod(sel_srv2(k),natlength_iold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_new(2,k,i)=elem_prod(sel_srv2(k),natlength_mnew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(2,k,i)=elem_prod(sel_srv2(k),natlength_mold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
    }
   if(i>1987)
   {
      pred_srv1(k,i) = q1*elem_prod(natlength(k,i),sel_srv3(k));
     if(k<2)
	 {
      pred_srv1_bioms(k,i) = q1*((natlength_inew(k,i)*elem_prod(sel_srv3(k),wtf(1)))+((natlength_mnew(k,i)+natlength_mold(k,i))*elem_prod(sel_srv3(k),wtf(2))));
      }
     else
	 {
      pred_srv1_bioms(k,i) = q1*(natlength(k,i)*elem_prod(sel_srv3(k),wtm));
      }
      pred_p_srv1_len_new(1,k,i)=elem_prod(sel_srv3(k),natlength_inew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(1,k,i)=elem_prod(sel_srv3(k),natlength_iold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_new(2,k,i)=elem_prod(sel_srv3(k),natlength_mnew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(2,k,i)=elem_prod(sel_srv3(k),natlength_mold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
    }
    //next line used to fix q1 to 1.0 - problem is if you start from a bin file, even if the bounds
    // are set different in the tpl file the program will take to value from the bin file and use that 
    //   pred_srv1(i)=1.0*(natage(i)*elem_prod(sel_srv1,wt));
    explbiom(i)=((natlength_old(2,i)*elem_prod(sel(2,i),wtm))+(natlength_new(2,i)*elem_prod(sel(1,i),wtm)));
    popn(i)+=sum(natlength(k,i));
   }  //end of k loop
    pred_bio(i)+=natlength_inew(1,i)*wtf(1)+(natlength_mnew(1,i)+natlength_mold(1,i))*wtf(2)+(natlength_inew(2,i)+natlength_mnew(2,i)+natlength_mold(2,i))*wtm;
  depletion = pred_bio(endyr) / pred_bio(styr);
  fspbios=fspbio;
  mspbios=mspbio_matetime;
    legal_males(i)=0.;
    legal_srv_males(i)=0.;
    legal_srv_males_n(i)=0.0;
    legal_srv_males_o(i)=0.0;
        legal_males(i)=0.5*natlength(2,i,16);
        legal_srv_males_n(i)=0.5*natlength_new(2,i,16);
        legal_srv_males_o(i)=0.5*natlength_old(2,i,16);
        legal_males_bio(i)=legal_males(i)*wtm(16);
        if(i<1982)
		{
        legal_srv_males(i)=0.5*natlength(2,i,16)*sel_srv1(2,16);
        legal_srv_males_n(i)=0.5*natlength_new(2,i,16)*sel_srv1(2,16);
        legal_srv_males_o(i)=0.5*natlength_old(2,i,16)*sel_srv1(2,16);
        }
        if(i>1981 && i<1988)
		{
         legal_srv_males(i)=0.5*natlength(2,i,16)*sel_srv2(2,16);
         legal_srv_males_n(i)=0.5*natlength_new(2,i,16)*sel_srv2(2,16);
         legal_srv_males_o(i)=0.5*natlength_old(2,i,16)*sel_srv2(2,16);
        }
        if(i>1987)
		{
         legal_srv_males(i)=0.5*natlength(2,i,16)*sel_srv3(2,16);
         legal_srv_males_n(i)=0.5*natlength_new(2,i,16)*sel_srv3(2,16);
         legal_srv_males_o(i)=0.5*natlength_old(2,i,16)*sel_srv3(2,16);
        }
        legal_srv_males_bio(i)=legal_srv_males(i)*wtm(16);
	for(j=17;j<=nlenm;j++)
    {
        legal_males(i)+=natlength(2,i,j);
        legal_males_bio(i)+=natlength(2,i,j)*wtm(j);
        if(i<1982)
		{
          legal_srv_males(i)+=natlength(2,i,j)*sel_srv1(2,j);
          legal_srv_males_n(i)+=natlength_new(2,i,j)*sel_srv1(2,j);
          legal_srv_males_o(i)+=natlength_old(2,i,j)*sel_srv1(2,j);
          legal_srv_males_bio(i)+=natlength(2,i,j)*sel_srv1(2,j)*wtm(j);
         }
        if(i>1981 && i<1988)
		{
          legal_srv_males(i)+=natlength(2,i,j)*sel_srv2(2,j);
          legal_srv_males_n(i)+=natlength_new(2,i,j)*sel_srv2(2,j);
          legal_srv_males_o(i)+=natlength_old(2,i,j)*sel_srv2(2,j);
          legal_srv_males_bio(i)+=natlength(2,i,j)*sel_srv2(2,j)*wtm(j);
          }
        if(i>1987)
		{
          legal_srv_males(i)+=natlength(2,i,j)*sel_srv3(2,j);
          legal_srv_males_n(i)+=natlength_new(2,i,j)*sel_srv3(2,j);
          legal_srv_males_o(i)+=natlength_old(2,i,j)*sel_srv3(2,j);
          legal_srv_males_bio(i)+=natlength(2,i,j)*sel_srv3(2,j)*wtm(j);
          }
    }
    // legal_malesd=legal_males;
    // recf_sd=mfexp(mean_log_rec(1)+rec_dev(1));
    // recm_sd=mfexp(mean_log_rec(2)+rec_dev(2));
}

void model_parameters::get_catch_at_len(void)
{
  ofstream& post= *pad_post;
 //cout<<" begin catch at len"<<endl;
    pred_catch.initialize();
    pred_catch_ret.initialize();
    pred_catch_gt101.initialize();
    pred_catch_no_gt101.initialize();
	// take out catch all at once - survey is start, catch_midpt has midpoint of the
    // catch(fishing season weighted by the catch)
    // then rest of nat mort and growth.
  for (i=styr;i<=endyr;i++)
  { 
   for(k=1;k<=2;k++)
    {      
    for (j = 1 ; j<= nlenm; j++)
    {
        if(k>1)
        {
        if(j>15)
		{
             pred_catch_no_gt101(i) += (Fimat(1,i,j)/(Fimat(1,i,j)+Fdisct(2,i,j)))*natl_inew_fishtime(2,i,j)*(1-Simat(2,1,i,j)) + 
                                         (Fmat(1,i,j)/(Fmat(1,i,j)+Fdisct(2,i,j)))*natl_mnew_fishtime(2,i,j)*(1-Smat(2,1,i,j))+ 
                                          (Fimat(2,i,j)/(Fimat(2,i,j)+Fdisct(2,i,j)))*natl_iold_fishtime(2,i,j)*(1-Simat(2,2,i,j))+
                                           (Fmat(2,i,j)/(Fmat(2,i,j)+Fdisct(2,i,j)))*natl_mold_fishtime(2,i,j)*(1-Smat(2,2,i,j));
             pred_catch_gt101(i)+= (Fimat(1,i,j)/(Fimat(1,i,j)+Fdisct(2,i,j)))*natl_inew_fishtime(2,i,j)*(1-Simat(2,1,i,j)) + 
                                         (Fmat(1,i,j)/(Fmat(1,i,j)+Fdisct(2,i,j)))*natl_mnew_fishtime(2,i,j)*(1-Smat(2,1,i,j))+ 
                                          (Fimat(2,i,j)/(Fimat(2,i,j)+Fdisct(2,i,j)))*natl_iold_fishtime(2,i,j)*(1-Simat(2,2,i,j))+
                                           (Fmat(2,i,j)/(Fmat(2,i,j)+Fdisct(2,i,j)))*natl_mold_fishtime(2,i,j)*(1-Smat(2,2,i,j)) * wtm(j);
           if(j<17)
		   {
               pred_catch_gt101(i)=pred_catch_gt101(i)*0.5;
               pred_catch_no_gt101(i)=pred_catch_no_gt101(i)*0.5;
            }
        }
        catch_lmale_new(i,j) = (Fimat(1,i,j)/(Fimat(1,i,j)+Fdisct(2,i,j)))*natl_inew_fishtime(2,i,j)*(1-Simat(2,1,i,j))+(Fmat(1,i,j)/(Fmat(1,i,j)+Fdisct(2,i,j)))*natl_mnew_fishtime(2,i,j)*(1-Smat(2,1,i,j)); 
        catch_lmale_old(i,j) = (Fimat(2,i,j)/(Fimat(2,i,j)+Fdisct(2,i,j)))*natl_iold_fishtime(2,i,j)*(1-Simat(2,2,i,j))+(Fmat(2,i,j)/(Fmat(2,i,j)+Fdisct(2,i,j)))*natl_mold_fishtime(2,i,j)*(1-Smat(2,2,i,j));
        catch_lmale(i,j)= catch_lmale_new(i,j)+catch_lmale_old(i,j);
        pred_catch(i) += catch_lmale(i,j)*wtm(j);
        catch_male_ret_new(i,j) = (Fimat_ret(1,i,j)/(Fimat(1,i,j)+Fdisct(2,i,j)))*natl_inew_fishtime(2,i,j)*(1-Simat(2,1,i,j)) + (Fmat_ret(1,i,j)/(Fmat(1,i,j)+Fdisct(2,i,j)))*natl_mnew_fishtime(2,i,j)*(1-Smat(2,1,i,j)); 
        catch_male_ret_old(i,j) = (Fimat_ret(2,i,j)/(Fimat(2,i,j)+Fdisct(2,i,j)))*natl_iold_fishtime(2,i,j)*(1-Simat(2,2,i,j))  + (Fmat_ret(2,i,j)/(Fmat(2,i,j)+Fdisct(2,i,j)))*natl_mold_fishtime(2,i,j)*(1-Smat(2,2,i,j));
        catch_male_ret(i,j)=catch_male_ret_new(i,j)+catch_male_ret_old(i,j);
        pred_catch_ret(i) += catch_male_ret(i,j)*wtm(j);
       } //end k
    catch_discp(1,i,j) = (Fdiscf(i,j)/(Fdiscf(i,j)+Fdisct(1,i,j))) *(natl_new_fishtime(1,i,j)+natl_old_fishtime(1,i,j))*(1-mfexp(-1.0*(Fdiscf(i,j)+Fdisct(1,i,j))));
    // pred_catch_trawl
    } 
	if(k<2)
	{
      pred_catch_disc(k,i)=catch_discp(k,i)*wtf(2);
    }
    else
	{
      pred_catch_disc(k,i)=catch_discp(k,i)*wtm;
    }
    //  cout<<" 3 catlen"<<endl;
   }
   pred_catch_trawl(i)= (natl_new_fishtime(1,i)+natl_old_fishtime(1,i))*elem_prod(1-mfexp(-1.0*Fdisct(1,i)),wtf(2))+ (natl_new_fishtime(2,i)+natl_old_fishtime(2,i))*elem_prod(1-mfexp(-1.0*Fdisct(2,i)),wtm);
   pred_p_fish_fit(1,i)=(elem_prod(sel_fit(1,i),natl_new_fishtime(2,i)))/(popn_fit(i));
   pred_p_fish_fit(2,i)=(elem_prod(sel_fit(2,i),natl_old_fishtime(2,i)))/(popn_fit(i));
   pred_p_fish(1,i)=(elem_prod(sel(1,i),natl_new_fishtime(2,i))/(popn_lmale(i)));
   pred_p_fish(2,i)=(elem_prod(sel(2,i),natl_old_fishtime(2,i))/(popn_lmale(i)));
   if(i>=yrs_fish_discf(1) && i<endyr)
   {
   pred_p_fish_discf(i)=elem_prod(sel_discf(i),(natl_new_fishtime(1,i)+natl_old_fishtime(1,i)))/popn_disc(1,i);
   }
   else
   {
   pred_p_fish_discf(i)=elem_prod(sel_discf_e,(natl_new_fishtime(1,i)+natl_old_fishtime(1,i)))/popn_disc(1,i);
   }
   pred_p_trawl(1,i)=elem_prod(sel_trawl(1),(natl_new_fishtime(1,i)+natl_old_fishtime(1,i)))/(totn_trawl(1,i)+totn_trawl(2,i));
   pred_p_trawl(2,i)=elem_prod(sel_trawl(2),(natl_new_fishtime(2,i)+natl_old_fishtime(2,i)))/(totn_trawl(1,i)+totn_trawl(2,i));
   preds_sexr(i)=totn_srv1(1,i)/(totn_srv1(1,i)+totn_srv1(2,i));
  }
}

void model_parameters::evaluate_the_objective_function(void)
{
  ofstream& post= *pad_post;
  len_likeyr.initialize();
  len_like.initialize();
  len_like_srv.initialize();
  len_like10_ind.initialize();
  len_like10_nmfs.initialize();
  sel_like=0.;
  fpen=.0;
  like_initsmo=0.0;
  ghl_like=0.0;
  rec_like=.0;
  surv_like=.0;
  surv2_like=0.0;
  surv3_like=0.0;
  catch_like1=.0;
  catch_like2=.0;
  catch_likef=.0;
  catch_liket=.0;
  like_q=0.0;
  like_natm.initialize();
  like_af.initialize();
  like_bf.initialize();
  like_bf1.initialize();
  like_mmat.initialize();
  like_am.initialize();
  like_bm.initialize();
  sexr_like.initialize();
  like_natm.initialize();
  sumrecf.initialize();
  sumrecm.initialize();
  sel_like_50m.initialize();
   sel_avg_like.initialize();
  like_srvsel.initialize();
  like_initnum.initialize();
  Fout.initialize();
  len_like_ex=0.0;
  selsmo_like=0.0;
  //  catch_likef=0.;
  growth_like=0.0;
  f.initialize();
  largemale_like=0.0;
 if (active(rec_devf))
   {
    rec_like = like_wght_recf*norm2(rec_devf)+like_wght_rec*norm2(rec_devm);
    Fout(1)=rec_like;
    f += rec_like;
   // if(active(rec_devm))
	// {
      // for(i=styr;i<endyr-retro_years;i++)
        // sexr_like += like_wght_sexr*square((mean_log_rec(1)+rec_devf(i))-(mean_log_rec(2)+rec_devm(i)));
    // }
    }  //end of active rec_devf
 //constraint on intial numbers in small length bins for old shell males
   for(i=1;i<=5;i++)
    {
     like_initnum += old_shell_constraint  * square(mfexp(mnatlen_styr(2,i)));
     }
     f+= like_initnum;
    Fout(2)=like_initnum;
 //jim said take the log--CSS if only adding the likecomp when active...why is it commented out?
  if(active(log_sel50_dev_mn))
   {
    sel_like_50m=like_wght_sel50*norm2(first_difference((log_avg_sel50_mn+log_sel50_dev_mn)(styr,endyr-1)));
    sel_like_50m+=like_wght_sel50*norm2(first_difference((log_avg_sel50_mo+log_sel50_dev_mo)(styr,endyr-1)));
 //   f +=sel_like_50m;
  }
  int ii;
  int ij;
  int ik;
  int mat;
    //==========Retained fishery lengths likelihood component===============
  for (i=1; i <= (nobs_fish - retro_years); i++)
  {
    ii=yrs_fish(i);
   for (j=1; j<=nlenm; j++)
    len_like(1)-=nsamples_fish(1,i)*(obs_p_fish_ret(1,i,j)+obs_p_fish_ret(2,i,j))*log(pred_p_fish_fit(2,ii,j)+pred_p_fish_fit(1,ii,j)+p_const); //old and new together
   for(k=1;k<=2;k++)
    effn_fish_ret(k,ii)=1/(norm2(pred_p_fish_fit(k,ii)-obs_p_fish_ret(k,i))/(pred_p_fish_fit(1,ii)*(1-pred_p_fish_fit(1,ii))+pred_p_fish_fit(2,ii)*(1-pred_p_fish_fit(2,ii))));
  }
 for (i=1; i <= (nobs_fish_discm - retro_years); i++)
  {
    ij=yrs_fish_discm(i);
   for (j=1; j<=nlenm; j++)
    len_like(2)-=disclen_mult(1)*nsamples_fish_discm(1,i)*(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j))*log(pred_p_fish(1,ij,j)+pred_p_fish(2,ij,j)+p_const);
   for(k=1;k<=2;k++)
    effn_fish_tot(k,ij)=1/(norm2(pred_p_fish(k,ij)-obs_p_fish_tot(k,i))/(pred_p_fish(1,ij)*(1-pred_p_fish(1,ij))+pred_p_fish(2,ij)*(1-pred_p_fish(2,ij))));
  }
  //==========Discards female lengths likelihood component=============== 
  for (i=1; i <= (nobs_fish_discf - retro_years); i++)
  {
    ik=yrs_fish_discf(i);
   for (j=1; j<=nlenm; j++)
	len_like(3)-=nsamples_fish_discf(i)*(obs_p_fish_discf(i,j))*log(pred_p_fish_discf(ik,j)+p_const);
  }
  //==========Trawl fishery lengths likelihood component===============
  for (i=1; i <= (nobs_trawl - retro_years); i++)
  {
    ij=yrs_trawl(i);
  //trawlfishery length likelihood sur
    for(k=1;k<=2;k++)
      for (j=1; j<=nlenm; j++)
       len_like(5)-=nsamples_trawl(k,i)*(obs_p_trawl(k,i,j))*log(pred_p_trawl(k,ij,j)+p_const);
  }
  //add the offset to the likelihood   
  len_like(1)-=offset(1);
  len_like(2)-=offset(2);
  len_like(3)-=offset(3);
  len_like(5)-=offset(5);
 for (j=1; j<=nlenm; j++)
   {
   	len_like(6)-=1.*nsamples_srv2_length(1,1,1)*((obs_p_srv2_len(1,1,1,j)+obs_p_srv2_len(1,1,2,j))) *log(pred_p_srv2_len_ind(1,j)+p_const);   			  	
    len_like(6)-=1.*nsamples_srv2_length(1,2,1)*((obs_p_srv2_len(1,2,1,j)+obs_p_srv2_len(1,2,2,j))) *log(pred_p_srv2_len_ind(2,j)+p_const);
    len_like(7)-=1.*nsamples_srv2_length(2,1,1)*(obs_p_srv2_len(2,1,1,j)+obs_p_srv2_len(2,1,2,j)) *log(pred_p_srv2_len_nmfs(1,j)+p_const);
    len_like(7)-=1.*nsamples_srv2_length(2,2,1)*(obs_p_srv2_len(2,2,1,j)+obs_p_srv2_len(2,2,2,j)) *log(pred_p_srv2_len_nmfs(2,j)+p_const);
   	len_like10_ind-=1.*nsamples_srv10_length(1,1,1)*((obs_p_srv10_len(1,1,1,j)+obs_p_srv10_len(1,1,2,j))) *log(pred_p_srv10_len_ind(1,j)+p_const);   			  	
    len_like10_ind-=1.*nsamples_srv10_length(1,2,1)*((obs_p_srv10_len(1,2,1,j)+obs_p_srv10_len(1,2,2,j))) *log(pred_p_srv10_len_ind(2,j)+p_const);
    len_like10_nmfs-=1.*nsamples_srv10_length(2,1,1)*(obs_p_srv10_len(2,1,1,j)+obs_p_srv10_len(2,1,2,j)) *log(pred_p_srv10_len_nmfs(1,j)+p_const);
    len_like10_nmfs-=1.*nsamples_srv10_length(2,2,1)*(obs_p_srv10_len(2,2,1,j)+obs_p_srv10_len(2,2,2,j)) *log(pred_p_srv10_len_nmfs(2,j)+p_const);
    }               
   len_like(6)-=offset_srv2(1);
   len_like(7)-=offset_srv2(2);
   len_like10_ind-=offset_srv10(1);
   len_like10_nmfs-=offset_srv10(2);
 //===========survey lengths likelihood components==================
   for(k=1;k<=2;k++)  //sex
   {
    for (i=1; i <=(nobs_srv1_length - retro_years); i++)
    {
      ii=yrs_srv1_length(i);
      //survey likelihood 
      for (j=1; j<=nlenm; j++)
      {
 //  obs(maturity, SC, sex, year), pred(maturity,sex, year)
 // this is for mature new and old shell together
       len_like(4)-=nsamples_srv1_length(2,1,k,i)*(obs_p_srv1_len(2,1,k,i,j)+obs_p_srv1_len(2,2,k,i,j))*log(pred_p_srv1_len_new(2,k,ii,j)+pred_p_srv1_len_old(2,k,ii,j)+p_const);
 //immature new and old together
        len_like(4)-=nsamples_srv1_length(1,1,k,i)*(obs_p_srv1_len(1,1,k,i,j)+obs_p_srv1_len(1,2,k,i,j))*log(pred_p_srv1_len_new(1,k,ii,j)+pred_p_srv1_len_old(1,k,ii,j)+p_const);
 //immature new and old together in likelihood indices are (mat,shell,sex,year,length)
        len_likeyr(1,1,k,i)-=nsamples_srv1_length(1,1,k,i)*(obs_p_srv1_len(1,1,k,i,j)+obs_p_srv1_len(1,2,k,i,j))*log(pred_p_srv1_len_new(1,k,ii,j)+pred_p_srv1_len_old(1,k,ii,j)+p_const);
        len_like_srv(1,2,k)=0.0;
        len_like_srv(1,1,k)-=nsamples_srv1_length(1,1,k,i)*(obs_p_srv1_len(1,1,k,i,j)+obs_p_srv1_len(1,2,k,i,j))*log(pred_p_srv1_len_new(1,k,ii,j)+pred_p_srv1_len_old(1,k,ii,j)+p_const);
     //mature
        len_likeyr(2,1,k,i)-=nsamples_srv1_length(2,1,k,i)*(obs_p_srv1_len(2,1,k,i,j))*log(pred_p_srv1_len_new(2,k,ii,j)+p_const);
        len_likeyr(2,2,k,i)-=nsamples_srv1_length(2,2,k,i)*(obs_p_srv1_len(2,2,k,i,j))*log(pred_p_srv1_len_old(2,k,ii,j)+p_const);
        len_like_srv(2,1,k)-=nsamples_srv1_length(2,1,k,i)*(obs_p_srv1_len(2,1,k,i,j))*log(pred_p_srv1_len_new(2,k,ii,j)+p_const);
        len_like_srv(2,2,k)-=nsamples_srv1_length(2,2,k,i)*(obs_p_srv1_len(2,2,k,i,j))*log(pred_p_srv1_len_old(2,k,ii,j)+p_const);
    }  //j loop     
    effn_srv1(1,1,k,ii)=1./(norm2(pred_p_srv1_len_new(1,k,ii)-obs_p_srv1_len(1,1,k,i))/(pred_p_srv1_len_new(1,1,ii)*(1-pred_p_srv1_len_new(1,1,ii))+pred_p_srv1_len_new(1,2,ii)*(1-pred_p_srv1_len_new(1,2,ii))+pred_p_srv1_len_old(1,1,ii)*(1-pred_p_srv1_len_old(1,1,ii))+pred_p_srv1_len_old(1,2,ii)*(1-pred_p_srv1_len_old(1,2,ii))));
  if(k>1) //what is this if statement about?? if and else have the same statement inside...
	{
    effn_srv1(1,2,k,ii)=0.0;
    }
    else
	{
    //effn is zero for old shell immature females
    effn_srv1(1,2,k,ii)=0.0;
    }
    effn_srv1(2,1,k,ii)=1./(norm2(pred_p_srv1_len_new(2,k,ii)-(obs_p_srv1_len(2,1,k,i)))/(pred_p_srv1_len_new(2,1,ii)*(1-pred_p_srv1_len_new(2,1,ii))+pred_p_srv1_len_new(2,2,ii)*(1-pred_p_srv1_len_new(2,2,ii))+pred_p_srv1_len_old(2,1,ii)*(1-pred_p_srv1_len_old(2,1,ii))+pred_p_srv1_len_old(2,2,ii)*(1-pred_p_srv1_len_old(2,2,ii))));
    effn_srv1(2,2,k,ii)=1./(norm2(pred_p_srv1_len_old(2,k,ii)-(obs_p_srv1_len(2,2,k,i)))/(pred_p_srv1_len_new(2,1,ii)*(1-pred_p_srv1_len_new(2,1,ii))+pred_p_srv1_len_new(2,2,ii)*(1-pred_p_srv1_len_new(2,2,ii))+pred_p_srv1_len_old(2,1,ii)*(1-pred_p_srv1_len_old(2,1,ii))+pred_p_srv1_len_old(2,2,ii)*(1-pred_p_srv1_len_old(2,2,ii))));
    } // year loop
  }  //k (sex) loop
   like_initsmo = init_yr_len_smooth_f *norm2(first_difference(mnatlen_styr(1)))+init_yr_len_smooth_m *norm2(first_difference(first_difference(mnatlen_styr(2))))+norm2(first_difference(fnatlen_styr(1)))+norm2(first_difference(fnatlen_styr(2)));
   Fout(30)=like_initsmo;
   f += like_initsmo;
  //===================extra weight for start year length comp.==========================   
  //==immature males of both shel conditions
   len_like(4)-=nsamples_srv1_length(1,1,2,1)*(obs_p_srv1_len(1,1,2,1))*log(pred_p_srv1_len_new(1,2,styr)+p_const);
   if(moltp(2,nlenm)<0.99)
    len_like(4)-=nsamples_srv1_length(1,2,2,1)*(obs_p_srv1_len(1,2,2,1))*log(pred_p_srv1_len_old(1,2,styr)+p_const);
 //  obs(maturity, SC, sex, year), pred(maturity,sex, year)
    multi=1.0;
	//==mature both sexes and shell conditions
    len_like_ex-=multi*nsamples_srv1_length(2,1,2,1)*(obs_p_srv1_len(2,1,2,1))*log(pred_p_srv1_len_new(2,2,styr)+p_const);
    len_like_ex-= multi*nsamples_srv1_length(2,2,2,1)*(obs_p_srv1_len(2,2,2,1))*log(pred_p_srv1_len_old(2,2,styr)+p_const);
    len_like_ex-=multi*nsamples_srv1_length(2,1,1,1)*(obs_p_srv1_len(2,1,1,1))*log(pred_p_srv1_len_new(2,1,styr)+p_const);
    len_like_ex-=multi*nsamples_srv1_length(2,2,1,1)*(obs_p_srv1_len(2,2,1,1))*log(pred_p_srv1_len_old(2,1,styr)+p_const);
    len_like_ex-=multi*nsamples_srv1_length(1,1,1,1)*(obs_p_srv1_len(1,1,1,1)+obs_p_srv1_len(1,2,1,1))*log(pred_p_srv1_len_new(1,1,styr)+pred_p_srv1_len_old(1,1,styr)+p_const);
   Fout(25)=len_like_ex;
   f+=len_like_ex;
  len_like(4)-=offset(4);
 for(mat=1;mat<=2;mat++)
  for(k=1;k<=1;k++)
    len_like_srv(mat,k,1)-=offset_srv(mat,k,1);
  if(active(selsmo10ind))
  {
    selsmo_like= selsmo_wght *norm2(first_difference(first_difference(selsmo10ind)));
    f+=selsmo_like;
    Fout(28)+=selsmo_like;
   }
  if(active(selsmo09ind))
  {
    selsmo_like= selsmo_wght *norm2(first_difference(first_difference(selsmo09ind)));
    f+=selsmo_like;
    Fout(28)+=selsmo_like;
   }
   if(active(log_dev_50f))
   {
    discf_like= femSel_wght *norm2(first_difference(first_difference(log_dev_50f)));
    f+=discf_like;
    Fout(29)=discf_like;
    }
  f+=like_wght(1)*len_like(1);     // Retained fishery
  Fout(3)=like_wght(1)*len_like(1);
  f+=like_wght(2)*len_like(2);     // Total (ret+disc)
  Fout(4)=like_wght(2)*len_like(2);
  f+=like_wght(3)*len_like(3);     // Female
  Fout(5)=like_wght(3)*len_like(3);
  f+=like_wght(4)*len_like(4);     // Survey
  Fout(6)=like_wght(4)*len_like(4);
  f+=like_wght(7)*len_like(5);     // Trawl
  Fout(7)=like_wght(7)*len_like(5);
  f+=len_like(6);                 // Industry Survey BSFRF length
  Fout(8)=len_like(6);
  f+=len_like(7);                // Industry survey NMFS length
  Fout(9)=len_like(7);
  f+=len_like10_ind;                 // Industry 2010 Survey BSFRF length
  Fout(26)=len_like10_ind;
  f+=len_like10_nmfs;                // Industry 2010 survey NMFS length
  Fout(27)=len_like10_nmfs;
 like_natm=0.0;
 if(active(Mmult))
 {  
  like_natm   = natm_mult_wght  * square((Mmult    - 1.0)    / natm_mult_var);
  f += like_natm;
  Fout(10)=like_natm;
 }
 if(active(Mmult_imat))
 {  
  like_natm   += natm_Immult_wght  * square((Mmult_imat    - 1.0)    / natm_Immult_var);
  f += like_natm;
  Fout(10)+=like_natm;
 }
  if(active(Mmultf))
 {  
  like_natm   += natm_mult_wght  * square((Mmultf    - 1.0)    / natm_mult_var);
  f += like_natm;
  Fout(10)+=like_natm;
 }
   like_mat=0.0;
  if(active(mateste))
   {
     like_mat = smooth_mat_wght *norm2(first_difference(first_difference(mateste)));
     // like_mat += mat_est_wght *norm2(maturity_est(2)-maturity_logistic)/(mat_est_vs_obs_sd *mat_est_vs_obs_sd );
     f += like_mat;
     Fout(11)+=like_mat;
    }
 if(active(matestfe))
    {
     like_mat += smooth_mat_wght_f *norm2(first_difference(first_difference(matestfe)));
     // like_mat += 1.0*norm2(maturity_est(1)-maturity_average(1))/(mat_est_vs_obs_sd_f *mat_est_vs_obs_sd_f );
     f += like_mat;
     Fout(11)+=like_mat;    
    }
   if(active(am))
   {
	 for(i=1;i<=nobs_growm ;i++)
	  like_am+= growth_data_wght_m *pow(malegrowdaty(i)-(am+bm*malegrowdatx(i)),2);       
      f += like_am;
      Fout(12)=like_am;
    }
   if(active(af))
   {
	for(i=1;i<=nobs_growf ;i++)
     like_af+= growth_data_wght_f *pow(femalegrowdaty(i)-((af+bf*femalegrowdatx(i))*(1-cumd_norm((femalegrowdatx(i)-deltaf)/st_gr))+((af+(bf-bf1)*deltaf)+bf1*femalegrowdatx(i))*(cumd_norm((femalegrowdatx(i)-deltaf)/st_gr))),2);
	 f += like_af;
     Fout(13)=like_af;
   }
  // ==================Fit to indices===============================
  // (lognormal) weight each years estimate by 1/(2*variance) - use cv of biomass in sqrt(log(cv^2+1)) as sd of log(biomass) 
 for(i=1;i<=nobs_srv1;i++)
  {
    biom_tmp(1,yrs_srv1(i))=fspbio_srv1(yrs_srv1(i));
    biom_tmp(2,yrs_srv1(i))=mspbio_srv1(yrs_srv1(i));
  }
  for(i=1;i<=nobs_srv1;i++)
    {
     cv_srv1(1,yrs_srv1(i))=like_wght(5)*cv_srv1o(1,i);
     cv_srv1(2,yrs_srv1(i))=like_wght_mbio*cv_srv1o(2,i);
     cv_srv1_nowt(1,yrs_srv1(i))=cv_srv1o(1,i);
     cv_srv1_nowt(2,yrs_srv1(i))=cv_srv1o(2,i);
    }
    // surv_like += norm2(elem_div( log(obs_srv1_spbiom(1)(yrs_srv1)+p_const )-log(biom_tmp(1)(yrs_srv1)+p_const ),sqrt(2)*sqrt(log(elem_prod(cv_srv1(1)(yrs_srv1),cv_srv1(1)(yrs_srv1))+1.0))));
    // surv_like += norm2(elem_div( log(obs_srv1_spbiom(2)(yrs_srv1)+p_const )-log(biom_tmp(2)(yrs_srv1)+p_const ),sqrt(2)*sqrt(log(elem_prod(cv_srv1(2)(yrs_srv1),cv_srv1(2)(yrs_srv1))+1.0))));
  for(i=1;i<=(nobs_srv1 - retro_years);i++)
  {
   surv_like += pow((log(obs_srv1_spbiom(1)(yrs_srv1(i))+p_const )-log(biom_tmp(1)(yrs_srv1(i))+p_const ))/(sqrt(2)*sqrt(log((cv_srv1(1)(yrs_srv1(i))*cv_srv1(1)(yrs_srv1(i)))+1.0))),2);
   surv_like += pow((log(obs_srv1_spbiom(2)(yrs_srv1(i))+p_const )-log(biom_tmp(2)(yrs_srv1(i))+p_const ))/(sqrt(2)*sqrt(log((cv_srv1(2)(yrs_srv1(i))*cv_srv1(2)(yrs_srv1(i)))+1.0))),2);
  }
    surv2_like = pow(( log(obs_srv2_spbiom(1,1)+p_const)-log(fspbio_srv2_ind+p_const))/ (sqrt(2)*sqrt(log((obs_srv2_cv(1,1)*obs_srv2_cv(1,1))+1.0))),2.0);
    surv2_like += extra_wght_ind_m *pow((log(obs_srv2_spbiom(1,2)+p_const)-log(mspbio_srv2_ind+p_const))/ (sqrt(2)*sqrt(log((obs_srv2_cv(1,2)*obs_srv2_cv(1,2))+1.0))),2.0);
    surv3_like = pow(( log(obs_srv2_spbiom(2,1)+p_const)-log(fspbio_srv2_nmfs+p_const))/ (sqrt(2)*sqrt(log((obs_srv2_cv(2,1)*obs_srv2_cv(2,1))+1.0))),2.0);
    surv3_like += pow(( log(obs_srv2_spbiom(2,2)+p_const)-log(mspbio_srv2_nmfs+p_const))/(sqrt(2)*sqrt(log((obs_srv2_cv(2,2)*obs_srv2_cv(2,2))+1.0))),2.0);
   f+=1.0*surv2_like;
   f+=1.0*surv3_like;			  	 
   Fout(14)=1.0*surv2_like;
   Fout(15)=1.0*surv3_like;
   surv10_like = pow(( log(obs_srv10_spbiom(1,1)+p_const)-log(fspbio_srv10_ind+ p_const))/(sqrt(2)*sqrt(log((obs_srv10_cv(1,1)*obs_srv10_cv(1,1))+1.0))),2.0);
   surv10_like += extra_wght_ind_m *pow((log(obs_srv10_spbiom(1,2)+p_const)-log(mspbio_srv10_ind+p_const))/(sqrt(2)*sqrt(log((obs_srv10_cv(1,2)*obs_srv10_cv(1,2))+1.0))),2.0);
   surv10nmfs_like = pow(( log(obs_srv10_spbiom(2,1)+p_const)-log(fspbio_srv10_nmfs+p_const))/(sqrt(2)*sqrt(log((obs_srv10_cv(2,1)*obs_srv10_cv(2,1))+1.0))),2.0);
   surv10nmfs_like += pow(( log(obs_srv10_spbiom(2,2)+p_const)-log(mspbio_srv10_nmfs+p_const))/(sqrt(2)*sqrt(log((obs_srv10_cv(2,2)*obs_srv10_cv(2,2))+1.0))),2.0);
   f+=1.0*surv10_like;
   f+=1.0*surv10nmfs_like;			  	 
   Fout(23)=1.0*surv10_like;
   Fout(24)=1.0*surv10nmfs_like;
   int yrc;
   cpue_pred = cpueq * legal_males;
   // calculate the observed large males
   obs_lmales.initialize();
   obs_lmales_bio.initialize();
   for(i=1;i<=nobs_srv1_length;i++)
    {
    //take 1/2 of the 100-104 bin, 
          obs_lmales(i)=0.5*obs_srv1_num(2,yrs_srv1_length(i),16);
          obs_lmales_bio(i)=obs_lmales(i)*wtm(16);
      for(j=17;j<=nlenm;j++)
        {
          obs_lmales(i)+=obs_srv1_num(2,yrs_srv1_length(i),j);
          obs_lmales_bio(i)+=obs_srv1_num(2,yrs_srv1_length(i),j)*wtm(j);
        }
     }
  //fishery cpue likelihood
    cpue_like=0.0;
   if(active(cpueq))
   { 
   for(yrc=styr;yrc<=(endyr-1- retro_years);yrc++)
    cpue_like += pow(((log(cpue(yrc+1)+1e-9)-log(cpue_pred(yrc)+1e-9))/(sqrt(2)* cpue_cv)),2.0);
    //cpue_like += norm2(elem_div( log(obs_srv1_spbiom(2)(yrs_srv1)+p_const )-log(biom_tmp(2)(yrs_srv1)+p_const ),sqrt(2)*sqrt(log(elem_prod(cv_srv1(2)(yrs_srv1),cv_srv1(2)(yrs_srv1))+1.0))));
    f+=cpue_like;
    Fout(16)=cpue_like;
   }
 //don't include last year as that would be endyr+1 fishery season
   catch_like1 = norm2((obs_catchtot_biom(1992,endyr-1-retro_years)-catch_ret(1992,endyr-1-retro_years)+p_const)-(pred_catch(1992,endyr-1-retro_years)-pred_catch_ret(1992,endyr-1- retro_years)+p_const));
    catch_likef = norm2((obs_catchdf_biom)-(pred_catch_disc(1)));
    catch_likef += smooth_disc_catch *norm2(first_difference(pred_catch_disc(1)));
    catch_like2 = norm2((catch_ret(styr,endyr-1-retro_years))-(pred_catch_ret(styr,endyr-1-retro_years)));
    catch_liket = norm2((obs_catcht_biom)-(pred_catch_trawl));
  f += like_wght(6)*1.*catch_like2;
  Fout(17)= like_wght(6)*1.*catch_like2;
  f += wght_total_catch*1.*catch_like1;
  Fout(18)= wght_total_catch*1.*catch_like1;
  f += like_wght(6)*catch_liket;
  Fout(19)= like_wght(6)*catch_liket;
  f += like_wght(6)*disc_catch_wght_fem *catch_likef;
  Fout(20)= like_wght(6)*disc_catch_wght_fem *catch_likef;
  f += surv_like;
  Fout(21)= surv_like;
  fpen.initialize();
 if(f_penalties>0)
 {
  if (current_phase()<2)
    fpen = like_wght_fph1*norm2(mfexp(fmort_dev+log_avg_fmort)-1.15);
  else
    fpen = like_wght_fph2*norm2(mfexp(fmort_dev+log_avg_fmort)-1.15);
 }
  if (active(fmortt_dev_era1))
  {
   fpen += norm2(fmortt_dev_era1(styr,1991));
   fpen += norm2(fmortt_dev_era2(1992,endyr-1));
  }
  if (active(fmort_dev))
  {
 // this keeps fmorts from going too high if like_wght_fdev is high...THIS ACTUALLY CONSTRAINS F TO ZERO
    fpen += like_wght_fdev*norm2(mfexp(fmort_dev(styr,endyr-1)+log_avg_fmort));
    fpen += 1.0*norm2(fmortdf_dev(styr,endyr-1));
   }
  Fout(22)= fpen;
  f+=fpen;
  call_no += 1;
  cout <<"Likes = "<< Fout << endl;
  cout <<"phase = "<< current_phase() << " call = " << call_no << " Total Like = " << f << endl;
}

void model_parameters::get_fut_mortality(void)
{
  ofstream& post= *pad_post;
 for(i=ipass;i<=endyr+Nproj-1;i++)
 {
	if(IsB0==1)
	{
      fmort(i)=0.0; 
	  fmortdf(i)=0.0;
      fmortt(i)=0.0;
	}
	else
	{
	fmort(i)=FutMort;  
	fmortdf(i)=mfexp(log_avg_fmortdf);
    fmortt(i)=mfexp(log_avg_fmortt);
	}	
 }
  for (i=ipass;i<=endyr+Nproj;i++)
 {
   if(i>=yrs_fish_discf(1) && i<endyr)
   {
    Fdiscf(i)=sel_discf(i)*fmortdf(i);
   }
   else
   {
    sel_discf_e=(1./(1.+mfexp(-1.*fish_disc_slope_f*(length_bins-mfexp(fish_disc_sel50_f)))));
    Fdiscf(i)=sel_discf_e*fmortdf(i);
    }
   Fdisct(1,i)=sel_trawl(1)*fmortt(i);   
   Fdisct(2,i)=sel_trawl(2)*fmortt(i);   
   for(k=1;k<=2;k++) //over new (k=1) and old (k=2) shell...
   { 
     F(k,i) = sel(k,i)*fmort(i);       
     F_ret(k,i)=sel_fit(k,i)*fmort(i);
     Fmat(k,i) = sel(k,i)*fmort(i);       
     Fmat_ret(k,i)=sel_fit(k,i)*fmort(i);
     Fimat(k,i) = sel(k,i)*fmort(i);      
     Fimat_ret(k,i)= sel_fit(k,i)*fmort(i);
     Ftot(1,k,i)=Fdiscf(i) + Fdisct(1,i);
     S(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Smat(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Simat(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Smat(2,k,i)=mfexp(-1.0*(Fmat(k,i)+Fdisct(2,i)));
     Simat(2,k,i)=mfexp(-1.0*(Fimat(k,i)+Fdisct(2,i)));
   } 
 }
  //===============================================================================
}

void model_parameters::Find_F35(void)
{
  ofstream& post= *pad_post;
 int icnt;
  IsB0 = 1;
  FutRec = 100000;
  FutMort = 0;
  ipass = endyr+1;
  get_fut_mortality();
  for (ipass=endyr+1;ipass<=endyr+Nproj-1;ipass++) 
	get_num_at_len_yr(); 
  Bzero =   mspbio_matetime(endyr+Nproj-2);  
  // Find F35%  
  Target = 0.35;
  IsB0 = 0;
  FutMort = 0.3;
  for (icnt=1;icnt<=20;icnt++)
   {
    ipass = endyr+1;
    get_fut_mortality();
   for (ipass=endyr+1;ipass<=endyr+Nproj-1;ipass++)
    get_num_at_len_yr(); 
    Btest = mspbio_matetime(endyr+Nproj-2) ;   
    Ratio = Btest/Bzero;
    cout << FutMort << " " << Ratio << endl;
    FutMort = FutMort * Ratio / Target;
   }
  F35 = FutMort; 
  SBPRF35 = Btest/FutRec;
}

void model_parameters::Find_OFL(void)
{
  ofstream& post= *pad_post;
  dvariable Fmsy,Rbar,nn,alpha,beta;
  int BMSY_Yr1, BMSY_Yr2,ii,Iyr,kk,jj;
 //Define time period for BMSY  
  BMSY_Yr1 = styr;BMSY_Yr2 = endyr-1;
  alpha = 0.1;
  beta = 0.25;
  Fmsy = F35;
 // Find Rbar 
  Rbar = 0; nn= 0;
  for (Iyr=BMSY_Yr1;Iyr<=BMSY_Yr2;Iyr++)
   {
    Rbar += mfexp(mean_log_rec(2)+rec_devm(Iyr));
    nn += 1;
   }
  Rbar = Rbar / nn;
 // Specify the BMSY proxy
  Bmsy = SBPRF35 * Rbar;
  //Begin projection
  // This code adjusts F so that the MMB at mate time next year is equal to the BMSY proxy
  // The resulting F is the FOFL
  ipass = endyr;
  if (ipass > endyr) 
	FutRec = Rbar;
  FutMort = Fmsy;
  get_fut_mortality();
  get_num_at_len_yr();
   if (mspbio_matetime(ipass) < Bmsy)
   {
     FutMort = 0;
     get_fut_mortality();
     get_num_at_len_yr();
     if (mspbio_matetime(ipass) > beta*Bmsy)
     {
       FutMort = Fmsy/2;
       get_fut_mortality();
       get_num_at_len_yr();
      for (ii=1;ii<=15;ii++)
       {
        FutMort = Fmsy*(mspbio_matetime(ipass)/Bmsy-alpha)/(1-alpha);
        get_fut_mortality();
        get_num_at_len_yr();
       }
      }
    }
   FOFL = FutMort;
   cout<<"Bzero"<<Bzero<<endl;     
   cout<<"SBPRF35"<<SBPRF35<<endl;   
   cout<<"FOFL"<<FOFL<<endl;
  //project under FOFL
   get_fut_mortality();
   get_num_at_len_yr();
   // find OFL from applying the FOFL
   get_catch_at_len();
   OFL = 0;
   OFL += pred_catch_ret(endyr);
   OFL += pred_catch_trawl(endyr);
   for(j=1;j<=2;j++)
    OFL += pred_catch_disc(j,endyr);
   cout<<"OFL"<<OFL<<endl; 
}

void model_parameters::Francis_weights(void)
{
  ofstream& post= *pad_post;
  int w,x,y,z;
  Lbar.initialize();
  Lbar_hat_new.initialize();
  Lbar_hat_old.initialize();
  SE_Lbar_hat_new.initialize();
  SE_Lbar_hat_old.initialize();
  Francis_var_temp_new.initialize();
  Francis_var_temp_old.initialize();
  // obs(maturity, SC, sex, year,length), pred(maturity,sex, year,length)
  // calculate observed and predicted average lengths by year (z) 
  for(w=1;w<=2;w++)
   for(x=1;x<=2;x++)
	for(y=1;y<=2;y++) 
     for(z=1;z<=nobs_srv1_length;z++)
	 {
	  Lbar(w,x,y,z) 		=sum(elem_prod(obs_p_srv1_len(w,x,y,z)/sum(obs_p_srv1_len(w,x,y,z)),length_bins));
	  Lbar_hat_new(x,y,z) 	=sum(elem_prod(pred_p_srv1_len_new(x,y,yrs_srv1_length(z))/sum(pred_p_srv1_len_new(x,y,yrs_srv1_length(z))),length_bins));
	  Lbar_hat_old(x,y,z) 	=sum(elem_prod(pred_p_srv1_len_old(x,y,yrs_srv1_length(z))/sum(pred_p_srv1_len_old(x,y,yrs_srv1_length(z))),length_bins));	
     }
  // calculate the predicted standard error of the mean length of the catch by year (z)	 
   for(x=1;x<=2;x++)
	for(y=1;y<=2;y++) 
     for(z=1;z<=nobs_srv1_length;z++)
	 {	
      SE_Lbar_hat_new(x,y,z) = sqrt(sum(elem_prod((pred_p_srv1_len_new(x,y,yrs_srv1_length(z))/sum(pred_p_srv1_len_new(x,y,yrs_srv1_length(z)))),pow((length_bins-Lbar_hat_new(x,y,z)),2))))/sqrt(nsamples_srv1_length(x,1,y,z)); 
      SE_Lbar_hat_old(x,y,z) = sqrt(sum(elem_prod((pred_p_srv1_len_old(x,y,yrs_srv1_length(z))/sum(pred_p_srv1_len_old(x,y,yrs_srv1_length(z)))),pow((length_bins-Lbar_hat_old(x,y,z)),2))))/sqrt(nsamples_srv1_length(x,2,y,z)); 
	  }
  // calculate the term for which the inverse variance is used to modify the input effective N for weighting len comps	
   for(x=1;x<=2;x++)
	for(y=1;y<=2;y++) 
     for(z=1;z<=nobs_srv1_length;z++)
	 {	
	  Francis_var_temp_new(x,y,z) = (Lbar(x,1,y,z)-Lbar_hat_new(x,y,z))/SE_Lbar_hat_new(x,y,z);
	  Francis_var_temp_old(x,y,z) = (Lbar(x,2,y,z)-Lbar_hat_old(x,y,z))/SE_Lbar_hat_old(x,y,z);
	 }
   // find the mena for females and males separately to calculate the variance
   // why isn't there a 'var()' function?! ugh.
   countFem = 0;
   countMal = 0;
   totalFem = 0;
   totalMal = 0;
  for(x=1;x<=2;x++)
   for(z=1;z<=nobs_srv1_length;z++)
   {
    // the if statements are because some of these are 'NaN' e.g. immature old shell crab
	// should find a better way to do this
	//females
	 if(Francis_var_temp_new(x,1,z)>0)
	 {
	  totalFem += Francis_var_temp_new(x,1,z);
	  countFem += 1;
	 }
	 if(Francis_var_temp_old(x,1,z)>0)
	 {
	  totalFem += Francis_var_temp_old(x,1,z);
	  countFem += 1;
	 }
	//males
	 if( Francis_var_temp_new(x,2,z)>0)
	 {
	  totalMal +=  Francis_var_temp_new(x,2,z);
	  countMal += 1;
	 }
 	 if( Francis_var_temp_old(x,2,z)>0)
	 {
	  totalMal +=  Francis_var_temp_old(x,2,z);
	  countMal += 1;
	 }
    }
	FemMeanVarTerm = 0.0;
	MaleMeanVarTerm = 0.0;
  // find the variance 
  for(x=1;x<=2;x++)
   for(z=1;z<=nobs_srv1_length;z++)
   {
	if(Francis_var_temp_new(x,1,z)>0)
     FemMeanVarTerm += square(totalFem/countFem - Francis_var_temp_new(x,1,z));
 	if(Francis_var_temp_old(x,1,z)>0)
	 FemMeanVarTerm += square(totalFem/countFem - Francis_var_temp_old(x,1,z));
 	if(Francis_var_temp_new(x,2,z)>0)
	 MaleMeanVarTerm += square(totalMal/countMal - Francis_var_temp_new(x,2,z));
 	if(Francis_var_temp_old(x,2,z)>0)
	 MaleMeanVarTerm += square(totalMal/countMal - Francis_var_temp_old(x,2,z));
   }
   Francis_weight_m = 1/(MaleMeanVarTerm/(countMal-1));
   Francis_weight_f = 1/(FemMeanVarTerm/(countFem-1));
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  Find_F35();
  Find_OFL();
  Francis_weights();
  report<<f<<" "<<Bmsy<<" "<< mspbio_srv1(endyr) << " " << mspbio_fishtime(endyr)<< " "<< mspbio_matetime(endyr-1)<< " " << F35 << " " << FOFL << " " << OFL << endl;
  report<<Fout<< " "<<endl;
  report<<"# mature male biomass at survey"<<endl;
  report<<mspbio_srv1<<endl;
  report<<"# Francis weight F"<<endl;
  report<<Francis_weight_f<<endl;
  report<<"# Francis weight M"<<endl;
  report<<Francis_weight_m<<endl;
  // this writes gradients for each parameter at each phase and at end "gradients.dat"
   save_gradients(gradients);
   //The report file holds all of the quantities used in Jack's projection program
   //The "Rout" below this writes to a file that is read for plotting in the safe document
  report<<"#number of length bins"<<endl;
  report<<nlenm<<endl;
  report<<"#Nat mort immature female/male"<<endl;
  report<<M<<endl;
  report<<"#nat mort mature new shell female/male"<<endl;
  report<<M_matn<<endl;
  report<<"#nat mort mature old shell female/male"<<endl;
  report<<M_mato<<endl;
  report<<"#constant recruitment"<<endl;
  report<<"1000000"<<endl;
  report<<"#average of last 4 years sel total male new old shell"<<endl;
  report<<(sel(1,endyr-4)+sel(1,endyr-3)+sel(1,endyr-2)+sel(1,endyr-1))/4.0<<endl;
  report<<(sel(1,endyr-4)+sel(2,endyr-3)+sel(2,endyr-2)+sel(2,endyr-1))/4.0<<endl;
  report<<"#average of last 4 years sel retained curve male new old shell"<<endl;
  report<<(sel_fit(1,endyr-3)+sel_fit(1,endyr-2)+sel_fit(1,endyr-1))/3.0<<endl;
  report<<(sel_fit(2,endyr-3)+sel_fit(2,endyr-2)+sel_fit(2,endyr-1))/3.0<<endl;
  report<<"#trawl selectivity female male"<<endl;
  report<<sel_trawl<<endl;
  report<<"#female pot discard selectivity"<<endl;
  report<<sel_discf(2009)<<endl;
  report<<"#maturity curve new shell female male"<<endl;
  report<<maturity_est(1)<<endl;
  report<<maturity_est(2)<<endl;
  report<<"#maturity curve old shell female male"<<endl;
  report<<maturity_old_average<<endl;
  report<<"#molting probability immature female male"<<endl;
  report<<moltp<<endl;
  report<<"#molting probability mature female male"<<endl;
  report<<moltp_mat<<endl;
  report<<"#prop recruits to new shell"<<endl;
  report<<proprecn<<endl;
  report<<"#distribution of recruits to length bins"<<endl;
  report<<rec_len<<endl;
  report<<"#time of catch in fraction of year from survey - 7 months"<<endl;
  report<<catch_midpt(endyr)<<endl;
  report<<"#number at length new shell females males at time of fishery endyr from model"<<endl;
  report<<natl_new_fishtime(1,endyr)<<endl;
  report<<natl_new_fishtime(2,endyr)<<endl;
  report<<"#number at length old shell females males at time of fishery endyr from model"<<endl;
  report<<natl_old_fishtime(1,endyr)<<endl;
  report<<natl_old_fishtime(2,endyr)<<endl;
  report<<"#last year male spawning biomass"<<endl;
  report<<mspbio(endyr)<<endl;
  report<<"#last year female spawning biomass"<<endl;
  report<<fspbio(endyr)<<endl;
  report<<"#last year male spawning biomass at matingtime"<<endl;
  report<<mspbio_matetime(endyr-1)<<endl;
  report<<"#last year female spawning biomass at matingtime"<<endl;
  report<<fspbio_matetime(endyr-1)<<endl;
  report<<"#numbers at length immature new shell female male last year"<<endl;
  report<<natlength_inew(1,endyr)<<endl;
  report<<natlength_inew(2,endyr)<<endl;
  report<<"#numbers at length immature old shell female male last year"<<endl;
  report<<natlength_iold(1,endyr)<<endl;
  report<<natlength_iold(2,endyr)<<endl;
  report<<"#numbers at length mature new shell female male last year"<<endl;
  report<<natlength_mnew(1,endyr)<<endl;
  report<<natlength_mnew(2,endyr)<<endl;
  report<<"#numbers at length mature old shell female male last year"<<endl;
  report<<natlength_mold(1,endyr)<<endl;
  report<<natlength_mold(2,endyr)<<endl;
  report<<"#weight at length female juvenile"<<endl;
  report<<wtf(1)<<endl;
  report<<"#weight at length female mature"<<endl;
  report<<wtf(2)<<endl;
  report<<"#weight at length male"<<endl;
  report<<wtm<<endl;
  report<<"#length-length transition matrix"<<endl;
  report<<len_len<<endl;
  report<<"#female discard pot fishing F"<<endl;
  report<<fmortdf(endyr-1)<<endl;
  report<<"#trawl fishing F female male"<<endl;
  report<<fmortt(endyr-1)<<endl;
  report<<"#number of recruits from the model styr to endyr-1"<<endl;
  report<<endyr-styr<<endl;
  report <<"#recruitments female, male start year to endyr-1 from model" << endl;
  for(i=styr; i<endyr; i++)
  {
    report << mfexp(mean_log_rec(1)+rec_dev(1,i))<<" ";
  }
  report <<endl<< "#recruitments male, male start year+1 to endyr-1 from model" << endl;
  for(i=styr; i<endyr; i++)
  {
    report << mfexp(mean_log_rec(2)+rec_dev(2,i))<<" ";
  }
  report<<endl;
  report<<"#male spawning biomass at matetime for endyr-5 to endyr-1 for spawner recruit curve to estimate recruitments"<<endl;
  report<<mspbio_matetime(endyr-5,endyr-1)<<endl;
  report<<"#male spawning biomass at matetime for str year to endyr-1 for spawner recruit curve to estimate recruitments"<<endl;
  report<<mspbio_matetime(styr,endyr-1)<<endl;
  report <<"#selectivity survey males 1989 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  report << sel_srv3(2) << endl;
    if (last_phase())
  {
    R_out << "$Estimated numbers of immature new shell female crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out <<  i<<" "<<natlength_inew(1,i) << endl;
       }
  R_out << "$Estimated numbers of immature old shell female crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out <<  i<<" "<<natlength_iold(1,i) << endl;
       }
  R_out << "$Estimated numbers of mature new shell female crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out <<  i<<" "<<natlength_mnew(1,i) << endl;
       }
  R_out << "$Estimated numbers of mature old shell female crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out <<  i<<" "<<natlength_mold(1,i) << endl;
       }
  R_out << "$Estimated numbers of immature new shell male crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
        R_out << i<<" "<<natlength_inew(2,i) << endl;
      }
  R_out << "$Estimated numbers of immature old shell male crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
        R_out << i<<" "<<natlength_iold(2,i) << endl;
      }
  R_out << "$Estimated numbers of mature new shell male crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
        R_out << i<<" "<<natlength_mnew(2,i) << endl;
      }
 //  R_out << "Estimated numbers of mature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Estimated numbers of mature old shell male crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
        R_out << i<<" "<<natlength_mold(2,i) << endl;
      }
 //  R_out << "Observed numbers of immature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of immature new shell female crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(1,1,1,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of mature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of mature new shell female crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(2,1,1,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of mature old shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of mature old shell female crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(2,2,1,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of immature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of immature new shell male crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(1,1,2,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of immature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of immature old shell male crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(1,2,2,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of mature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of mature new shell male crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(2,1,2,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of mature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of mature old shell male crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(2,2,2,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
  //  R_out << "Observed Survey Numbers by length females:  'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed survey numbers female" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
        R_out<<yrs_srv1_length(i)<<" " << obs_srv1_num(1,yrs_srv1_length(i)) << endl;
      }
  //  R_out << "Observed Survey Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed survey numbers male" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
        R_out<<yrs_srv1_length(i)<<" " << obs_srv1_num(2,yrs_srv1_length(i))<< endl;
      }
  //  R_out << "Predicted Survey Numbers by length females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted survey numbers female" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
       R_out<<yrs_srv1_length(i)<<" "  << pred_srv1(1,yrs_srv1_length(i)) << endl;
      }
  //  R_out << "Predicted Survey Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted survey numbers male" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
         R_out<<yrs_srv1_length(i)<<" "  << pred_srv1(2,yrs_srv1_length(i)) << endl;
       }
  //  R_out << "Predicted pop Numbers by length females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted population numbers female" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out<<i<<" "<< natlength(1,i)<< endl;
      }
  //  R_out << "Predicted pop Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted population numbers male" << endl;
     for(i=styr;i<=endyr;i++)
      {
         R_out<<i<<" "<< natlength(2,i)<< endl;
       }
  //  R_out<<"observed number of males greater than 101 mm: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed number males greater than 101mm" << endl;
  R_out<<obs_lmales<<endl;
  //  R_out<<"observed biomass of males greater than 101 mm: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed biomass males greater than 101mm" << endl;
  R_out<<obs_lmales_bio<<endl;
  //  R_out<<"Predicted number of males greater than 101 mm: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted number males greater than 101mm" << endl;
  R_out<<num_males_gt101<<endl;
  //  R_out<<"Predicted biomass of males greater than 101 mm: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted biomass males greater than 101mm" << endl;
  R_out<<bio_males_gt101<<endl;
  //  R_out<<"pop estimate numbers of males >101: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Population numbers male" << endl;
        R_out<<legal_males<<endl;
      //  R_out<<"estimated population biomass of males > 101: seq(1978,"<<endyr<<") "<<endl;
  R_out << "$Population biomass male" << endl;
      R_out<<legal_males_bio<<endl;
      //  R_out<<"estimated survey numbers of males > 101: seq(1978,"<<endyr<<") "<<endl;
  R_out << "$Estimated survey numbers male" << endl;
      R_out<<legal_srv_males<<endl;
      //  R_out<<"estimated survey biomass of males > 101: seq(1978,"<<endyr<<") "<<endl;
  R_out << "$Estimated survey biomass male" << endl;
      R_out<<legal_srv_males_bio<<endl;
  //  R_out << "Observed survey biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey biomass" << endl;
  R_out << obs_srv1_biom(styr,endyr)<<endl;
  //  R_out << "predicted survey biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted survey biomass" << endl;
  R_out << pred_srv1_bioms(1)+pred_srv1_bioms(2)<<endl;
  //survey numbers
    for(k=1;k<=2;k++)
    {
     for(i=styr;i<=endyr;i++)
      {
       tmpo(k,i)=sum(obs_srv1_num(k,i));
       tmpp(k,i)=sum(pred_srv1(k,i));
      }
     }
  //  R_out << "Observed survey numbers female: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey numbers female" << endl;
  R_out << tmpo(1)<<endl;
  //  R_out << "Observed survey numbers male: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey numbers male" << endl;
  R_out << tmpo(2)<<endl;
  //  R_out << "predicted survey numbers female: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted survey numbers female" << endl;
  R_out << tmpp(1)<<endl;
  //  R_out << "predicted survey numbers male: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted survey numbers male" << endl;
  R_out << tmpp(2)<<endl;
  //  R_out << "Observed survey female spawning biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey female spawning biomass" << endl;
  R_out << obs_srv1_spbiom(1)<<endl;
  //  R_out << "Observed survey male spawning biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey male spawning biomass" << endl;
  R_out << obs_srv1_spbiom(2)<<endl;
  //  R_out << "Observed survey female new spawning numbers: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$bserved survey female new spawning numbers" << endl;
  R_out << obs_srv1_spnum(1,1)<<endl;
  //  R_out << "Observed survey female old spawning numbers: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey female old spawning numbers" << endl;
  R_out << obs_srv1_spnum(2,1)<<endl;
  //  R_out << "Observed survey male new spawning numbers: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey male new spawning numbers" << endl;
  R_out << obs_srv1_spnum(1,2)<<endl;
  //  R_out << "Observed survey male old spawning numbers: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey male old spawning numbers" << endl;
  R_out << obs_srv1_spnum(2,2)<<endl;
  //  R_out << "Observed survey female biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey female biomass" << endl;
  R_out << obs_srv1_bioms(1)<<endl;
  //  R_out << "Observed survey male biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey male biomass" << endl;
  R_out << obs_srv1_bioms(2)<<endl;
  //  R_out << "natural mortality immature females, males: 'FemM','MaleM'" << endl;
  R_out << "$natural mortality immature" << endl;
  R_out << M << endl;
  //  R_out << "natural mortality mature females, males: 'FemMm','MaleMm'" << endl;
  R_out << "$atural mortality mature" << endl;
  R_out << M_matn << endl;
  //  R_out << "natural mortality mature old shell females, males: 'FemMmo','MaleMmo'" << endl;
  R_out << "$natural mortality mature old shell" << endl;
  R_out << M_mato << endl;
  //  R_out << "Predicted Biomass: seq(1978,"<<endyr<<")" << endl;
  R_out << "$Predicted Biomass" << endl;
  R_out << pred_bio << endl;
  //  R_out << "Predicted total population numbers: seq(1978,"<<endyr<<") "<<endl;
  R_out << "$Predicted total population numbers" << endl;
  R_out <<popn<<endl;
  R_out << "$Female Spawning Biomass" << endl;
  R_out << fspbio << endl;
  //  R_out << "Male Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mature Male Biomass" << endl;
  R_out << mspbio << endl;
  //  R_out << "Total Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Total Spawning Biomass" << endl;
  R_out << fspbio+mspbio << endl;
  //  R_out << "Female Spawning Biomass at fish time: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Female Spawning Biomass at fish time" << endl;
  R_out << fspbio_fishtime << endl;
  //  R_out << "Male Spawning Biomass at fish time: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mature Male Biomass at fish time" << endl;
  R_out << mspbio_fishtime << endl;
  //  R_out << "Total Spawning Biomass at fish time: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Total Spawning Biomass at fish time" << endl;
  R_out << fspbio_fishtime+mspbio_fishtime << endl;
  //  R_out << "Mating time Female Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mating time Female Spawning Biomass" << endl;
  R_out << fspbio_matetime << endl;
  //  R_out << "Mating time Male Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mature male biomass at mating" << endl;
  R_out << mspbio_matetime << endl;
  //  R_out << "Mating time Male old shell Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Old shell mature male biomass at mating" << endl;
  R_out << mspbio_old_matetime << endl;
  //  R_out << "Mating time female new shell Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$New shell female spawning biomass at mating time" << endl;
  R_out << fspbio_new_matetime << endl;
  //  R_out << "Mating time Total Spawning Biomass : seq(1978,"<<endyr<<") " << endl;
  R_out << "$Total mature biomass at mating time" << endl;
  R_out << fspbio_matetime+mspbio_matetime << endl;
  //  R_out << "Mating time effective Female Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Effective female spawning biomass at mating time" << endl;
  R_out << efspbio_matetime << endl;
  //  R_out << "Mating time effective Male Spawning Biomass(old shell only): seq(1978,"<<endyr<<") " << endl;
  R_out << "$Effective mature male biomass at mating time" << endl;
  R_out << emspbio_matetime << endl;
  //  R_out << "Mating time Total effective Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Total effective spawning biomass at mating time" << endl;
  R_out << efspbio_matetime+emspbio_matetime << endl;
  //  R_out << "Mating time male Spawning numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Spawning numbers male at mating time" << endl;
  R_out << mspnum_matetime << endl;
  //  R_out << "Mating time Female Spawning numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Spawning number female at mating time" << endl;
  R_out << efspnum_matetime << endl;
  //  R_out << "Mating time Male Spawning numbers(old shell only): seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mating time Male Spawning numbers old " << endl;
  R_out << emspnum_old_matetime << endl;
  //  R_out << "ratio Mating time Female Spawning numbers to male old shell mature numbers : seq(1978,"<<endyr<<") " << endl;
  R_out << "$ratio Mating time Female Spawning numbers to male old shell mature numbers" << endl;
  R_out << elem_div(efspnum_matetime,emspnum_old_matetime) << endl;
  //  R_out << "Mating time effective Female new shell Spawning biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mating time effective Female new shell Spawning biomass" << endl;
  R_out <<efspbio_new_matetime << endl;
  //  R_out << "Mating time Female new shell Spawning numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mating time Female new shell Spawning numbers" << endl;
  R_out << fspnum_new_matetime << endl;
  //  R_out << "ratio Mating time Female new shell Spawning numbers to male old shell mature numbers : seq(1978,"<<endyr<<") " << endl;
  R_out << "$ratio Mating time Female new shell Spawning numbers to male old shell mature numbers" << endl;
  R_out << elem_div(fspnum_new_matetime,emspnum_old_matetime) << endl;
  //  R_out << "Predicted Female survey Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Female survey Biomass" << endl;
  R_out << pred_srv1_bioms(1) << endl;
  //  R_out << "Predicted Male survey Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Male survey Biomass" << endl;
  R_out << pred_srv1_bioms(2)<< endl;
  //  R_out << "Predicted Female survey mature Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Female survey mature Biomass" << endl;
  R_out << fspbio_srv1 << endl;
  //  R_out << "Predicted Male survey mature Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Male survey mature Biomass" << endl;
  R_out << mspbio_srv1<< endl;
  //  R_out << "Predicted total survey mature Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted total survey mature Biomass" << endl;
  R_out << fspbio_srv1+mspbio_srv1<< endl;
  //  R_out << "Predicted Female survey new mature numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Female survey new mature numbers" << endl;
  R_out << fspbio_srv1_num(1) << endl;
  //  R_out << "Predicted Female survey old mature numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$redicted Female survey old mature numbers" << endl;
  R_out << fspbio_srv1_num(2) << endl;
  //  R_out << "Predicted Male survey new mature numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Male survey new mature numbers" << endl;
  R_out << mspbio_srv1_num(1)<< endl;
  //  R_out << "Predicted Male survey old mature numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Male survey old mature numbers" << endl;
  R_out << mspbio_srv1_num(2)<< endl;
  //  R_out << "Observed industry survey mature biomass: seq(1,4) " << endl;
  R_out << "$Observed industry survey mature biomass" << endl;
  R_out << obs_srv2_spbiom(1,1)<<" "<<obs_srv2_spbiom(1,2)<<" "<<obs_srv2_spbiom(2,1)<<" "<<obs_srv2_spbiom(2,2)<<endl;
  //  R_out << "Predicted industry survey mature biomass: seq(1,4) " << endl;
  R_out << "$Predicted industry survey mature biomass" << endl;
  R_out << fspbio_srv2_ind<<" "<<mspbio_srv2_ind<<" "<<fspbio_srv2_nmfs<<" "<<mspbio_srv2_nmfs<<endl;
  //  R_out << "Observed Prop industry survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$bserved Prop industry survey female" << endl;
  R_out <<(obs_p_srv2_len(1,1,1)+obs_p_srv2_len(1,1,2))<<endl;
  //  R_out << "Observed Prop industry survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop industry survey male" << endl;
  R_out <<(obs_p_srv2_len(1,2,1)+obs_p_srv2_len(1,2,2))<<endl;
  //  R_out << "Observed Prop industry nmfs survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop industry nmfs survey female" << endl;
  R_out <<obs_p_srv2_len(2,1,1)+obs_p_srv2_len(2,1,2)<<endl;
  //  R_out << "Observed Prop industry nmfs survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop industry nmfs survey male" << endl;
  R_out <<obs_p_srv2_len(2,2,1)+obs_p_srv2_len(2,2,2)<<endl;
  //  R_out << "Predicted Prop industry survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop industry survey female" << endl;
  R_out <<pred_p_srv2_len_ind(1)<<endl;
  //  R_out << "Predicted Prop industry survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop industry survey male" << endl;
  R_out <<pred_p_srv2_len_ind(2)<<endl;
  //  R_out << "Predicted Prop industry nmfs survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop industry nmfs survey female" << endl;
  R_out <<pred_p_srv2_len_nmfs(1)<<endl;
  //  R_out << "Predicted Prop industry nmfs survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$redicted Prop industry nmfs survey male" << endl;
  R_out <<pred_p_srv2_len_nmfs(2)<<endl;
  //  R_out << "Observed 2010 industry survey mature biomass: seq(1,4) " << endl;
  R_out << "$Observed 2010 industry survey mature biomass" << endl;
  R_out << obs_srv10_spbiom(1,1)<<" "<<obs_srv10_spbiom(1,2)<<" "<<obs_srv10_spbiom(2,1)<<" "<<obs_srv10_spbiom(2,2)<<endl;
  //  R_out << "Predicted 2010 industry survey mature biomass: seq(1,4) " << endl;
  R_out << "$Predicted 2010 industry survey mature biomass" << endl;
  R_out << fspbio_srv10_ind<<" "<<mspbio_srv10_ind<<" "<<fspbio_srv10_nmfs<<" "<<mspbio_srv10_nmfs<<endl;
  //  R_out << "Observed Prop 2010 industry survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop 2010 industry survey female" << endl;
  R_out <<(obs_p_srv10_len(1,1,1)+obs_p_srv10_len(1,1,2))<<endl;
  //  R_out << "Observed Prop 2010 industry survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop 2010 industry survey male" << endl;
  R_out <<(obs_p_srv10_len(1,2,1)+obs_p_srv10_len(1,2,2))<<endl;
  //  R_out << "Observed Prop 2010 industry nmfs survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop 2010 industry nmfs survey female" << endl;
  R_out <<obs_p_srv10_len(2,1,1)+obs_p_srv10_len(2,1,2)<<endl;
  //  R_out << "Observed Prop 2010 industry nmfs survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop 2010 industry nmfs survey male" << endl;
  R_out <<obs_p_srv10_len(2,2,1)+obs_p_srv10_len(2,2,2)<<endl;
  //  R_out << "Predicted Prop 2010 industry survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop 2010 industry survey female" << endl;
  R_out <<pred_p_srv10_len_ind(1)<<endl;
  //  R_out << "Predicted Prop 2010 industry survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop 2010 industry survey male" << endl;
  R_out <<pred_p_srv10_len_ind(2)<<endl;
  //  R_out << "Predicted Prop 2010 industry nmfs survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop 2010 industry nmfs survey female" << endl;
  R_out <<pred_p_srv10_len_nmfs(1)<<endl;
  //  R_out << "Predicted Prop 2010 industry nmfs survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop 2010 industry nmfs survey male" << endl;
  R_out <<pred_p_srv10_len_nmfs(2)<<endl;
   //  R_out << "Observed Prop fishery ret new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery retained new male" << endl;
      for (i=1; i<=nobs_fish; i++)
        {
          R_out << yrs_fish(i) << " " << obs_p_fish_ret(1,i)<< endl;
         }
    //  R_out << "Predicted length prop fishery ret new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery retained new male" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish; i++)
        {
          ii=yrs_fish(i);  
          R_out <<  ii  <<  " "  <<  pred_p_fish_fit(1,ii)  << endl;
         }
  //  R_out << "Observed Prop fishery ret old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery retained old male" << endl;
      for (i=1; i<=nobs_fish; i++)
        {
          R_out << yrs_fish(i) << " " << obs_p_fish_ret(2,i)<< endl;
         }
    //  R_out << "Predicted length prop fishery ret old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery retained old male" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish; i++)
        {
          ii=yrs_fish(i);  
          R_out <<  ii  <<  " "  <<  pred_p_fish_fit(2,ii)  << endl;
         }
  //  R_out << "Observed Prop fishery total new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery total new male" << endl;
      for (i=1; i<=nobs_fish_discm; i++)
        {
          R_out << yrs_fish_discm(i) << " " << obs_p_fish_tot(1,i) << endl;
         }
    //  R_out << "Predicted length prop fishery total new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery total new male" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish_discm; i++)
        {
          ii=yrs_fish_discm(i);  
          R_out <<  ii  <<  " "  <<  pred_p_fish(1,ii)  << endl;
         }
  //  R_out << "Observed Prop fishery total old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery total old male" << endl;
      for (i=1; i<=nobs_fish_discm; i++)
        {
          R_out << yrs_fish_discm(i) << " " << obs_p_fish_tot(2,i) << endl;
         }
    //  R_out << "Predicted length prop fishery total old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery total old male" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish_discm; i++)
        {
          ii=yrs_fish_discm(i);  
          R_out <<  ii  <<  " "  <<  pred_p_fish(2,ii)  << endl;
         }
  //  R_out << "Observed Prop fishery discard new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery discard new male" << endl;
      for (i=1; i<=nobs_fish_discm; i++)
        {
          R_out << yrs_fish_discm(i) << " " << obs_p_fish_discm(1,i) << endl;
         }
 // R_out << "Observed Prop fishery discard old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery discard old male" << endl;
      for (i=1; i<=nobs_fish_discm; i++)
        {
          R_out << yrs_fish_discm(i) << " " << obs_p_fish_discm(2,i)<< endl;
         }
    //  R_out << "Observed length prop fishery discard all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion fishery discard all female" << endl;
    for (i=1; i<=nobs_fish_discf; i++)
    {
      R_out <<  yrs_fish_discf(i)  <<  " "  <<  obs_p_fish_discf(i)  << endl;
    }
    //  R_out << "Predicted length prop fishery discard all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery discard all female" << endl;
    for (i=1; i<=nobs_fish_discf; i++)
    {
      ii=yrs_fish_discf(i);  
      R_out <<  ii  <<  " "  <<  pred_p_fish_discf(ii)  << endl;
    }
    //  R_out << "Predicted length prop trawl females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion trawl female" << endl;
       for (i=1; i<=nobs_trawl; i++)
        {
          ii=yrs_trawl(i);  
          R_out <<  ii  <<  " "  <<  pred_p_trawl(1,ii)  << endl;
         }
    //  R_out << "Observed length prop trawl females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion trawl female" << endl;
       for (i=1; i<=nobs_trawl; i++)
        {  
           R_out <<  yrs_trawl(i)  <<  " "  <<  obs_p_trawl(1,i)  << endl;
         }
    //  R_out << "Predicted length prop trawl males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion trawl male" << endl;
       for (i=1; i<=nobs_trawl; i++)
        {
          ii=yrs_trawl(i);  
          R_out <<  ii  <<  " "  <<  pred_p_trawl(2,ii)  << endl;
         }
    //  R_out << "Observed length prop trawl males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion trawl male" << endl;
       for (i=1; i<=nobs_trawl; i++)
        {  
           R_out <<  yrs_trawl(i)  <<  " "  <<  obs_p_trawl(2,i)  << endl;
         }
  //    for (i=1; i<=nobs_fish; i++)
  //      {
  //        R_out << yrs_fish(i) << " " << obs_p_fish(2,i) << endl;
  //       }
   //R_out<<pred_p_fish<<endl;
   //    for (i=1; i<=nobs_fish; i++)
  //      {
    //      ii=yrs_fish(i);  
    //      R_out <<  ii  <<  " "  <<  pred_p_fish(2,ii)  << endl;
    //     }
 //  R_out << "Observed Length Prop survey immature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey immature new female" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(1,1,1,i) << endl;
           }
  //  R_out << "Predicted length prop survey immature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey immature new female" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_new(1,1,ii) << endl;
          }
  //  R_out << "Observed Length Prop survey immature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey immature old female" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(1,2,1,i) << endl;
           }
  //  R_out << "Predicted length prop survey immature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey immature old female" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_old(1,1,ii) << endl;
          }
  //  R_out << "Observed Length Prop survey immature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey immature new male" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             R_out << yrs_srv1_length(i) <<" " <<obs_p_srv1_len(1,1,2,i) << endl;
          }
  //  R_out << "Predicted length prop survey immature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey immature new male" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_new(1,2,ii) << endl;
         }
 //  R_out << "Observed Length Prop survey immature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey immature old male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   R_out << yrs_srv1_length(i) <<" " <<obs_p_srv1_len(1,2,2,i) << endl;
 }
 //  R_out << "Predicted length prop survey immature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey immature old male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   ii=yrs_srv1_length(i);  
   R_out << ii << " " << pred_p_srv1_len_old(1,2,ii) << endl;
 }
 //  R_out << "Observed Length Prop survey mature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey mature new female" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(2,1,1,i) << endl;
           }
  //  R_out << "Predicted length prop survey mature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey mature new female" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_new(2,1,ii) << endl;
          }
  //  R_out << "Observed Length Prop survey mature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey mature old female" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(2,2,1,i) << endl;
           }
  //  R_out << "Predicted length prop survey mature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey mature old female" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_old(2,1,ii) << endl;
          }
  //  R_out << "Observed Length Prop survey mature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey mature new male" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             R_out << yrs_srv1_length(i) <<" " <<obs_p_srv1_len(2,1,2,i) << endl;
          }
  //  R_out << "Predicted length prop survey mature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey mature new male" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_new(2,2,ii) << endl;
         }
 //  R_out << "Observed Length Prop survey mature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey mature old male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   R_out << yrs_srv1_length(i) <<" " <<obs_p_srv1_len(2,2,2,i) << endl;
 }
 //  R_out << "Predicted length prop survey mature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey mature old male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   ii=yrs_srv1_length(i);  
   R_out << ii << " " << pred_p_srv1_len_old(2,2,ii) << endl;
 }
 //  R_out << "Observed Length Prop survey all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey all females" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(1,1,1,i)+obs_p_srv1_len(2,1,1,i)+obs_p_srv1_len(2,2,1,i)<< endl;
              tmpp4+=obs_p_srv1_len(1,1,1,i)+obs_p_srv1_len(2,1,1,i)+obs_p_srv1_len(2,2,1,i);
                         }
 //  R_out << "Predicted length prop survey all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey all female" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   ii=yrs_srv1_length(i);  
   R_out << ii << " " << pred_p_srv1_len_new(1,1,ii)+pred_p_srv1_len_new(2,1,ii)+pred_p_srv1_len_old(2,1,ii) << endl;
    tmpp1+=pred_p_srv1_len_new(1,1,ii)+pred_p_srv1_len_new(2,1,ii)+pred_p_srv1_len_old(2,1,ii);
    }
 //  R_out << "Observed Length Prop survey all males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey all male" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(1,1,2,i)+obs_p_srv1_len(1,2,2,i)+obs_p_srv1_len(2,1,2,i)+obs_p_srv1_len(2,2,2,i)<< endl;
         tmpp2+=obs_p_srv1_len(1,1,2,i)+obs_p_srv1_len(1,2,2,i)+obs_p_srv1_len(2,1,2,i)+obs_p_srv1_len(2,2,2,i);
                        }
 //  R_out << "Predicted length prop survey all males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey all male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   ii=yrs_srv1_length(i);  
   R_out << ii << " " << pred_p_srv1_len_new(1,2,ii)+pred_p_srv1_len_new(2,2,ii)+pred_p_srv1_len_old(2,2,ii) << endl;
  tmpp3+=pred_p_srv1_len_new(1,2,ii)+pred_p_srv1_len_new(2,2,ii)+pred_p_srv1_len_old(2,2,ii);
    }
  //  R_out << "Sum of predicted prop survey all females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Sum of predicted proportion survey all female" << endl;
            R_out <<tmpp1<<endl;
  //  R_out << "Sum of predicted prop survey all males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Sum of predicted proportion survey all male" << endl;
            R_out <<tmpp3<<endl;
  //  R_out << "Sum of Observed prop survey all females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Sum of observed proportion survey all female" << endl;
            R_out <<tmpp4<<endl;
  //  R_out << "Sum of Observed prop survey all males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Sum of observed proportion survey all male" << endl;
            R_out <<tmpp2<<endl;
  //  R_out << "Predicted mean postmolt length females:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted mean postmolt length female" << endl;
  R_out << mean_length(1) << endl;
  //  R_out << "Predicted mean postmolt length males:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted mean postmolt length male" << endl;
  R_out << mean_length(2)<<endl; 
  //  R_out<< "sd mean length females:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<<endl;
  R_out << "$sd mean length female" << endl;
  R_out<<sd_mean_length(1)<<endl;
  //  R_out<< "sd mean length males:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<<endl;
  R_out << "$sd mean length male" << endl;
  R_out<<sd_mean_length(2)<<endl;
  //  R_out << "af: 'females'" << endl;
  R_out << "$af " << endl;
  R_out << af << endl;
  //  R_out << "am: 'males'" << endl;
  R_out << "$am" << endl;
  R_out << am << endl;
  R_out << "$af2" << endl;
  R_out << af +(bf-bf1)*deltaf << endl;
  //  R_out << "bf: 'females'" << endl;
  R_out << "$bf " << endl;
  R_out << bf << endl;
  R_out << "$bf1" << endl;
  R_out << bf1 << endl;
  R_out << "$EN4" << endl;
    R_out << deltaf << endl;
  R_out << "$bm " << endl;
  R_out << bm << endl;
  R_out << "$st_gr" << endl;
  R_out << st_gr << endl;
  R_out << "$st_gr" << endl;
  R_out << st_gr << endl;
    //  R_out << "Predicted probability of maturing females:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted probability of maturing female" << endl;
  R_out << maturity_est(1)<<endl; 
  //  R_out << "Predicted probability of maturing males:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted probability of maturing male" << endl;
  R_out << maturity_est(2)<<endl; 
  //  R_out<<"molting probs female: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<<endl;
  R_out << "$Molting probability female" << endl;
  R_out<<moltp(1)<<endl;
  //  R_out<<"molting probs male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Molting probability male" << endl;
  R_out<<moltp(2)<<endl;
  //  R_out <<"Molting probability mature males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Molting probability mature males" << endl;
  R_out <<moltp_mat(2)<<endl;
  //  R_out << "observed pot fishery cpue 1979 fishery to endyr fishery: seq(1979,"<<endyr<<")" << endl;
  R_out << "$observed pot fishery cpue" << endl;
  R_out <<cpue(styr,endyr)<<endl;
  //  R_out << "predicted pot fishery cpue 1978 to endyr-1 survey: seq(1978,"<<endyr-1<<")" << endl;
  R_out << "$predicted pot fishery cpue" << endl;
  R_out <<cpue_pred(styr,endyr-1)<<endl;
  //  R_out << "observed retained catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$Observed retained catch biomass" << endl;
  R_out << catch_ret(styr,endyr-1) << endl;
  //  R_out << "predicted retained catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$Predicted retained catch biomas" << endl;
  R_out << pred_catch_ret(styr,endyr-1)<<endl;
  //  R_out << "predicted retained new catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained new catch biomass" << endl;
  R_out << (catch_male_ret_new*wtm)(styr,endyr-1)<<endl;
  //  R_out << "predicted retained old catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained old catch biomass" << endl;
  R_out << (catch_male_ret_old*wtm)(styr,endyr-1)<<endl;
  //  R_out << "observed retained+discard male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$observed retained plus discard male catch biomass" << endl;
  R_out << obs_catchtot_biom(styr,endyr-1) << endl;
  //  R_out << "predicted retained+discard male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained plus discard male catch biomass" << endl;
  R_out << pred_catch(styr,endyr-1) << endl;
  //  R_out << "predicted retained+discard new male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained plus discard new male catch biomass" << endl;
  R_out << (catch_lmale_new*wtm)(styr,endyr-1) << endl;
  //  R_out << "predicted retained+discard old male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained plus discard old male catch biomass" << endl;
  R_out << (catch_lmale_old*wtm)(styr,endyr-1) << endl;
  //  R_out << "observed discard male mortality biomass: seq(1979,"<<endyr<<")"<<endl;
  // R_out << "$FD" << endl;
  // R_out << (obs_catchtot_biom-catch_ret)(styr,endyr-1) <<endl;
  //  R_out << "predicted discard male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted discard male catch biomass" << endl;
  R_out << pred_catch(styr,endyr-1) -pred_catch_ret(styr,endyr-1)<< endl;
  //  R_out << "observed female discard mortality biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$observed female discard mortality biomass" << endl;
  R_out << obs_catchdf_biom(styr,endyr-1) << endl;
  //  R_out << "predicted female discard mortality biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted female discard mortality biomass" << endl;
  R_out << pred_catch_disc(1)(styr,endyr-1) << endl;
  //  R_out << "observed male discard mortality biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$observed male discard mortality biomass" << endl;
  R_out << obs_catchdm_biom(styr,endyr-1) << endl;
  //  R_out << "observed trawl catch biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$observed trawl catch biomass" << endl;
  R_out << obs_catcht_biom<<endl;
  //  R_out << "predicted trawl catch biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$predicted trawl catch biomass" << endl;
  R_out <<pred_catch_trawl<<endl;
  //  R_out << "estimated retained catch div. by male spawning biomass at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated retained catch div by male spawning biomass at fishtime" << endl;
  R_out <<elem_div(pred_catch_ret,mspbio_fishtime)(styr,endyr-1) << endl;
  //  R_out << "estimated total catch div. by male spawning biomass at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch div by male spawning biomass at fishtime" << endl;
  R_out <<elem_div(pred_catch,mspbio_fishtime)(styr,endyr-1) << endl;
  //  R_out << "estimated total catch of males >101 div. by males >101 at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch of males 101 div by males 101 at fishtime" << endl;
  R_out <<elem_div(pred_catch_gt101(styr,endyr-1),bio_males_gt101(styr,endyr-1)) << endl;
  //  R_out << "estimated total catch numbers of males >101 div. by males numbers >101 at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch numbers of males 101 div by males numbers 101 at fishtime" << endl;
  R_out <<elem_div(pred_catch_no_gt101(styr,endyr-1),num_males_gt101(styr,endyr-1)) << endl;
  //  R_out << "estimated total catch numbers of males >101 div. by survey estimate males numbers >101 at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch numbers of males 101 div by survey estimate males numbers 101 at fishtime" << endl;
     for(i=styr;i<endyr;i++)
        {
         obs_tmp(i) = obs_lmales(i-(styr-1));
        }
  R_out <<elem_div(pred_catch_no_gt101(styr,endyr-1),obs_tmp(styr,endyr-1)*mfexp(-M_matn(2)*(7/12))) << endl;
     for(i=styr;i<endyr;i++)
        {
         obs_tmp(i) = obs_lmales_bio(i-(styr-1));
        }
  //  R_out << "estimated total catch biomass of males >101 div. by survey estimate male biomass >101 at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch biomass of males 101 div by survey estimate male biomass 101 at fishtime" << endl;
  R_out <<elem_div(pred_catch_gt101(styr,endyr-1),obs_tmp(styr,endyr-1)*mfexp(-M_matn(2)*(7/12)) ) << endl;
  //  R_out << "estimated total catch biomass div. by survey estimate male mature biomass at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch biomass div. by survey estimate male mature biomass at fishtime" << endl;
  R_out <<elem_div(pred_catch(styr,endyr-1),((obs_srv1_spbiom(2))(styr,endyr-1))*mfexp(-M_matn(2)*(7/12))) << endl;
  //  R_out << "estimated annual total fishing mortality: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated annual total fishing mortality" << endl;
  R_out << mfexp(log_avg_fmort+fmort_dev)(styr,endyr-1) << endl;
  //  R_out <<"retained F: seq(1978,"<<endyr<<")" << endl;
  R_out << "$retained F" << endl;
         for(i=styr;i<=endyr;i++){
          R_out <<F_ret(1,i)(22)<<" ";
         }
  R_out<<endl;
  //  R_out <<"ghl: seq(1978,"<<endyr<<")" << endl;
  R_out << "$ghl" << endl;
  R_out <<catch_ghl/2.2<<endl;
    //  R_out << "estimated annual fishing mortality females pot: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated annual fishing mortality females pot" << endl;
  R_out << fmortdf(styr,endyr-1) <<endl;
  //  R_out << "estimated annual fishing mortality trawl bycatch: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated annual fishing mortality trawl bycatch" << endl;
  R_out << fmortt(styr,endyr-1) <<endl;
  R_out << "$estimated number of recruits female" << endl;
  for(i=styr; i<=endyr-1; i++)
  {
    R_out << mfexp(mean_log_rec(1)+rec_dev(1,i))<<" ";
  }
  R_out<<endl;
  //  R_out <<endl<< "estimated number of recruitments male: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated number of recruits male" << endl;
  for(i=styr; i<=endyr-1; i++)
  {
    R_out << mfexp(mean_log_rec(2)+rec_dev(2,i))<<" ";
  }
  for(i=1;i<=median_rec_yrs;i++)R_out<<2*median_rec<<" ";
    R_out<<endl;
 // R_out<<"distribution of recruits to length bins: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<<endl;
  R_out << "$distribution of recruits to length bins" << endl;
  R_out<<rec_len<<endl;
 // R_out<<"fishery total selectivity new shell 50% parameter: seq(1979,"<<endyr<<")"<<endl;
  R_out << "$fishery total selectivity new shell 50 parameter" << endl;
  R_out <<mfexp(log_avg_sel50_mn+log_sel50_dev_mn)(styr,endyr-1)<<endl;
 // R_out <<"fishery total selectivity old shell 50% parameter: seq(1979,"<<endyr<<")"<<endl;
  R_out << "$fishery total selectivity old shell 50 parameter" << endl;
  R_out <<mfexp(log_avg_sel50_mo+log_sel50_dev_mo)(styr,endyr-1)<<endl;
  R_out << "$selectivity fishery total new male" << endl;
  R_out << sel(1) << endl;
 // R_out << "selectivity fishery total old males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity fishery total old male" << endl;
  R_out << sel(2) << endl;
 // R_out << "selectivity fishery ret new males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity fishery retained new male" << endl;
  R_out << sel_fit(1) << endl;
 // R_out << "selectivity fishery ret old males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity fishery retained old male" << endl;
  R_out << sel_fit(2) << endl;
 // R_out <<"retention curve males new: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$retention curve new male" << endl;
  R_out <<sel_ret(1,endyr-1)<<endl;
 // R_out <<"retention curve males old: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$retention curve old male" << endl;
  R_out <<sel_ret(2,endyr-1)<<endl;
 // R_out << "selectivity discard females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity discard female" << endl;
  R_out <<sel_discf<<endl;
 // R_out << "selectivity trawl females:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity trawl female" << endl;
  R_out <<sel_trawl(1)<<endl;
 // R_out << "selectivity trawl males:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity trawl male" << endl;
  R_out <<sel_trawl(2)<<endl;
 // R_out << "selectivity survey females 1978 1981: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey female Era 1" << endl;
  R_out << sel_srv1(1) << endl;
 // R_out << "selectivity survey males 1978 1981: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey male Era 1" << endl;
  R_out << sel_srv1(2) << endl;
 // R_out << "selectivity survey females 1982 to 1988: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey female Era 2" << endl;
  R_out << sel_srv2(1) << endl;
 // R_out << "selectivity survey males 1982 to 1988: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey male Era 2" << endl;
  R_out << sel_srv2(2) << endl;
 // R_out << "selectivity survey females 1989 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey female Era 3" << endl;
  R_out << sel_srv3(1) << endl;
 // R_out << "selectivity survey males 1989 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey male Era 3" << endl;
  R_out << sel_srv3(2) << endl;
 // R_out << "selectivity industry survey females 2009: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity industry survey females 2009" << endl;
  R_out << sel_srvind(1) << endl;
 // R_out << "selectivity industry survey males 2009: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity industry survey males 2009" << endl;
  R_out << sel_srvind(2) << endl;
 // R_out << "selectivity nmfs industry survey females 2009: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity nmfs industry survey females 2009" << endl;
  R_out << sel_srvnmfs(1) << endl;
 // R_out << "selectivity nmfs industry survey males 2009: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity nmfs industry survey males 2009" << endl;
  R_out << sel_srvnmfs(2) << endl;
 // R_out << "selectivity industry survey females 2010: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity industry survey females 2010" << endl;
  R_out << sel_srv10ind(1) << endl;
 // R_out << "selectivity industry survey males 2010: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity industry survey males 2010" << endl;
  R_out << sel_srv10ind(2) << endl;
 // R_out << "selectivity nmfs industry survey females 2010: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity nmfs industry survey females 2010" << endl;
  R_out << sel_srv10nmfs(1) << endl;
 // R_out << "selectivity nmfs industry survey males 2010: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity nmfs industry survey males 2010" << endl;
  R_out << sel_srv10nmfs(2) << endl;
  R_out << "$survey CV" << endl;
  R_out << cv_srv1o << endl;
  R_out << "$Observed Lbar" <<endl;
  R_out << Lbar <<endl;
  R_out << "$Predicted Lbar new shell" <<endl;
  R_out << Lbar_hat_new <<endl;
  R_out << "$Predicted Lbar old shell" <<endl;
  R_out << Lbar_hat_old <<endl;
 }
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{500,1000,1000,1000,1000,1000,3000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1,1,1,1,.01,.001,1e-3,1e-3}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::final_calcs()
{
 time(&finish); 
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_post;
  pad_post = NULL;
}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 3500000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  // the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  time(&start);
  CheckFile.open("Check.Out");
  R_out.open("R_input.txt");
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
