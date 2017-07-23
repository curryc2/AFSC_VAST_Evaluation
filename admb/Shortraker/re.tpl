 // Random walk model for survey averaging using log-scale process and observation errors
DATA_SECTION

  !!CLASS ofstream evalout("evalout.prj");

  init_int styr
  init_int endyr
  ivector yrs(styr,endyr);
  !! yrs.fill_seqadd(styr,1);
  init_int num_indx;
  init_int n_PE
  init_ivector PE_vec(1,num_indx)
  init_int nobs;
  init_ivector yrs_srv(1,nobs);
  init_matrix srv_est(1,nobs,1,num_indx)
  init_matrix srv_cv(1,nobs,1,num_indx)
  matrix srv_sd(1,nobs,1,num_indx)
  !! if (mean(srv_cv)>5) srv_cv = elem_div(srv_cv,srv_est+0.0001);
  !! srv_sd = elem_prod(srv_cv,srv_cv)+1;
  !! srv_sd = sqrt(log(srv_sd));
  matrix yvar(1,nobs,1,num_indx)
  matrix yconst(1,nobs,1,num_indx)
  !! yvar = elem_prod(srv_sd,srv_sd);
  !! yconst = log(2.0*M_PI*yvar);

PARAMETER_SECTION
  init_bounded_vector logSdLam(1,n_PE,-5,2)
  sdreport_matrix biomsd(styr,endyr,1,num_indx);
  sdreport_matrix biomA(styr,endyr,1,num_indx);
  random_effects_matrix biom(styr,endyr,1,num_indx);
  sdreport_number Like;
  objective_function_value jnll;

PROCEDURE_SECTION
  jnll=0.0;

  for(int j=1; j<=num_indx; ++j)
  {
  for(int i=styr+1; i<=endyr; ++i)
  {
    step(biom(i-1,j),biom(i,j),logSdLam(PE_vec(j)));
  }
  for(int i=1; i<=nobs; ++i)
  {
    if(srv_est(i,j)>-1) obs(biom(yrs_srv(i),j),i,j);
  }
  }

  if (sd_phase()) 
  {
    biomA = exp(biom);
    biomsd = biom;
  }

  if (mceval_phase()){
   evalout<<logSdLam<<" "<<jnll<<endl;}


SEPARABLE_FUNCTION void step(const dvariable& biom1, const dvariable& biom2, const dvariable& logSdLam)
  dvariable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);

SEPARABLE_FUNCTION void obs(const dvariable& biom, int i, int j)
  jnll+=0.5*(yconst(i,j) + square(biom-log(srv_est(i,j)+0.0001))/yvar(i,j));

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(777000);
  arrmblsize = 777000;

REPORT_SECTION
  biomsd = biom;
  Like = jnll;
  report << srv_sd <<endl;

GLOBALS_SECTION
  #include <admodel.h>
  #undef REPORT
  #define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("rwout.rep");
  adstring sppname;

FINAL_SECTION
  dvar_vector srv_est_TOT = rowsum(srv_est);
  dvar_vector biom_TOT = rowsum(biomA);
  dvar_vector SD_numer = rowsum(elem_prod(exp(2*biomsd+square(biomsd.sd)),(exp(square(biomsd.sd))-1)));
  dvar_vector SD_denom = square(rowsum(exp(biomsd+0.5*square(biomsd.sd))));
  dvar_vector SD_biom_TOT = sqrt(log(elem_div(SD_numer,SD_denom)+1));
  dvar_vector biom_TOT_UCI = exp(log(biom_TOT)+1.96*SD_biom_TOT);
  dvar_vector biom_TOT_LCI = exp(log(biom_TOT)-1.96*SD_biom_TOT);
  dvar_matrix UCI = exp(biomsd+1.96*biomsd.sd);
  dvar_matrix LCI = exp(biomsd-1.96*biomsd.sd);

  write_R(yrs_srv);
  write_R(srv_est_TOT);
  write_R(yrs);
  write_R(biom_TOT);
  write_R(SD_biom_TOT);
  write_R(biom_TOT_UCI);
  write_R(biom_TOT_LCI);
  write_R(yrs_srv);
  write_R(srv_est);
  write_R(srv_sd);
  write_R(yrs);
  write_R(LCI);
  write_R(biomA);
  write_R(UCI);
  write_R(biomsd);
  write_R(biomsd.sd);

  mysum.close();
