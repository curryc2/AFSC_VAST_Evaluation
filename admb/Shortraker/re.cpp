  #include <admodel.h>
  #undef REPORT
  #define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("rwout.rep");
  adstring sppname;
#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <re.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_evalout = new ofstream("evalout.prj");;
  styr.allocate("styr");
  endyr.allocate("endyr");
  yrs.allocate(styr,endyr);
 yrs.fill_seqadd(styr,1);
  num_indx.allocate("num_indx");
  n_PE.allocate("n_PE");
  PE_vec.allocate(1,num_indx,"PE_vec");
  nobs.allocate("nobs");
  yrs_srv.allocate(1,nobs,"yrs_srv");
  srv_est.allocate(1,nobs,1,num_indx,"srv_est");
  srv_cv.allocate(1,nobs,1,num_indx,"srv_cv");
  srv_sd.allocate(1,nobs,1,num_indx);
 if (mean(srv_cv)>5) srv_cv = elem_div(srv_cv,srv_est+0.0001);
 srv_sd = elem_prod(srv_cv,srv_cv)+1;
 srv_sd = sqrt(log(srv_sd));
  yvar.allocate(1,nobs,1,num_indx);
  yconst.allocate(1,nobs,1,num_indx);
 yvar = elem_prod(srv_sd,srv_sd);
 yconst = log(2.0*M_PI*yvar);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  logSdLam.allocate(1,n_PE,-5,2,"logSdLam");
  biomsd.allocate(styr,endyr,1,num_indx,"biomsd");
  biomA.allocate(styr,endyr,1,num_indx,"biomA");
  biom.allocate(styr,endyr,1,num_indx,"biom");
  Like.allocate("Like");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  jnll =0.0;
  ofstream& evalout= *pad_evalout;
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
}

void SEPFUN1  model_parameters::step(const dvariable& biom1, const dvariable& biom2, const dvariable& logSdLam)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  dvariable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);
  end_df1b2_funnel();
}

void SEPFUN1  model_parameters::obs(const dvariable& biom, int i, int j)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  jnll+=0.5*(yconst(i,j) + square(biom-log(srv_est(i,j)+0.0001))/yvar(i,j));
  end_df1b2_funnel();
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
  biomsd = biom;
  Like = jnll;
  report << srv_sd <<endl;
}

void model_parameters::final_calcs()
{
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
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
  ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_MAX_NVAR_OFFSET(777000);
  arrmblsize = 777000;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
      if (!arrmblsize) arrmblsize=150000;
    df1b2variable::noallocate=1;
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    initial_df1b2params::separable_flag=1;
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

void model_parameters::preliminary_calculations(void){
  #if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

  #endif

}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

void df1b2_parameters::user_function(void)
{
  jnll =0.0;
  ofstream& evalout= *pad_evalout;
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
}

void   df1b2_pre_parameters::step(const funnel_init_df1b2variable& biom1, const funnel_init_df1b2variable& biom2, const funnel_init_df1b2variable& logSdLam)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  df1b2variable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);
  end_df1b2_funnel();
}

void   df1b2_pre_parameters::obs(const funnel_init_df1b2variable& biom, int i, int j)
{
  begin_df1b2_funnel();
  ofstream& evalout= *pad_evalout;
  jnll+=0.5*(yconst(i,j) + square(biom-log(srv_est(i,j)+0.0001))/yvar(i,j));
  end_df1b2_funnel();
}
   
void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
  df1b2_gradlist::set_no_derivatives(); 
  quadratic_prior::in_qp_calculations=1; 
}  
  
void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
  (*re_objective_function_value::pobjfun)=0; 
  other_separable_stuff_begin(); 
  f1b2gradlist->reset();  
  if (!quadratic_prior::in_qp_calculations) 
  { 
    df1b2_gradlist::set_yes_derivatives();  
  } 
  funnel_init_var::allocate_all();  
}  
 
void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
  lapprox->do_separable_stuff(); 
  other_separable_stuff_end(); 
} 
  
void model_parameters::begin_df1b2_funnel(void) 
{ 
  if (lapprox)  
  {  
    {  
      begin_funnel_stuff();  
    }  
  }  
}  
 
void model_parameters::end_df1b2_funnel(void) 
{  
  if (lapprox)  
  {  
    end_df1b2_funnel_stuff();  
  }  
} 

void df1b2_parameters::allocate(void) 
{
  logSdLam.allocate(1,n_PE,-5,2,"logSdLam");
  biomsd.allocate(styr,endyr,1,num_indx,"biomsd");
  biomA.allocate(styr,endyr,1,num_indx,"biomA");
  biom.allocate(styr,endyr,1,num_indx,"biom");
  Like.allocate("Like");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  jnll.allocate("jnll");  /* ADOBJECTIVEFUNCTION */
}
