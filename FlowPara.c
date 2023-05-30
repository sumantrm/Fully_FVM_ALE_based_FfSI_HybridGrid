#include "DesignAndFlowData.h"

void define_flow_para()
{

  delt	= Dt; // time-step-size
  rho	= 1.0; // density
  mu	= 0.01; // viscosity
  count	= 0;

  u_inlet_max = 1.0;

  wu_vel= 1.0;
  wv_vel= 1.0;
  w_p	= 1.0;
  Sm_st_new = 0.0;

  w_s_u	= 1.0;
  w_s_v	= 1.0;
  
  L_zeta= L;
  L_eta	= D;
  it	= 1;
  
  per 	= 0.3; // percentage of SEMM domain during mesh generation

}
