#include "DesignAndFlowData.h"

void solve_fluid()
{

  /* calculatjon of u* value */
  for(j=1;j<=Nx-2;j++)
  for(i=1;i<=Ny-2;i++)
   {
     p_ne[i][j]=(del_v[i+1][j+1]*p[i][j]+del_v[i][j]*p[i+1][j+1]+del_v[i+1][j]*p[i][j+1]+del_v[i][j+1]*p[i+1][j])/(del_v[i+1][j+1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j+1]);
     p_se[i][j]=(del_v[i-1][j+1]*p[i][j]+del_v[i][j]*p[i-1][j+1]+del_v[i-1][j]*p[i][j+1]+del_v[i][j+1]*p[i-1][j])/(del_v[i-1][j+1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j+1]);
     p_nw[i][j]=(del_v[i+1][j-1]*p[i][j]+del_v[i][j]*p[i+1][j-1]+del_v[i+1][j]*p[i][j-1]+del_v[i][j-1]*p[i+1][j])/(del_v[i+1][j-1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j-1]);
     p_sw[i][j]=(del_v[i-1][j-1]*p[i][j]+del_v[i][j]*p[i-1][j-1]+del_v[i-1][j]*p[i][j-1]+del_v[i][j-1]*p[i-1][j])/(del_v[i-1][j-1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j-1]);
   }

maxerror_u=1.0;
while(maxerror_u>1e-6)
 {
   maxerror_u=0.0;
   for(j=1;j<=Nx-2;j++)
   for(i=Ny-2;i>=1;i--)
    {

     u_ne[i][j]=(del_v[i+1][j+1]*u_b[i][j]+del_v[i][j]*u_b[i+1][j+1]+del_v[i+1][j]*u_b[i][j+1]+del_v[i][j+1]*u_b[i+1][j])/(del_v[i+1][j+1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j+1]);
     u_se[i][j]=(del_v[i-1][j+1]*u_b[i][j]+del_v[i][j]*u_b[i-1][j+1]+del_v[i-1][j]*u_b[i][j+1]+del_v[i][j+1]*u_b[i-1][j])/(del_v[i-1][j+1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j+1]);
     u_nw[i][j]=(del_v[i+1][j-1]*u_b[i][j]+del_v[i][j]*u_b[i+1][j-1]+del_v[i+1][j]*u_b[i][j-1]+del_v[i][j-1]*u_b[i+1][j])/(del_v[i+1][j-1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j-1]);
     u_sw[i][j]=(del_v[i-1][j-1]*u_b[i][j]+del_v[i][j]*u_b[i-1][j-1]+del_v[i-1][j]*u_b[i][j-1]+del_v[i][j-1]*u_b[i-1][j])/(del_v[i-1][j-1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j-1]);

     /* defining mass terms */
     m_e=rho*((del_v[i][j]*(u_b[i][j+1]-umesh[i][j+1])+del_v[i][j+1]*(u_b[i][j]-umesh[i][j]))*dels_xe[i][j]/(del_v[i][j]+del_v[i][j+1]) + (del_v[i][j]*(v_b[i][j+1]-vmesh[i][j+1])+del_v[i][j+1]*(v_b[i][j]-vmesh[i][j]))*dels_ye[i][j]/(del_v[i][j]+del_v[i][j+1]))/dels_e[i][j];
     m_w=rho*((del_v[i][j]*(u_b[i][j-1]-umesh[i][j-1])+del_v[i][j-1]*(u_b[i][j]-umesh[i][j]))*dels_xw[i][j]/(del_v[i][j]+del_v[i][j-1]) + (del_v[i][j]*(v_b[i][j-1]-vmesh[i][j-1])+del_v[i][j-1]*(v_b[i][j]-vmesh[i][j]))*dels_yw[i][j]/(del_v[i][j]+del_v[i][j-1]))/dels_w[i][j];
     m_n=rho*((del_v[i][j]*(u_b[i+1][j]-umesh[i+1][j])+del_v[i+1][j]*(u_b[i][j]-umesh[i][j]))*dels_xn[i][j]/(del_v[i][j]+del_v[i+1][j]) + (del_v[i][j]*(v_b[i+1][j]-vmesh[i+1][j])+del_v[i+1][j]*(v_b[i][j]-vmesh[i][j]))*dels_yn[i][j]/(del_v[i][j]+del_v[i+1][j]))/dels_n[i][j];
     m_s=rho*((del_v[i][j]*(u_b[i-1][j]-umesh[i-1][j])+del_v[i-1][j]*(u_b[i][j]-umesh[i][j]))*dels_xs[i][j]/(del_v[i][j]+del_v[i-1][j]) + (del_v[i][j]*(v_b[i-1][j]-vmesh[i-1][j])+del_v[i-1][j]*(v_b[i][j]-vmesh[i][j]))*dels_ys[i][j]/(del_v[i][j]+del_v[i-1][j]))/dels_s[i][j];

     we1_p=del_v[i][j]/(del_v[i][j+1]+del_v[i][j]);
     we1_m=del_v[i][j+1]/(del_v[i][j]+del_v[i][j+1]);
     we2_p=del_v[i][j+1]/(del_v[i][j]+del_v[i][j+1]);
     we2_m=del_v[i][j]/(del_v[i][j]+del_v[i][j+1]);
     we3_p=0.0;
     we3_m=0.0;

     ww1_p=del_v[i][j-1]/(del_v[i][j-1]+del_v[i][j]);
     ww1_m=del_v[i][j]/(del_v[i][j]+del_v[i][j-1]);
     ww2_p=del_v[i][j]/(del_v[i][j]+del_v[i][j-1]);
     ww2_m=del_v[i][j-1]/(del_v[i][j]+del_v[i][j-1]);
     ww3_p=0.0;
     ww3_m=0.0;

     wn1_p=del_v[i][j]/(del_v[i+1][j]+del_v[i][j]);
     wn1_m=del_v[i+1][j]/(del_v[i][j]+del_v[i+1][j]);
     wn2_p=del_v[i+1][j]/(del_v[i][j]+del_v[i+1][j]);
     wn2_m=del_v[i][j]/(del_v[i][j]+del_v[i+1][j]);
     wn3_p=0.0;
     wn3_m=0.0;

     ws1_p=del_v[i-1][j]/(del_v[i-1][j]+del_v[i][j]);
     ws1_m=del_v[i][j]/(del_v[i][j]+del_v[i-1][j]);
     ws2_p=del_v[i][j]/(del_v[i][j]+del_v[i-1][j]);
     ws2_m=del_v[i-1][j]/(del_v[i][j]+del_v[i-1][j]);
     ws3_p=0.0;
     ws3_m=0.0;

     ue_p=we1_p*u_b[i][j+1]+we2_p*u_b[i][j]+we3_p*u_b[i][j-1];
     ue_m=we1_m*u_b[i][j]+we2_m*u_b[i][j+1]+we3_m*u_b[j+2][i];

     uw_p=ww1_p*u_b[i][j]+ww2_p*u_b[i][j-1]+ww3_p*u_b[j-2][i];
     uw_m=ww1_m*u_b[i][j-1]+ww2_m*u_b[i][j]+ww3_m*u_b[i][j+1];

     un_p=wn1_p*u_b[i+1][j]+wn2_p*u_b[i][j]+wn3_p*u_b[i-1][j];
     un_m=wn1_m*u_b[i][j]+wn2_m*u_b[i+1][j]+wn3_m*u_b[j][i+2];

     us_p=ws1_p*u_b[i][j]+ws2_p*u_b[i-1][j]+ws3_p*u_b[j][i-2];
     us_m=ws1_m*u_b[i-1][j]+ws2_m*u_b[i][j]+ws3_m*u_b[i+1][j];

     /*  defining advection term */
     au=((m_e > 0 ? m_e : 0)*ue_p + (m_e < 0 ? m_e : 0)*ue_m)*dels_e[i][j] - ((m_w > 0 ? m_w : 0)*uw_p + (m_w < 0 ? m_w : 0)*uw_m)*dels_w[i][j] + ((m_n > 0 ? m_n : 0)*un_p +(m_n < 0 ? m_n : 0)*un_m)*dels_n[i][j] - ((m_s > 0 ? m_s : 0)*us_p +(m_s < 0 ? m_s : 0)*us_m)*dels_s[i][j] ;

     ap=rho*del_v[i][j]/delt + mu*(dels_se[i][j]/ds_e[i][j]+dels_sw[i][j]/ds_w[i][j]+dels_sn[i][j]/ds_n[i][j]+dels_ss[i][j]/ds_s[i][j]);

     ub=(1.0/ap)*( rho*del_v[i][j]*u[i][j]/delt - au + mu*(u_b[i][j+1]*dels_se[i][j]/ds_e[i][j]+u_b[i][j-1]*dels_sw[i][j]/ds_w[i][j]+u_b[i+1][j]*dels_sn[i][j]/ds_n[i][j]+u_b[i-1][j]*dels_ss[i][j]/ds_s[i][j]+(u_ne[i][j]-u_se[i][j])*dels_te[i][j]/dt_e[i][j]-(u_nw[i][j]-u_sw[i][j])*dels_tw[i][j]/dt_w[i][j]+(u_ne[i][j]-u_nw[i][j])*dels_tn[i][j]/dt_n[i][j]-(u_se[i][j]-u_sw[i][j])*dels_ts[i][j]/dt_s[i][j]) );

     if(fabs(ub-u_b[i][j])>maxerror_u)
       maxerror_u=fabs(ub-u_b[i][j]);

     u_b[i][j]=ub; 
   
    }
}

maxerror_v=1.0;
while(maxerror_v>1e-6)
 {
  maxerror_v=0.0;
  for(j=1;j<=Nx-2;j++)
  for(i=Ny-2;i>=1;i--)
   {

     v_ne[i][j]=(del_v[i+1][j+1]*v_b[i][j]+del_v[i][j]*v_b[i+1][j+1]+del_v[i+1][j]*v_b[i][j+1]+del_v[i][j+1]*v_b[i+1][j])/(del_v[i+1][j+1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j+1]);
     v_se[i][j]=(del_v[i-1][j+1]*v_b[i][j]+del_v[i][j]*v_b[i-1][j+1]+del_v[i-1][j]*v_b[i][j+1]+del_v[i][j+1]*v_b[i-1][j])/(del_v[i-1][j+1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j+1]);
     v_nw[i][j]=(del_v[i+1][j-1]*v_b[i][j]+del_v[i][j]*v_b[i+1][j-1]+del_v[i+1][j]*v_b[i][j-1]+del_v[i][j-1]*v_b[i+1][j])/(del_v[i+1][j-1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j-1]);
     v_sw[i][j]=(del_v[i-1][j-1]*v_b[i][j]+del_v[i][j]*v_b[i-1][j-1]+del_v[i-1][j]*v_b[i][j-1]+del_v[i][j-1]*v_b[i-1][j])/(del_v[i-1][j-1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j-1]);

     /* defining mass terms */    
     m_e=rho*((del_v[i][j]*(u_b[i][j+1]-umesh[i][j+1])+del_v[i][j+1]*(u_b[i][j]-umesh[i][j]))*dels_xe[i][j]/(del_v[i][j]+del_v[i][j+1]) + (del_v[i][j]*(v_b[i][j+1]-vmesh[i][j+1])+del_v[i][j+1]*(v_b[i][j]-vmesh[i][j]))*dels_ye[i][j]/(del_v[i][j]+del_v[i][j+1]))/dels_e[i][j];
     m_w=rho*((del_v[i][j]*(u_b[i][j-1]-umesh[i][j-1])+del_v[i][j-1]*(u_b[i][j]-umesh[i][j]))*dels_xw[i][j]/(del_v[i][j]+del_v[i][j-1]) + (del_v[i][j]*(v_b[i][j-1]-vmesh[i][j-1])+del_v[i][j-1]*(v_b[i][j]-vmesh[i][j]))*dels_yw[i][j]/(del_v[i][j]+del_v[i][j-1]))/dels_w[i][j];
     m_n=rho*((del_v[i][j]*(u_b[i+1][j]-umesh[i+1][j])+del_v[i+1][j]*(u_b[i][j]-umesh[i][j]))*dels_xn[i][j]/(del_v[i][j]+del_v[i+1][j]) + (del_v[i][j]*(v_b[i+1][j]-vmesh[i+1][j])+del_v[i+1][j]*(v_b[i][j]-vmesh[i][j]))*dels_yn[i][j]/(del_v[i][j]+del_v[i+1][j]))/dels_n[i][j];
     m_s=rho*((del_v[i][j]*(u_b[i-1][j]-umesh[i-1][j])+del_v[i-1][j]*(u_b[i][j]-umesh[i][j]))*dels_xs[i][j]/(del_v[i][j]+del_v[i-1][j]) + (del_v[i][j]*(v_b[i-1][j]-vmesh[i-1][j])+del_v[i-1][j]*(v_b[i][j]-vmesh[i][j]))*dels_ys[i][j]/(del_v[i][j]+del_v[i-1][j]))/dels_s[i][j];

     we1_p=del_v[i][j]/(del_v[i][j+1]+del_v[i][j]);
     we1_m=del_v[i][j+1]/(del_v[i][j]+del_v[i][j+1]);
     we2_p=del_v[i][j+1]/(del_v[i][j]+del_v[i][j+1]);
     we2_m=del_v[i][j]/(del_v[i][j]+del_v[i][j+1]);
     we3_p=0.0;
     we3_m=0.0;

     ww1_p=del_v[i][j-1]/(del_v[i][j-1]+del_v[i][j]);
     ww1_m=del_v[i][j]/(del_v[i][j]+del_v[i][j-1]);
     ww2_p=del_v[i][j]/(del_v[i][j]+del_v[i][j-1]);
     ww2_m=del_v[i][j-1]/(del_v[i][j]+del_v[i][j-1]);
     ww3_p=0.0;
     ww3_m=0.0;

     wn1_p=del_v[i][j]/(del_v[i+1][j]+del_v[i][j]);
     wn1_m=del_v[i+1][j]/(del_v[i][j]+del_v[i+1][j]);
     wn2_p=del_v[i+1][j]/(del_v[i][j]+del_v[i+1][j]);
     wn2_m=del_v[i][j]/(del_v[i][j]+del_v[i+1][j]);
     wn3_p=0.0;
     wn3_m=0.0;

     ws1_p=del_v[i-1][j]/(del_v[i-1][j]+del_v[i][j]);
     ws1_m=del_v[i][j]/(del_v[i][j]+del_v[i-1][j]);
     ws2_p=del_v[i][j]/(del_v[i][j]+del_v[i-1][j]);
     ws2_m=del_v[i-1][j]/(del_v[i][j]+del_v[i-1][j]);
     ws3_p=0.0;
     ws3_m=0.0;

     ve_p=we1_p*v_b[i][j+1]+we2_p*v_b[i][j]+we3_p*v_b[i][j-1];
     ve_m=we1_m*v_b[i][j]+we2_m*v_b[i][j+1]+we3_m*v_b[j+2][i];

     vw_p=ww1_p*v_b[i][j]+ww2_p*v_b[i][j-1]+ww3_p*v_b[j-2][i];
     vw_m=ww1_m*v_b[i][j-1]+ww2_m*v_b[i][j]+ww3_m*v_b[i][j+1];

     vn_p=wn1_p*v_b[i+1][j]+wn2_p*v_b[i][j]+wn3_p*v_b[i-1][j];
     vn_m=wn1_m*v_b[i][j]+wn2_m*v_b[i+1][j]+wn3_m*v_b[j][i+2];

     vs_p=ws1_p*v_b[i][j]+ws2_p*v_b[i-1][j]+ws3_p*v_b[j][i-2];
     vs_m=ws1_m*v_b[i-1][j]+ws2_m*v_b[i][j]+ws3_m*v_b[i+1][j];

     /*  defining advectjon term */
     av=((m_e > 0 ? m_e : 0)*ve_p + (m_e < 0 ? m_e : 0)*ve_m)*dels_e[i][j] - ((m_w > 0 ? m_w : 0)*vw_p + (m_w < 0 ? m_w : 0)*vw_m)*dels_w[i][j] + ((m_n > 0 ? m_n : 0)*vn_p +(m_n < 0 ? m_n : 0)*vn_m)*dels_n[i][j] - ((m_s > 0 ? m_s : 0)*vs_p +(m_s < 0 ? m_s : 0)*vs_m)*dels_s[i][j] ;       

     ap=rho*del_v[i][j]/delt + mu*(dels_se[i][j]/ds_e[i][j]+dels_sw[i][j]/ds_w[i][j]+dels_sn[i][j]/ds_n[i][j]+dels_ss[i][j]/ds_s[i][j]);

     vb=(1.0/ap)*( rho*del_v[i][j]*v[i][j]/delt - av + mu*(v_b[i][j+1]*dels_se[i][j]/ds_e[i][j]+v_b[i][j-1]*dels_sw[i][j]/ds_w[i][j]+v_b[i+1][j]*dels_sn[i][j]/ds_n[i][j]+v_b[i-1][j]*dels_ss[i][j]/ds_s[i][j]+(v_ne[i][j]-v_se[i][j])*dels_te[i][j]/dt_e[i][j]-(v_nw[i][j]-v_sw[i][j])*dels_tw[i][j]/dt_w[i][j]+(v_ne[i][j]-v_nw[i][j])*dels_tn[i][j]/dt_n[i][j]-(v_se[i][j]-v_sw[i][j])*dels_ts[i][j]/dt_s[i][j]) );

     if(fabs(vb-v_b[i][j])>maxerror_v)
	maxerror_v=fabs(vb-v_b[i][j]);

     v_b[i][j]=vb;
   }
}

  for(j=0;j<=Nx-1;j++)
  for(i=0;i<=Ny-1;i++)
   {
   
     p_new[i][j]=p[i][j];
     p_old[i][j]=p[i][j];
     
   }

  max_Sm_st=0.0;
  for(j=1;j<=Nx-2;j++)
  for(i=Ny-2;i>=1;i--)
    {

      ue=(del_v[i][j]*u_b[i][j+1]+del_v[i][j+1]*u_b[i][j])/(del_v[i][j]+del_v[i][j+1]);
      uw=(del_v[i][j]*u_b[i][j-1]+del_v[i][j-1]*u_b[i][j])/(del_v[i][j]+del_v[i][j-1]);
      un=(del_v[i][j]*u_b[i+1][j]+del_v[i+1][j]*u_b[i][j])/(del_v[i][j]+del_v[i+1][j]);
      us=(del_v[i][j]*u_b[i-1][j]+del_v[i-1][j]*u_b[i][j])/(del_v[i][j]+del_v[i-1][j]);

      ve=(del_v[i][j]*v_b[i][j+1]+del_v[i][j+1]*v_b[i][j])/(del_v[i][j]+del_v[i][j+1]);
      vw=(del_v[i][j]*v_b[i][j-1]+del_v[i][j-1]*v_b[i][j])/(del_v[i][j]+del_v[i][j-1]);
      vn=(del_v[i][j]*v_b[i+1][j]+del_v[i+1][j]*v_b[i][j])/(del_v[i][j]+del_v[i+1][j]);
      vs=(del_v[i][j]*v_b[i-1][j]+del_v[i-1][j]*v_b[i][j])/(del_v[i][j]+del_v[i-1][j]);

      m_st_e[i][j]=(rho*((ue)*dels_xe[i][j]+(ve)*dels_ye[i][j]))/dels_e[i][j];
      m_st_w[i][j]=(rho*((uw)*dels_xw[i][j]+(vw)*dels_yw[i][j]))/dels_w[i][j];
      m_st_n[i][j]=(rho*((un)*dels_xn[i][j]+(vn)*dels_yn[i][j]))/dels_n[i][j];
      m_st_s[i][j]=(rho*((us)*dels_xs[i][j]+(vs)*dels_ys[i][j]))/dels_s[i][j];

      Sm_st[i][j]=m_st_e[i][j]*dels_e[i][j]+m_st_n[i][j]*dels_n[i][j]-m_st_s[i][j]*dels_s[i][j]-m_st_w[i][j]*dels_w[i][j];

      if(fabs(Sm_st[i][j])>max_Sm_st)
        max_Sm_st=fabs(Sm_st[i][j]);

     p_new_ne[i][j]=(del_v[i+1][j+1]*p_new[i][j]+del_v[i][j]*p_new[i+1][j+1]+del_v[i+1][j]*p_new[i][j+1]+del_v[i][j+1]*p_new[i+1][j])/(del_v[i+1][j+1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j+1]);
     p_new_se[i][j]=(del_v[i-1][j+1]*p_new[i][j]+del_v[i][j]*p_new[i-1][j+1]+del_v[i-1][j]*p_new[i][j+1]+del_v[i][j+1]*p_new[i-1][j])/(del_v[i-1][j+1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j+1]);
     p_new_nw[i][j]=(del_v[i+1][j-1]*p_new[i][j]+del_v[i][j]*p_new[i+1][j-1]+del_v[i+1][j]*p_new[i][j-1]+del_v[i][j-1]*p_new[i+1][j])/(del_v[i+1][j-1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j-1]);
     p_new_sw[i][j]=(del_v[i-1][j-1]*p_new[i][j]+del_v[i][j]*p_new[i-1][j-1]+del_v[i-1][j]*p_new[i][j-1]+del_v[i][j-1]*p_new[i-1][j])/(del_v[i-1][j-1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j-1]);
    }

    
  maxerror_p=1.0;
  while(maxerror_p>=1e-4)
  {

  maxerror_p=0.0;
  for(j=1;j<=Nx-2;j++)
  for(i=1;i<=Ny-2;i++)
   {

     ap = -( dels_se[i][j]/ds_e[i][j]+dels_sw[i][j]/ds_w[i][j]+dels_sn[i][j]/ds_n[i][j]+dels_ss[i][j]/ds_s[i][j] );
     p_new[i][j]=Sm_st[i][j]/delt-( p_new[i][j+1]*dels_se[i][j]/ds_e[i][j]+p_new[i][j-1]*dels_sw[i][j]/ds_w[i][j]+p_new[i+1][j]*dels_sn[i][j]/ds_n[i][j]+p_new[i-1][j]*dels_ss[i][j]/ds_s[i][j] )-((p_new_ne[i][j]-p_new_se[i][j])*dels_te[i][j]/dt_e[i][j] - (p_new_nw[i][j]-p_new_sw[i][j])*dels_tw[i][j]/dt_w[i][j] + (p_new_ne[i][j]-p_new_nw[i][j])*dels_tn[i][j]/dt_n[i][j] - (p_new_se[i][j]-p_new_sw[i][j])*dels_ts[i][j]/dt_s[i][j] );
     p_new[i][j]=p_new[i][j]/ap;

     error=fabs(p_new[i][j]-p_old[i][j]);
     if(error>maxerror_p)
     maxerror_p=error;
      
     p_new_ne[i][j]=(del_v[i+1][j+1]*p_new[i][j]+del_v[i][j]*p_new[i+1][j+1]+del_v[i+1][j]*p_new[i][j+1]+del_v[i][j+1]*p_new[i+1][j])/(del_v[i+1][j+1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j+1]);
     p_new_se[i][j]=(del_v[i-1][j+1]*p_new[i][j]+del_v[i][j]*p_new[i-1][j+1]+del_v[i-1][j]*p_new[i][j+1]+del_v[i][j+1]*p_new[i-1][j])/(del_v[i-1][j+1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j+1]);
     p_new_nw[i][j]=(del_v[i+1][j-1]*p_new[i][j]+del_v[i][j]*p_new[i+1][j-1]+del_v[i+1][j]*p_new[i][j-1]+del_v[i][j-1]*p_new[i+1][j])/(del_v[i+1][j-1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j-1]);
     p_new_sw[i][j]=(del_v[i-1][j-1]*p_new[i][j]+del_v[i][j]*p_new[i-1][j-1]+del_v[i-1][j]*p_new[i][j-1]+del_v[i][j-1]*p_new[i-1][j])/(del_v[i-1][j-1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j-1]);

   }

    //   BC
    // B.C 
    for(i=0;i<=Ny-1;i++)
     {	   

       if(y[i][0]>=D_f)
        {
          p_new[i][0]=0.0;
          p_new[i][Nx-1]=0.0;

        }
       else
        {	   
          p_new[i][0]=p_new[i][1];
          p_new[i][Nx-1]=p_new[i][Nx-2];

        }
        
     }

    for(j=0;j<=Nx-1;j++)
     {

        p_new[0][j]=p_new[1][j];//-ds_s[1][j]*((v_s_new[Ny_s-1][j]-2.0*v_s[Ny_s-1][j]+v_s_old[Ny_s-1][j])*ny +(u_s_new[Ny_s-1][j]-2.0*u_s[Ny_s-1][j]+u_s_old[Ny_s-1][j])*nx)/pow(delt,2);
        p_new[Ny-1][j]=p_new[Ny-2][j];

     }


    
    for(j=0;j<=Nx-1;j++)
    for(i=0;i<=Ny-1;i++)
      p_old[i][j] = p_new[i][j];


   }

  
  /* while(max_Sm_st>1e-6)
   {

   max_Sm_st=0.0;
   for(j=1;j<=Nx-2;j++)
   for(i=Ny-2;i>=1;i--)
    {
     
     // Sm_dash[i][j]=-delt*( (p_old[i][j+1])*dels_se[i][j]/ds_e[i][j] - (-p_old[i][j-1])*dels_sw[i][j]/ds_w[i][j] + (p_old[i+1][j])*dels_sn[i][j]/ds_n[i][j] - (-p_old[i-1][j])*dels_ss[i][j]/ds_s[i][j] + (p_ne[i][j]-p_se[i][j])*dels_te[i][j]/dt_e[i][j] - (p_nw[i][j]-p_sw[i][j])*dels_tw[i][j]/dt_w[i][j] + (p_ne[i][j]-p_nw[i][j])*dels_tn[i][j]/dt_n[i][j] - (p_se[i][j]-p_sw[i][j])*dels_ts[i][j]/dt_s[i][j] );

     Sm_dash[i][j]=-delt*( (p_new[i][j+1]-p_new[i][j])*dels_se[i][j]/ds_e[i][j] - (p_new[i][j]-p_new[i][j-1])*dels_sw[i][j]/ds_w[i][j] + (p_new[i+1][j]-p_new[i][j])*dels_sn[i][j]/ds_n[i][j] - (p_new[i][j]-p_new[i-1][j])*dels_ss[i][j]/ds_s[i][j] + (p_new_ne[i][j]-p_new_se[i][j])*dels_te[i][j]/dt_e[i][j] - (p_new_nw[i][j]-p_new_sw[i][j])*dels_tw[i][j]/dt_w[i][j] + (p_new_ne[i][j]-p_new_nw[i][j])*dels_tn[i][j]/dt_n[i][j] - (p_new_se[i][j]-p_new_sw[i][j])*dels_ts[i][j]/dt_s[i][j] );
     ap=delt*( dels_se[i][j]/ds_e[i][j]+dels_sw[i][j]/ds_w[i][j]+dels_sn[i][j]/ds_n[i][j]+dels_ss[i][j]/ds_s[i][j] );
    
     // defining pressure correctjon p_dash term 
     p_new[i][j] =p_new[i][j]-(1.0/ap)*(Sm_st[i][j]+Sm_dash[i][j]);

     p_new_ne[i][j]=(del_v[i+1][j+1]*p_new[i][j]+del_v[i][j]*p_new[i+1][j+1]+del_v[i+1][j]*p_new[i][j+1]+del_v[i][j+1]*p_new[i+1][j])/(del_v[i+1][j+1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j+1]);
     p_new_se[i][j]=(del_v[i-1][j+1]*p_new[i][j]+del_v[i][j]*p_new[i-1][j+1]+del_v[i-1][j]*p_new[i][j+1]+del_v[i][j+1]*p_new[i-1][j])/(del_v[i-1][j+1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j+1]);
     p_new_nw[i][j]=(del_v[i+1][j-1]*p_new[i][j]+del_v[i][j]*p_new[i+1][j-1]+del_v[i+1][j]*p_new[i][j-1]+del_v[i][j-1]*p_new[i+1][j])/(del_v[i+1][j-1]+del_v[i][j]+del_v[i+1][j]+del_v[i][j-1]);
     p_new_sw[i][j]=(del_v[i-1][j-1]*p_new[i][j]+del_v[i][j]*p_new[i-1][j-1]+del_v[i-1][j]*p_new[i][j-1]+del_v[i][j-1]*p_new[i-1][j])/(del_v[i-1][j-1]+del_v[i][j]+del_v[i-1][j]+del_v[i][j-1]);

     m_st_e[i][j]=m_st_e[i][j]-delt*((p_new[i][j+1]-p_new[i][j])*dels_se[i][j]/ds_e[i][j] + (p_new_ne[i][j]-p_new_se[i][j])*dels_te[i][j]/dt_e[i][j])/dels_e[i][j];
     m_st_w[i][j]=m_st_w[i][j]-delt*((p_new[i][j]-p_new[i][j-1])*dels_sw[i][j]/ds_w[i][j] + (p_new_nw[i][j]-p_new_sw[i][j])*dels_tw[i][j]/dt_w[i][j])/dels_w[i][j];
     m_st_n[i][j]=m_st_n[i][j]-delt*((p_new[i+1][j]-p_new[i][j])*dels_sn[i][j]/ds_n[i][j] + (p_new_ne[i][j]-p_new_nw[i][j])*dels_tn[i][j]/dt_n[i][j])/dels_n[i][j];
     m_st_s[i][j]=m_st_s[i][j]-delt*((p_new[i][j]-p_new[i-1][j])*dels_ss[i][j]/ds_s[i][j] + (p_new_se[i][j]-p_new_sw[i][j])*dels_ts[i][j]/dt_s[i][j])/dels_s[i][j];

     Sm_st[i][j]=m_st_e[i][j]*dels_e[i][j]+m_st_n[i][j]*dels_n[i][j]-m_st_s[i][j]*dels_s[i][j]-m_st_w[i][j]*dels_w[i][j];

    if(fabs(Sm_st[i][j])>max_Sm_st)
	max_Sm_st=fabs(Sm_st[i][j]);

     p_old[i][j]=p_old[i][j]+p_new[i][j];

    }

    // B.C 
    for(i=0;i<=Ny-1;i++)
     {	   

       if(y[i][0]>=D_f)
        {
          p_new[i][0]=0.0;
          p_new[i][Nx-1]=0.0;

          p_old[i][0]=0.0;
          p_old[i][Nx-1]=0.0;
        }
       else
        {	   
          p_new[i][0]=p_new[i][1];
          p_new[i][Nx-1]=p_new[i][Nx-2];

          p_old[i][0]=p_old[i][1];
          p_old[i][Nx-1]=p_old[i][Nx-2];
        }
     }

    for(j=0;j<=Nx-1;j++)
     {

        ny = fabs(cos(alpha_e[Ny_s-2][j]));
        nx = fabs(sin(alpha_e[Ny_s-2][j]));

        p_new[0][j]=p_new[1][j];//-ds_s[1][j]*((v_s_new[Ny_s-1][j]-2.0*v_s[Ny_s-1][j]+v_s_old[Ny_s-1][j])*ny +(u_s_new[Ny_s-1][j]-2.0*u_s[Ny_s-1][j]+u_s_old[Ny_s-1][j])*nx)/pow(delt,2);
        p_new[Ny-1][j]=p_new[Ny-2][j];

        p_old[0][j]=p_old[1][j];//-ds_s[1][j]*((v_s_new[Ny_s-1][j]-2.0*v_s[Ny_s-1][j]+v_s_old[Ny_s-1][j])*ny +(u_s_new[Ny_s-1][j]-2.0*u_s[Ny_s-1][j]+u_s_old[Ny_s-1][j])*nx)/pow(delt,2);
        p_old[Ny-1][j]=p_old[Ny-2][j];
     }

   }*/
   
   

   /* Calculating velocities */
   for(j=1;j<=Nx-2;j++)
   for(i=1;i<=Ny-2;i++)
    {   

     pe=(del_v[i][j]*p_old[i][j+1]+del_v[i][j+1]*p_old[i][j])/(del_v[i][j]+del_v[i][j+1]);
     pw=(del_v[i][j]*p_old[i][j-1]+del_v[i][j-1]*p_old[i][j])/(del_v[i][j]+del_v[i][j-1]);
     pn=(del_v[i][j]*p_old[i+1][j]+del_v[i+1][j]*p_old[i][j])/(del_v[i][j]+del_v[i+1][j]);
     ps=(del_v[i][j]*p_old[i-1][j]+del_v[i-1][j]*p_old[i][j])/(del_v[i][j]+del_v[i-1][j]);

     u_new[i][j]=u_b[i][j]+(delt/(rho*del_v[i][j]))*(-pe*dels_xe[i][j]-pn*dels_xn[i][j]+pw*dels_xw[i][j]+ps*dels_xs[i][j]);
     v_new[i][j]=v_b[i][j]+(delt/(rho*del_v[i][j]))*(-pe*dels_ye[i][j]-pn*dels_yn[i][j]+pw*dels_yw[i][j]+ps*dels_ys[i][j]);

    }

   	for(j=1;j<=Nx-2;j++)
   	for(i=1;i<=Ny-2;i++)
    	 {   
		error=fabs(u_new[i][j]-u[i][j]);
     		if(error>maxerror)
			maxerror=error;

		error=fabs(v_new[i][j]-v[i][j]);
     		if(error>maxerror)
			maxerror=error;

		error=fabs(p_old[i][j]-p[i][j]);
     		if(error>maxerror)
			maxerror=error;
         }

	for(i=1;i<=(Ny-2);i++)
	for(j=1;j<=(Nx-2);j++)
	  u[i][j]=wu_vel*u_new[i][j]+(1.0-wu_vel)*u[i][j];

	for(i=1;i<=(Ny-2);i++)
	for(j=1;j<=(Nx-2);j++)
	  v[i][j]=wv_vel*v_new[i][j]+(1.0-wv_vel)*v[i][j];

   	for(i=1;i<=(Ny-2);i++)
   	for(j=1;j<=(Nx-2);j++)
     	  p[i][j]=w_p*p_old[i][j]+(1.0-w_p)*p[i][j];

}
