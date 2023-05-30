//#ifndef _DESIGNANDFLOWDATA_H_
//#define _DESIGNANDFLOWDATA_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define Nx 121
#define Ny 121

#define D 1.0
#define L 1.0

#define Dt 0.001
#define dT 0.00001

#define D_f 0.9
#define PI 3.14159265

#define Nx_s 121
#define Ny_s 5

#define Le 1.0
#define H 0.002

#define rho_s 500.0
#define nu 0.0
#define cd 0.0
#define E 25000.0

/* Fluid Variables */

double eta[Ny][Nx],zeta[Ny][Nx],x[Ny][Nx],y[Ny][Nx],xnew[Ny][Nx],ynew[Ny][Nx],J[Ny][Nx];
double u[Ny][Nx],u_new[Ny][Nx],v[Ny][Nx],v_new[Ny][Nx],p[Ny][Nx],p_new[Ny][Nx],p_old[Ny][Nx],u_dash[Ny][Nx],v_dash[Ny][Nx],p_dash[Ny][Nx];
double Sm_st[Ny][Nx],Sm_dash[Ny][Nx],m_st_e[Ny][Nx],m_st_w[Ny][Nx],m_st_n[Ny][Nx],m_st_s[Ny][Nx],m_dash,max_Sm_st,ap,Sm,rad;
double del_eta,del_zeta,alpha,s;
double zeta_x,zeta_xx,zeta_y,zeta_yy,eta_x,eta_xx,eta_y,eta_yy,d1,d2;
double x_zeta,x_eta,y_zeta,y_eta;

double dt_e[Ny][Nx],dt_w[Ny][Nx],dt_n[Ny][Nx],dt_s[Ny][Nx],ds_e[Ny][Nx],ds_w[Ny][Nx],ds_n[Ny][Nx],ds_s[Ny][Nx],del_v[Ny][Nx];
double x_ne[Ny][Nx],x_nw[Ny][Nx],x_se[Ny][Nx],x_sw[Ny][Nx],y_ne[Ny][Nx],y_nw[Ny][Nx],y_se[Ny][Nx],y_sw[Ny][Nx];
double dels_sn[Ny][Nx],dels_se[Ny][Nx],dels_sw[Ny][Nx],dels_ss[Ny][Nx],dels_tn[Ny][Nx],dels_ts[Ny][Nx],dels_tw[Ny][Nx],dels_te[Ny][Nx],dels_e[Ny][Nx],dels_w[Ny][Nx],dels_n[Ny][Nx],dels_s[Ny][Nx];
double dels_xn[Ny][Nx],dels_xe[Ny][Nx],dels_xw[Ny][Nx],dels_xs[Ny][Nx],dels_yn[Ny][Nx],dels_ys[Ny][Nx],dels_yw[Ny][Nx],dels_ye[Ny][Nx];
double u_ne[Ny][Nx],u_nw[Ny][Nx],u_sw[Ny][Nx],u_se[Ny][Nx],v_ne[Ny][Nx],v_nw[Ny][Nx],v_sw[Ny][Nx],v_se[Ny][Nx],u_b[Ny][Nx],v_b[Ny][Nx],au,du,av,dv;
double c_th_tx,c_th_ty,c_th_sx,c_th_sy,m_e,m_w,m_n,m_s,ue_p,uw_p,un_p,us_p,ue_m,uw_m,un_m,us_m,ve_p,vw_p,vn_p,vs_p,ve_m,vw_m,vn_m,vs_m,pe,pw,pn,ps,p_ne[Ny][Nx],p_nw[Ny][Nx],p_se[Ny][Nx],p_sw[Ny][Nx];
double p_new_ne[Ny][Nx],p_new_nw[Ny][Nx],p_new_se[Ny][Nx],p_new_sw[Ny][Nx];
double we1_p,we1_m,ww1_p,ww1_m,wn1_p,wn1_m,ws1_p,ws1_m,we2_p,we2_m,ww2_p,ww2_m,wn2_p,wn2_m,ws2_p,ws2_m,we3_p,we3_m,ww3_p,ww3_m,wn3_p,wn3_m,ws3_p,ws3_m;
double ue,uw,un,us,ve,vw,vn,vs;
double umesh_e,umesh_w,umesh_n,umesh_s,vmesh_e,vmesh_w,vmesh_n,vmesh_s;

double delt,rho,mu,u_inlet_max,alp_u,alp_v,alp_p,wu_vel,wv_vel,w_p,error,Sm_st_new,maxerror_u,maxerror_p,ub,maxerror_v,vb,theta;
double L_eta,L_zeta,maxerror,rho,mu,u_inlet_max,per,maxerror_disp,x_oldo[Nx];


/* Solid Variables */
double x_s[Ny_s][Nx_s],y_s[Ny_s][Nx_s],B[2*Nx_s*Ny_s],X[2*Nx_s*Ny_s],a0,matrix[2*Nx_s*Ny_s][2*Nx_s*Ny_s],d[2*Nx_s*Ny_s],sum,maxerror,del_eta_s,del_zeta_s;
double eta_s[Ny_s][Nx_s],zeta_s[Ny_s][Nx_s],x_init[Ny_s][Nx_s],y_init[Ny_s][Nx_s];
double u_s[Ny_s][Nx_s],u_s_new[Ny_s][Nx_s],v_s[Ny_s][Nx_s],v_s_new[Ny_s][Nx_s];
double alpha,error,w_s_u,w_s_v;
double time,mu_s,lamda,d1,d2;
double dt_s_e[Ny_s][Nx_s],dt_s_w[Ny_s][Nx_s],dt_s_n[Ny_s][Nx_s],dt_s_s[Ny_s][Nx_s],ds_s_e[Ny_s][Nx_s],ds_s_w[Ny_s][Nx_s],ds_s_n[Ny_s][Nx_s],ds_s_s[Ny_s][Nx_s],del_v_s[Ny_s][Nx_s];
double x_s_ne[Ny_s][Nx_s],x_s_nw[Ny_s][Nx_s],x_s_se[Ny_s][Nx_s],x_s_sw[Ny_s][Nx_s],y_s_ne[Ny_s][Nx_s],y_s_nw[Ny_s][Nx_s],y_s_se[Ny_s][Nx_s],y_s_sw[Ny_s][Nx_s];
double u_s_ne,u_s_nw,u_s_sw,u_s_se,v_s_ne,v_s_nw,v_s_sw,v_s_se,u_s_old[Ny_s][Nx_s],v_s_old[Ny_s][Nx_s],au,du,av,dv;
double alpha_e[Ny_s][Nx_s],alpha_w[Ny_s][Nx_s],alpha_n[Ny_s][Nx_s],alpha_s[Ny_s][Nx_s],beta_e[Ny_s][Nx_s],beta_w[Ny_s][Nx_s],beta_n[Ny_s][Nx_s],beta_s[Ny_s][Nx_s];
double Fx,Fy,dv_dx_e,dv_dx_w,dv_dx_n,dv_dx_s,dv_dy_e,dv_dy_w,dv_dy_n,dv_dy_s;
double ap_e_u_s[Ny_s][Nx_s],ap_w_u_s[Ny_s][Nx_s],ap_n_u_s[Ny_s][Nx_s],ap_s_u_s[Ny_s][Nx_s],f_u_s[Ny_s][Nx_s],ae_u_s[Ny_s][Nx_s],aw_u_s[Ny_s][Nx_s],as_u_s[Ny_s][Nx_s],an_u_s[Ny_s][Nx_s],ap_e_v_s[Ny_s][Nx_s],ap_w_v_s[Ny_s][Nx_s],ap_n_v_s[Ny_s][Nx_s],ap_s_v_s[Ny_s][Nx_s],f_v_s[Ny_s][Nx_s],ae_v_s[Ny_s][Nx_s],aw_v_s[Ny_s][Nx_s],as_v_s[Ny_s][Nx_s],an_v_s[Ny_s][Nx_s];
double ane_u_s[Ny_s][Nx_s],ase_u_s[Ny_s][Nx_s],anw_u_s[Ny_s][Nx_s],asw_u_s[Ny_s][Nx_s],ane_v_s[Ny_s][Nx_s],ase_v_s[Ny_s][Nx_s],anw_v_s[Ny_s][Nx_s],asw_v_s[Ny_s][Nx_s];
double du_dx_e,du_dx_w,du_dx_n,du_dx_s,du_dy_e,du_dy_w,du_dy_n,du_dy_s,fe,fw,fn,fs;
double e_xx_e,e_xx_w,e_xx_n,e_xx_s,e_yy_e,e_yy_w,e_yy_n,e_yy_s,e_xy_e,e_xy_w,e_xy_n,e_xy_s;
double x_s_old[Ny_s][Nx_s],y_s_old[Ny_s][Nx_s];
double dx_ne,dx_nw,dx_sw,dx_se,dy_ne,dy_nw,dy_sw,dy_se;
double s_xx_e,s_xx_w,s_xx_n,s_xx_s,s_yy_e,s_yy_w,s_yy_n,s_yy_s,s_xy_e,s_xy_w,s_xy_n,s_xy_s;
double dudt[Ny_s][Nx_s],dudt_new[Ny_s][Nx_s],dudt_old[Ny_s][Nx_s],dvdt[Ny_s][Nx_s],dvdt_new[Ny_s][Nx_s],dvdt_old[Ny_s][Nx_s],dxdt[Ny_s][Nx_s],dxdt_new[Ny_s][Nx_s],dydt[Ny_s][Nx_s],dydt_new[Ny_s][Nx_s];
double dx[Ny_s][Nx_s],dy[Ny_s][Nx_s],dx_bar[Ny_s][Nx_s],dy_bar[Ny_s][Nx_s],dx_old[Ny_s][Nx_s],dy_old[Ny_s][Nx_s],dx_new[Ny_s][Nx_s],dy_new[Ny_s][Nx_s];

/* Common Variables */
int    i,j,l1,l2,lst,it,it2,it3,ct,count,iter,k,j_f,j_s;
double umesh[Ny][Nx],vmesh[Ny][Nx],x_old[Ny][Nx],y_old[Ny][Nx],nx,ny,maxtol,tol,a1;

/* Fluid Subroutines */
void gridgen(int ,int , int ,int );
void geometric_parameters();
void define_flow_para();
void initial_cond_fluid();
void solve_fluid();
void bc_fluid();
void mesh_vel();

/* Solid Subroutines */
//void inverse();
void gridgen_solid(int ,int ,int ,int );
void initial_cond_solid();
void solve_solid();


