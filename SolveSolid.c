#include "DesignAndFlowData.h"

void solve_solid()
{

char* name = "solid_t";
char* extension = ".dat";
char fileSpec[strlen(name)+strlen(extension)+5];

char* name_stress = "stress_t";
char fileSpec2[strlen(name_stress)+strlen(extension)+5];

double wss;

if(iter==0)
{
for(i=0;i<Ny_s;i++)
for(j=0;j<Nx_s;j++)
 {

     u_s_old[i][j]=u_s[i][j];
     v_s_old[i][j]=v_s[i][j];

     u_s[i][j]=dx_new[i][j];
     v_s[i][j]=dy_new[i][j];

     dudt_old[i][j]=dudt[i][j];
     dvdt_old[i][j]=dvdt[i][j];

     dudt[i][j]=dxdt_new[i][j];
     dvdt[i][j]=dydt_new[i][j];

 }

 for(i=0;i<Ny_s;i++)
 for(j=0;j<Nx_s;j++)
  {

     dx[i][j]=u_s[i][j];
     dy[i][j]=v_s[i][j];

     dxdt[i][j]=dudt[i][j];
     dydt[i][j]=dvdt[i][j];

  }


}

it=0;
maxerror=1.0;
while(maxerror>1e-4)
{
maxerror=0.0;

for(i=0;i<Ny_s-0;i++)
for(j=0;j<Nx_s-0;j++)
  {
      dx_bar[i][j]=dx[i][j]+dT*(dxdt[i][j]-(3.0*dx[i][j]-4.0*u_s[i][j]+u_s_old[i][j])/(2.0*delt));
      dy_bar[i][j]=dy[i][j]+dT*(dydt[i][j]-(3.0*dy[i][j]-4.0*v_s[i][j]+v_s_old[i][j])/(2.0*delt));
  }

for(i=1;i<Ny_s-1;i++)
for(j=1;j<Nx_s-1;j++)
  {	
      	a0 = pow(dT,2.0)/(rho_s*del_v_s[i][j]);
      	a1 = pow(dT,1.0)/(rho_s*del_v_s[i][j]);

	/* Calculation of non-linear strains */
	dx_ne=0.25*(dx_bar[i][j]+dx_bar[i][j+1]+dx_bar[i+1][j]+dx_bar[i+1][j+1]);
	dy_ne=0.25*(dy_bar[i][j]+dy_bar[i][j+1]+dy_bar[i+1][j]+dy_bar[i+1][j+1]);

	dx_se=0.25*(dx_bar[i][j]+dx_bar[i][j+1]+dx_bar[i-1][j]+dx_bar[i-1][j+1]);
	dy_se=0.25*(dy_bar[i][j]+dy_bar[i][j+1]+dy_bar[i-1][j]+dy_bar[i-1][j+1]);

	dx_nw=0.25*(dx_bar[i][j]+dx_bar[i][j-1]+dx_bar[i+1][j]+dx_bar[i+1][j-1]);
	dy_nw=0.25*(dy_bar[i][j]+dy_bar[i][j-1]+dy_bar[i+1][j]+dy_bar[i+1][j-1]);

	dx_sw=0.25*(dx_bar[i][j]+dx_bar[i][j-1]+dx_bar[i-1][j]+dx_bar[i-1][j-1]);
	dy_sw=0.25*(dy_bar[i][j]+dy_bar[i][j-1]+dy_bar[i-1][j]+dy_bar[i-1][j-1]);

	dv_dx_e = ((dy_bar[i][j+1]-dy_bar[i][j])*cos(alpha_e[i][j])/ds_s_e[i][j] - (dy_ne-dy_se)*sin(beta_e[i][j])/dt_s_e[i][j])/cos(alpha_e[i][j]-beta_e[i][j]);
	dv_dy_e = ((dy_bar[i][j+1]-dy_bar[i][j])*sin(alpha_e[i][j])/ds_s_e[i][j] + (dy_ne-dy_se)*cos(beta_e[i][j])/dt_s_e[i][j])/cos(alpha_e[i][j]-beta_e[i][j]);
	du_dx_e = ((dx_bar[i][j+1]-dx_bar[i][j])*cos(alpha_e[i][j])/ds_s_e[i][j] - (dx_ne-dx_se)*sin(beta_e[i][j])/dt_s_e[i][j])/cos(alpha_e[i][j]-beta_e[i][j]);
	du_dy_e = ((dx_bar[i][j+1]-dx_bar[i][j])*sin(alpha_e[i][j])/ds_s_e[i][j] + (dx_ne-dx_se)*cos(beta_e[i][j])/dt_s_e[i][j])/cos(alpha_e[i][j]-beta_e[i][j]);

	dv_dx_w = ((dy_bar[i][j-1]-dy_bar[i][j])*cos(alpha_w[i][j])/ds_s_w[i][j] - (dy_sw-dy_nw)*sin(beta_w[i][j])/dt_s_w[i][j])/cos(alpha_w[i][j]-beta_w[i][j]);
	dv_dy_w = ((dy_bar[i][j-1]-dy_bar[i][j])*sin(alpha_w[i][j])/ds_s_w[i][j] + (dy_sw-dy_nw)*cos(beta_w[i][j])/dt_s_w[i][j])/cos(alpha_w[i][j]-beta_w[i][j]);
	du_dx_w = ((dx_bar[i][j-1]-dx_bar[i][j])*cos(alpha_w[i][j])/ds_s_w[i][j] - (dx_sw-dx_nw)*sin(beta_w[i][j])/dt_s_w[i][j])/cos(alpha_w[i][j]-beta_w[i][j]);
	du_dy_w = ((dx_bar[i][j-1]-dx_bar[i][j])*sin(alpha_w[i][j])/ds_s_w[i][j] + (dx_sw-dx_nw)*cos(beta_w[i][j])/dt_s_w[i][j])/cos(alpha_w[i][j]-beta_w[i][j]);

	dv_dx_n = ((dy_bar[i+1][j]-dy_bar[i][j])*cos(alpha_n[i][j])/ds_s_n[i][j] - (dy_nw-dy_ne)*sin(beta_n[i][j])/dt_s_n[i][j])/cos(alpha_n[i][j]-beta_n[i][j]);
	dv_dy_n = ((dy_bar[i+1][j]-dy_bar[i][j])*sin(alpha_n[i][j])/ds_s_n[i][j] + (dy_nw-dy_ne)*cos(beta_n[i][j])/dt_s_n[i][j])/cos(alpha_n[i][j]-beta_n[i][j]);
	du_dx_n = ((dx_bar[i+1][j]-dx_bar[i][j])*cos(alpha_n[i][j])/ds_s_n[i][j] - (dx_nw-dx_ne)*sin(beta_n[i][j])/dt_s_n[i][j])/cos(alpha_n[i][j]-beta_n[i][j]);
	du_dy_n = ((dx_bar[i+1][j]-dx_bar[i][j])*sin(alpha_n[i][j])/ds_s_n[i][j] + (dx_nw-dx_ne)*cos(beta_n[i][j])/dt_s_n[i][j])/cos(alpha_n[i][j]-beta_n[i][j]);

	dv_dx_s = ((dy_bar[i-1][j]-dy_bar[i][j])*cos(alpha_s[i][j])/ds_s_s[i][j] - (dy_se-dy_sw)*sin(beta_s[i][j])/dt_s_s[i][j])/cos(alpha_s[i][j]-beta_s[i][j]);
	dv_dy_s = ((dy_bar[i-1][j]-dy_bar[i][j])*sin(alpha_s[i][j])/ds_s_s[i][j] + (dy_se-dy_sw)*cos(beta_s[i][j])/dt_s_s[i][j])/cos(alpha_s[i][j]-beta_s[i][j]);
	du_dx_s = ((dx_bar[i-1][j]-dx_bar[i][j])*cos(alpha_s[i][j])/ds_s_s[i][j] - (dx_se-dx_sw)*sin(beta_s[i][j])/dt_s_s[i][j])/cos(alpha_s[i][j]-beta_s[i][j]);
	du_dy_s = ((dx_bar[i-1][j]-dx_bar[i][j])*sin(alpha_s[i][j])/ds_s_s[i][j] + (dx_se-dx_sw)*cos(beta_s[i][j])/dt_s_s[i][j])/cos(alpha_s[i][j]-beta_s[i][j]);

	e_xx_e=du_dx_e+0.5*(du_dx_e*du_dx_e+dv_dx_e*dv_dx_e);
	e_yy_e=dv_dy_e+0.5*(du_dy_e*du_dy_e+dv_dy_e*dv_dy_e);
	e_xy_e=du_dy_e+dv_dx_e+(du_dx_e*du_dy_e+dv_dx_e*dv_dy_e);

         e_xx_w=du_dx_w+0.5*(du_dx_w*du_dx_w+dv_dx_w*dv_dx_w);
	e_yy_w=dv_dy_w+0.5*(du_dy_w*du_dy_w+dv_dy_w*dv_dy_w);
	e_xy_w=du_dy_w+dv_dx_w+(du_dx_w*du_dy_w+dv_dx_w*dv_dy_w);

         e_xx_n=du_dx_n+0.5*(du_dx_n*du_dx_n+dv_dx_n*dv_dx_n);
	e_yy_n=dv_dy_n+0.5*(du_dy_n*du_dy_n+dv_dy_n*dv_dy_n);
	e_xy_n=du_dy_n+dv_dx_n+(du_dx_n*du_dy_n+dv_dx_n*dv_dy_n); 

         e_xx_s=du_dx_s+0.5*(du_dx_s*du_dx_s+dv_dx_s*dv_dx_s);
	e_yy_s=dv_dy_s+0.5*(du_dy_s*du_dy_s+dv_dy_s*dv_dy_s);
	e_xy_s=du_dy_s+dv_dx_s+(du_dx_s*du_dy_s+dv_dx_s*dv_dy_s);

         s_xx_e = lamda*(1.0-nu)*e_xx_e+lamda*nu*e_yy_e;
         s_xx_w = lamda*(1.0-nu)*e_xx_w+lamda*nu*e_yy_w;
         s_xx_n = lamda*(1.0-nu)*e_xx_n+lamda*nu*e_yy_n;
	s_xx_s = lamda*(1.0-nu)*e_xx_s+lamda*nu*e_yy_s;  

         s_yy_e = lamda*nu*e_xx_e+lamda*(1.0-nu)*e_yy_e; 
	s_yy_w = lamda*nu*e_xx_w+lamda*(1.0-nu)*e_yy_w;   
	s_yy_n = lamda*nu*e_xx_n+lamda*(1.0-nu)*e_yy_n; 
	s_yy_s = lamda*nu*e_xx_s+lamda*(1.0-nu)*e_yy_s; 

	s_xy_e = mu_s*e_xy_e;
	s_xy_w = mu_s*e_xy_w;
	s_xy_n = mu_s*e_xy_n;
	s_xy_s = mu_s*e_xy_s;

	fe = ((1.0+du_dx_e)*s_xx_e+du_dy_e*s_xy_e)*dt_s_e[i][j]*cos(alpha_e[i][j])+((1.0+du_dx_e)*s_xy_e+du_dy_e*s_yy_e)*dt_s_e[i][j]*sin(alpha_e[i][j]);
	fw = ((1.0+du_dx_w)*s_xx_w+du_dy_w*s_xy_w)*dt_s_w[i][j]*cos(alpha_w[i][j])+((1.0+du_dx_w)*s_xy_w+du_dy_w*s_yy_w)*dt_s_w[i][j]*sin(alpha_w[i][j]);
	fn = ((1.0+du_dx_n)*s_xx_n+du_dy_n*s_xy_n)*dt_s_n[i][j]*cos(alpha_n[i][j])+((1.0+du_dx_n)*s_xy_n+du_dy_n*s_yy_n)*dt_s_n[i][j]*sin(alpha_n[i][j]);
         fs = ((1.0+du_dx_s)*s_xx_s+du_dy_s*s_xy_s)*dt_s_s[i][j]*cos(alpha_s[i][j])+((1.0+du_dx_s)*s_xy_s+du_dy_s*s_yy_s)*dt_s_s[i][j]*sin(alpha_s[i][j]);	

	dxdt_new[i][j] = ( dxdt[i][j] + a1*(fe+fw+fn+fs) - (0.5*dT/delt)*(-4.0*dudt[i][j]+dudt_old[i][j]) )/(1.0+1.5*dT/delt);
	
	/* update displacements */
	dx_new[i][j] = ( dx[i][j] + dT*( dxdt_new[i][j] - 0.5*(-4.0*u_s[i][j]+u_s_old[i][j])/delt )  + 0.5*a0*( fe+fw+fn+fs - rho_s*del_v_s[i][j]*(3.0*dxdt_new[i][j]-4.0*dudt[i][j]+dudt_old[i][j])/(2.0*delt) ) )/(1.0 + 1.5*dT/delt);

	/*     ********************* U solved *************************      */

	fe = (dv_dx_e*s_xx_e+(1.0+dv_dy_e)*s_xy_e)*dt_s_e[i][j]*cos(alpha_e[i][j])+(dv_dx_e*s_xy_e+(1.0+dv_dy_e)*s_yy_e)*dt_s_e[i][j]*sin(alpha_e[i][j]);
	fw = (dv_dx_w*s_xx_w+(1.0+dv_dy_w)*s_xy_w)*dt_s_w[i][j]*cos(alpha_w[i][j])+(dv_dx_w*s_xy_w+(1.0+dv_dy_w)*s_yy_w)*dt_s_w[i][j]*sin(alpha_w[i][j]);
	fn = (dv_dx_n*s_xx_n+(1.0+dv_dy_n)*s_xy_n)*dt_s_n[i][j]*cos(alpha_n[i][j])+(dv_dx_n*s_xy_n+(1.0+dv_dy_n)*s_yy_n)*dt_s_n[i][j]*sin(alpha_n[i][j]);
	fs = (dv_dx_s*s_xx_s+(1.0+dv_dy_s)*s_xy_s)*dt_s_s[i][j]*cos(alpha_s[i][j])+(dv_dx_s*s_xy_s+(1.0+dv_dy_s)*s_yy_s)*dt_s_s[i][j]*sin(alpha_s[i][j]);         

	dydt_new[i][j] = ( dydt[i][j] + a1*(fe+fw+fn+fs) - (0.5*dT/delt)*(-4.0*dvdt[i][j]+dvdt_old[i][j]) )/(1.0 + 1.5*dT/delt);
    
	/* update displacements */
	dy_new[i][j] = ( dy[i][j] + dT*( dydt_new[i][j] - 0.5*(-4.0*v_s[i][j]+v_s_old[i][j])/delt ) + 0.5*a0*( fe+fw+fn+fs - rho_s*del_v_s[i][j]*(3.0*dydt_new[i][j]-4.0*dvdt[i][j]+dvdt_old[i][j])/(2.0*delt) ) )/(1.0 + 1.5*dT/delt);

  }

/* Apply B.C  */
for(i=0;i<Ny_s;i++)
  {
     dx_new[i][0]=0.0;
     dy_new[i][0]=0.0;
  }

// Top and Bottom v
for(j=1;j<Nx_s-1;j++)
  {
     i=Ny_s-1;

     dv_dx_s = ((dy_bar[i-1][j]-dy_bar[i][j])*cos(alpha_s[i-1][j])/ds_s_s[i-1][j] - 0.5*(dy_bar[i][j+1]-dy_bar[i][j-1])*sin(beta_s[i-1][j])/dt_s_s[i-1][j])/cos(alpha_s[i-1][j]-beta_s[i-1][j]);
     dv_dy_s = ((dy_bar[i-1][j]-dy_bar[i][j])*sin(alpha_s[i-1][j])/ds_s_s[i-1][j] + 0.5*(dy_bar[i][j+1]-dy_bar[i][j-1])*cos(beta_s[i-1][j])/dt_s_s[i-1][j])/cos(alpha_s[i-1][j]-beta_s[i-1][j]);
     du_dx_s = ((dx_bar[i-1][j]-dx_bar[i][j])*cos(alpha_s[i-1][j])/ds_s_s[i-1][j] - 0.5*(dx_bar[i][j+1]-dx_bar[i][j-1])*sin(beta_s[i-1][j])/dt_s_s[i-1][j])/cos(alpha_s[i-1][j]-beta_s[i-1][j]);
     du_dy_s = ((dx_bar[i-1][j]-dx_bar[i][j])*sin(alpha_s[i-1][j])/ds_s_s[i-1][j] + 0.5*(dx_bar[i][j+1]-dx_bar[i][j-1])*cos(beta_s[i-1][j])/dt_s_s[i-1][j])/cos(alpha_s[i-1][j]-beta_s[i-1][j]);    

     e_xx_s=du_dx_s+0.5*(du_dx_s*du_dx_s+dv_dx_s*dv_dx_s);
     e_yy_s=dv_dy_s+0.5*(du_dy_s*du_dy_s+dv_dy_s*dv_dy_s);
     e_xy_s=du_dy_s+dv_dx_s+(du_dx_s*du_dy_s+dv_dx_s*dv_dy_s);

     s_xy_s = mu_s*e_xy_s;
     s_yy_s = lamda*(1.0-nu)*e_yy_s+lamda*nu*e_xx_s;
     s_xx_s = lamda*(1.0-nu)*e_xx_s+lamda*nu*e_yy_s; 

     j_f = j;//(Nx-1 + 0.00001)*(1.0-(Nx_s-1.0-j)/(Nx_s-1.0));

     theta = acos( (y_s[i][j]-y_s[i][j+1])/(sqrt( pow(x_s_old[i][j+1]-x_s_old[i][j],2) + pow(y_s_old[i][j+1]-y_s_old[i][j],2) )) );     
     Fy = p_old[0][j_f];

     fs = -(dv_dx_s*s_xy_s+(1.0+dv_dy_s)*s_yy_s)*dt_s_s[i-1][j] + ( -p_old[0][j_f]*sin(theta) + mu*(v_b[1][j_f]-v_b[0][j_f])/(y[1][j_f]-y[0][j_f]) )*dt_s_s[i-1][j];

     dydt_new[i][j] = ( dydt[i][j] + a1*(fs) - (0.5*dT/delt)*(-4.0*dvdt[i][j]+dvdt_old[i][j]) )/(1.0 + 1.5*dT/delt);
     dy_new[i][j] = ( dy[i][j] + dT*( dydt_new[i][j] - 0.5*(-4.0*v_s[i][j]+v_s_old[i][j])/delt ) + 0.5*a0*( fs - rho_s*del_v_s[i][j]*(3.0*dydt_new[i][j]-4.0*dvdt[i][j]+dvdt_old[i][j])/(2.0*delt) ) )/(1.0 + 1.5*dT/delt);

     fs = -((1.0+du_dx_s)*s_xy_s+du_dy_s*s_yy_s)*dt_s_s[i-1][j] + ( -p_old[0][j_f]*cos(theta) + mu*(u_b[1][j_f]-u_b[0][j_f])/(y[1][j_f]-y[0][j_f]) )*dt_s_s[i-1][j];	
     dxdt_new[i][j] = ( dxdt[i][j] + a1*(fs) - (0.5*dT/delt)*(-4.0*dudt[i][j]+dudt_old[i][j]) )/(1.0+1.5*dT/delt);
     dx_new[i][j] = ( dx[i][j] + dT*( dxdt_new[i][j] - 0.5*(-4.0*u_s[i][j]+u_s_old[i][j])/delt )  + 0.5*a0*( fs - rho_s*del_v_s[i][j]*(3.0*dxdt_new[i][j]-4.0*dudt[i][j]+dudt_old[i][j])/(2.0*delt) ) )/(1.0 + 1.5*dT/delt);

     i=0;
     dv_dx_n = ((dy_bar[i+1][j]-dy_bar[i][j])*cos(alpha_n[i+1][j])/ds_s_n[i+1][j] - 0.5*(dy_bar[i][j-1]-dy_bar[i][j+1])*sin(beta_n[i+1][j])/dt_s_n[i+1][j])/cos(alpha_n[i+1][j]-beta_n[i+1][j]);
     dv_dy_n = ((dy_bar[i+1][j]-dy_bar[i][j])*sin(alpha_n[i+1][j])/ds_s_n[i+1][j] + 0.5*(dy_bar[i][j-1]-dy_bar[i][j+1])*cos(beta_n[i+1][j])/dt_s_n[i+1][j])/cos(alpha_n[i+1][j]-beta_n[i+1][j]);
     du_dx_n = ((dx_bar[i+1][j]-dx_bar[i][j])*cos(alpha_n[i+1][j])/ds_s_n[i+1][j] - 0.5*(dx_bar[i][j-1]-dx_bar[i][j+1])*sin(beta_n[i+1][j])/dt_s_n[i+1][j])/cos(alpha_n[i+1][j]-beta_n[i+1][j]);
     du_dy_n = ((dx_bar[i+1][j]-dx_bar[i][j])*sin(alpha_n[i+1][j])/ds_s_n[i+1][j] + 0.5*(dx_bar[i][j-1]-dx_bar[i][j+1])*cos(beta_n[i+1][j])/dt_s_n[i+1][j])/cos(alpha_n[i+1][j]-beta_n[i+1][j]);

     e_xx_n=du_dx_n+0.5*(du_dx_n*du_dx_n+dv_dx_n*dv_dx_n);
     e_yy_n=dv_dy_n+0.5*(du_dy_n*du_dy_n+dv_dy_n*dv_dy_n);
     e_xy_n=du_dy_n+dv_dx_n+(du_dx_n*du_dy_n+dv_dx_n*dv_dy_n); 

     s_xx_n = lamda*(1.0-nu)*e_xx_n+lamda*nu*e_yy_n;
     s_yy_n = lamda*nu*e_xx_n+lamda*(1.0-nu)*e_yy_n;
     s_xy_n = mu_s*e_xy_n;

     fn = (dv_dx_n*s_xy_n+(1.0+dv_dy_n)*s_yy_n)*dt_s_n[i+1][j] ;
     dydt_new[i][j] = ( dydt[i][j] + a1*(fn) - (0.5*dT/delt)*(-4.0*dvdt[i][j]+dvdt_old[i][j]) )/(1.0 + 1.5*dT/delt);
     dy_new[i][j] = ( dy[i][j] + dT*( dydt_new[i][j] - 0.5*(-4.0*v_s[i][j]+v_s_old[i][j])/delt ) + 0.5*a0*( fn - rho_s*del_v_s[i][j]*(3.0*dydt_new[i][j]-4.0*dvdt[i][j]+dvdt_old[i][j])/(2.0*delt) ) )/(1.0 + 1.5*dT/delt);

     fn = ((1.0+du_dx_n)*s_xy_n+du_dy_n*s_yy_n)*dt_s_n[i+1][j];
     dxdt_new[i][j] = ( dxdt[i][j] + a1*(fn) - (0.5*dT/delt)*(-4.0*dudt[i][j]+dudt_old[i][j]) )/(1.0+1.5*dT/delt);
     dx_new[i][j] = ( dx[i][j] + dT*( dxdt_new[i][j] - 0.5*(-4.0*u_s[i][j]+u_s_old[i][j])/delt )  + 0.5*a0*( fn - rho_s*del_v_s[i][j]*(3.0*dxdt_new[i][j]-4.0*dudt[i][j]+dudt_old[i][j])/(2.0*delt) ) )/(1.0 + 1.5*dT/delt);
   
  }

// Side BC
for(i=0;i<Ny_s;i++)
  {
     dx_new[i][Nx_s-1]=0.0;
     dy_new[i][Nx_s-1]=0.0;
  }

for(i=1;i<Ny_s-1;i++)
for(j=1;j<Nx_s-1;j++)
{

  error=fabs(dxdt_new[i][j]-dxdt[i][j])/dT;
  if(error>maxerror)
	maxerror=error;

  error=fabs(dydt_new[i][j]-dydt[i][j])/dT;
  if(error>maxerror)
	maxerror=error;

}

for(i=0;i<Ny_s;i++)
for(j=0;j<Nx_s;j++)
 {

   dx_old[i][j]=dx[i][j];
   dy_old[i][j]=dy[i][j];

   dx[i][j] = dx_new[i][j];
   dy[i][j] = dy_new[i][j];

   dxdt[i][j] = dxdt_new[i][j];
   dydt[i][j] = dydt_new[i][j];

 }
it++;
} // end of psuedo time loop

/* update interface variables */
for(j=0;j<Nx;j++)
   {

     j_s = j;//(Nx_s-1 + 0.00001)*(1.0-(Nx-1.0-j)/(Nx-1.0));

     x_s[Ny_s-1][j_s]=x_init[Ny_s-1][j_s]+dx_new[Ny_s-1][j_s];
     y_s[Ny_s-1][j_s]=y_init[Ny_s-1][j_s]+dy_new[Ny_s-1][j_s];

     tol = fabs(x_s[Ny_s-1][j_s] - x[0][j]);
     if(tol>maxtol)
       maxtol = tol;

     tol = fabs(y_s[Ny_s-1][j_s] - y[0][j]);
     if(tol>maxtol)
       maxtol = tol;

     x[0][j]=w_s_u*x_s[Ny_s-1][j_s] + (1.0-w_s_u)*x_s_old[Ny_s-1][j_s];
     y[0][j]=w_s_v*y_s[Ny_s-1][j_s] + (1.0-w_s_v)*y_s_old[Ny_s-1][j_s];

     x_s_old[Ny_s-1][j_s]=x_s[Ny_s-1][j_s];
     y_s_old[Ny_s-1][j_s]=y_s[Ny_s-1][j_s];
     
     umesh[0][j]=(x[0][j]-x_old[0][j])/delt;
     vmesh[0][j]=(y[0][j]-y_old[0][j])/delt;

   }

if(time>=1.0)
printf("%lf\t%lf\n",y_s[Ny_s-1][(Nx_s-1)/2],y[0][(Nx-1)/2]);

}
