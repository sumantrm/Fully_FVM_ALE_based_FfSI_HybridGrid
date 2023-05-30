#include "DesignAndFlowData.h"

void bc_fluid()
{

    /* B.C */
    /* TOP AND BOTTOM */
    for(j=1;j<=Nx-2;j++)
     {
      
      u[Ny-1][j] = u_inlet_max*(1.0-cos(2*PI*time/5.0));
      u_b[Ny-1][j] = u_inlet_max*(1.0-cos(2*PI*time/5.0));

      ny = fabs(cos(alpha_e[Ny_s-2][j]));
      nx = fabs(sin(alpha_e[Ny_s-2][j]));

      p[0][j]=p[1][j];//-ds_s[1][j]*((v_s_new[Ny_s-1][j]-2.0*v_s[Ny_s-1][j]+v_s_old[Ny_s-1][j])*ny +(u_s_new[Ny_s-1][j]-2.0*u_s[Ny_s-1][j]+u_s_old[Ny_s-1][j])*nx)/pow(delt,2);
      p[Ny-1][j]=p[Ny-2][j];

      u[0][j]=umesh[0][j];
      v[0][j]=vmesh[0][j];

      u_b[0][j]=umesh[0][j];
      v_b[0][j]=vmesh[0][j];

     }

    /* INLET and OUTLET */
    for(i=0;i<=Ny-1;i++)
     {   	
       if(y[i][0]>=D_f)
        {

          u[i][0]=u[i][1];
	 v[i][0]=v[i][1];
	  
	 u_b[i][0]=u_b[i][1];
	 v_b[i][0]=v_b[i][1];

          p[i][0]=0.0;
          p[i][Nx-1]=0.0;

          v[i][Nx-1]=v[i][Nx-2];
          u[i][Nx-1]=u[i][Nx-2];

          v_b[i][Nx-1]=v_b[i][Nx-2];
          u_b[i][Nx-1]=u_b[i][Nx-2];

        }

       else
        {
          p[i][0]=p[i][1];
          p[i][Nx-1]=p[i][Nx-2];
        }

     }
   
}
