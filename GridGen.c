#include "DesignAndFlowData.h"

// Novel hybrid SEMM-MRR-AMM-based Grid Generation 
void gridgen(int i_start,int i_end, int j_start,int j_end)
{

double a[Ny][Nx],b[Ny][Nx],c[Ny][Nx],d[Ny][Nx],e[Ny][Nx];

for(i=i_start;i<=i_end;i++)
for(j=j_start;j<=j_end;j++)
 {
   x_old[i][j]=x[i][j];
   y_old[i][j]=y[i][j];
 }

// Direct Mesh Motion of fluid cells near solid region (20% considered here)
for(i=i_start+1;i<i_start+per*(Ny-1);i++)
for(j=j_start+1;j<j_end;j++)
 {
   j_s = j;
   x[i][j] = x_s[Ny_s-1][j_s];
   y[i][j] = i*del_eta + dy_new[Ny_s-1][j_s];     
 }

maxerror_disp = 0.0;
for(j=j_start+1;j<j_end;j++)	
 if(abs(x[0][j]-x_oldo[j])>maxerror_disp)
    maxerror_disp = abs(x[0][j]-x_oldo[j]);
         
if(maxerror_disp>0.01) // 1% cut-off
 {
     
	for(j=j_start+1;j<j_end;j++)
	x_oldo[j] = x[0][j];
	              
	/* Finding x and y by solving laplacian equation */
	maxerror=1.0;
	
	while(maxerror>1e-8)
	{

		maxerror=0.0;
		for(i=i_start+per*(Ny-1);i<i_end;i++)
			for(j=j_start+1;j<j_end;j++)
			{

				x_zeta=0.5*(x[i+1][j]-x[i-1][j])/del_zeta;
				x_eta=0.5*(x[i][j+1]-x[i][j-1])/del_eta;
				y_zeta=0.5*(y[i+1][j]-y[i-1][j])/del_zeta;
				y_eta=0.5*(y[i][j+1]-y[i][j-1])/del_eta;

				a[i][j]=y_eta*y_eta+x_eta*x_eta;
				b[i][j]=-(y_zeta*y_eta + x_zeta*x_eta);
				c[i][j]=x_zeta*x_zeta+y_zeta*y_zeta;
				d[i][j]=0.0;
				e[i][j]=0.0;
				J[i][j]=1.0/(x_zeta*y_eta-x_eta*y_zeta);

				xnew[i][j]= (a[i][j]/pow(del_zeta,2)-d[i][j]/(2.0*del_zeta))*x[i-1][j] + (a[i][j]/pow(del_zeta,2)+d[i][j]/(2.0*del_zeta))*x[i+1][j] + (c[i][j]/pow(del_eta,2)-e[i][j]/(2.0*del_eta))*x[i][j-1] + (c[i][j]/pow(del_eta,2)+e[i][j]/(2.0*del_eta))*x[i][j+1] + (0.5*b[i][j]/(del_zeta*del_eta))*(x[i+1][j+1]+x[i-1][j-1]-x[i-1][j+1]-x[i+1][j-1]);
				ynew[i][j]= (a[i][j]/pow(del_zeta,2)-d[i][j]/(2.0*del_zeta))*y[i-1][j] + (a[i][j]/pow(del_zeta,2)+d[i][j]/(2.0*del_zeta))*y[i+1][j] + (c[i][j]/pow(del_eta,2)-e[i][j]/(2.0*del_eta))*y[i][j-1] + (c[i][j]/pow(del_eta,2)+e[i][j]/(2.0*del_eta))*y[i][j+1] + (0.5*b[i][j]/(del_zeta*del_eta))*(y[i+1][j+1]+y[i-1][j-1]-y[i-1][j+1]-y[i+1][j-1]);

				xnew[i][j]=0.5*xnew[i][j]/( a[i][j]/pow(del_zeta,2) + c[i][j]/pow(del_eta,2) );
				ynew[i][j]=0.5*ynew[i][j]/( a[i][j]/pow(del_zeta,2) + c[i][j]/pow(del_eta,2) );

				if(fabs(xnew[i][j]-x[i][j])>maxerror)
					maxerror=fabs(xnew[i][j]-x[i][j]);
				x[i][j]=xnew[i][j];

				if(fabs(ynew[i][j]-y[i][j])>maxerror)
					maxerror=fabs(ynew[i][j]-y[i][j]);
				y[i][j]=ynew[i][j];
			}
			
	}
	
  }
else  // Interpolation along Y-direction
  {
     
    for(i=i_start+per*(Ny-1);i<i_end;i++)
    for(j=j_start+1;j<j_end;j++)
      {
   	j_s = j;
   	x[i][j] = x[i][j];
   	y[i][j] = y[i-1][j] + ( y[i_end][j] - y[(int) (i_start+per*(Ny-1)-1)][j] )/(i_end-i_start-per*(Ny-1)+1);	         
      }   
  
  }  
    
}

void gridgen_solid(int i_start,int i_end, int j_start,int j_end)
{

double 	a[Nx_s][Ny_s],b[Nx_s][Ny_s],c[Nx_s][Ny_s],d[Nx_s][Ny_s],e[Nx_s][Ny_s],xnew[Nx_s][Ny_s],ynew[Nx_s][Ny_s],J[Nx_s][Ny_s];
double	zeta_s_x_s,zeta_s_xx,zeta_s_y_s,zeta_s_yy,eta_s_x_s,eta_s_xx,eta_s_y_s,eta_s_yy,x_s_zeta_s,x_s_eta_s,y_s_zeta_s,y_s_eta_s;

	/* Finding x_s and y_s by solving laplacian equation */
	maxerror=1.0;
	while(maxerror>1e-8)
	{
		maxerror=0.0;
		for(i=i_start+1;i<i_end;i++)
			for(j=j_start+1;j<j_end;j++)
			{
				x_s_zeta_s=0.5*(x_s[i+1][j]-x_s[i-1][j])/del_zeta_s;
				x_s_eta_s=0.5*(x_s[i][j+1]-x_s[i][j-1])/del_eta_s;
				y_s_zeta_s=0.5*(y_s[i+1][j]-y_s[i-1][j])/del_zeta_s;
				y_s_eta_s=0.5*(y_s[i][j+1]-y_s[i][j-1])/del_eta_s;

				a[i][j]=y_s_eta_s*y_s_eta_s+x_s_eta_s*x_s_eta_s;
				b[i][j]=-(y_s_zeta_s*y_s_eta_s + x_s_zeta_s*x_s_eta_s);
				c[i][j]=x_s_zeta_s*x_s_zeta_s+y_s_zeta_s*y_s_zeta_s;
				d[i][j]=0.0;
				e[i][j]=0.0;
				J[i][j]=1.0/(x_s_zeta_s*y_s_eta_s-x_s_eta_s*y_s_zeta_s);

				xnew[i][j]= (a[i][j]/pow(del_zeta_s,2)-d[i][j]/(2.0*del_zeta_s))*x_s[i-1][j] + (a[i][j]/pow(del_zeta_s,2)+d[i][j]/(2.0*del_zeta_s))*x_s[i+1][j] + (c[i][j]/pow(del_eta_s,2)-e[i][j]/(2.0*del_eta_s))*x_s[i][j-1] + (c[i][j]/pow(del_eta_s,2)+e[i][j]/(2.0*del_eta_s))*x_s[i][j+1] + (0.5*b[i][j]/(del_zeta_s*del_eta_s))*(x_s[i+1][j+1]+x_s[i-1][j-1]-x_s[i-1][j+1]-x_s[i+1][j-1]);
				ynew[i][j]= (a[i][j]/pow(del_zeta_s,2)-d[i][j]/(2.0*del_zeta_s))*y_s[i-1][j] + (a[i][j]/pow(del_zeta_s,2)+d[i][j]/(2.0*del_zeta_s))*y_s[i+1][j] + (c[i][j]/pow(del_eta_s,2)-e[i][j]/(2.0*del_eta_s))*y_s[i][j-1] + (c[i][j]/pow(del_eta_s,2)+e[i][j]/(2.0*del_eta_s))*y_s[i][j+1] + (0.5*b[i][j]/(del_zeta_s*del_eta_s))*(y_s[i+1][j+1]+y_s[i-1][j-1]-y_s[i-1][j+1]-y_s[i+1][j-1]);

				xnew[i][j]=0.5*xnew[i][j]/( a[i][j]/pow(del_zeta_s,2) + c[i][j]/pow(del_eta_s,2) );
				ynew[i][j]=0.5*ynew[i][j]/( a[i][j]/pow(del_zeta_s,2) + c[i][j]/pow(del_eta_s,2) );

				if(fabs(xnew[i][j]-x_s[i][j])>maxerror)
					maxerror=fabs(xnew[i][j]-x_s[i][j]);
				x_s[i][j]=xnew[i][j];

				if(fabs(ynew[i][j]-y_s[i][j])>maxerror)
					maxerror=fabs(ynew[i][j]-y_s[i][j]);
				y_s[i][j]=ynew[i][j];
			}
	}
}

void mesh_vel()
 {
   for(i=1;i<Ny;i++)
   for(j=0;j<Nx;j++)
    {
      umesh[i][j]=(x[i][j]-x_old[i][j])/delt;
      vmesh[i][j]=(y[i][j]-y_old[i][j])/delt;
    }

}
