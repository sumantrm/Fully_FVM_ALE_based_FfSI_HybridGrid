#include "DesignAndFlowData.h"

void initial_cond_solid()
{

	del_eta_s=H/(Ny_s-1);
	del_zeta_s=Le/(Nx_s-1);

	mu_s=E/(2.0*(1.0+nu));
	lamda=E/((1.0+nu)*(1.0-nu));

	it=Nx_s;
	it3=0;

	printf("******************************\n");
	for(j=0;j<Nx_s;j++)
		for(i=0;i<Ny_s;i++)
		{
			eta_s[i][j]  = i*del_eta_s;
			zeta_s[i][j] = j*del_zeta_s;

			x_s[i][j]    = zeta_s[i][j];
			y_s[i][j]    = eta_s[i][j] - H ;
			
			x_s_old[i][j]    = zeta_s[i][j];
			y_s_old[i][j]    = eta_s[i][j] - H ;
			
			x_init[i][j]    = zeta_s[i][j];
			y_init[i][j]    = eta_s[i][j] - H ;

			ds_s_e[i][j] = del_zeta_s;
			ds_s_w[i][j] = del_zeta_s;
			ds_s_n[i][j] = del_eta_s;
			ds_s_s[i][j] = del_eta_s;
			dt_s_e[i][j] = del_eta_s;
			dt_s_w[i][j] = del_eta_s;
			dt_s_n[i][j] = del_zeta_s;
			dt_s_s[i][j] = del_zeta_s;

		}

         //for(j=0;j<Nx;j++)
	//printf("%lf\n",x_s[Ny_s-1][j]);

         gridgen_solid(0,Ny_s-1,0,Nx_s-1);

/* Setting up Boundary x and y value by reading data from csv files */
for(j=1;j<=Nx_s-2;j++)
for(i=1;i<=Ny_s-2;i++)
{

	x_s_ne[i][j]=0.25*(x_s[i][j]+x_s[i][j+1]+x_s[i+1][j]+x_s[i+1][j+1]);
	y_s_ne[i][j]=0.25*(y_s[i][j]+y_s[i][j+1]+y_s[i+1][j]+y_s[i+1][j+1]);
	x_s_se[i][j]=0.25*(x_s[i][j]+x_s[i][j+1]+x_s[i-1][j]+x_s[i-1][j+1]);
	y_s_se[i][j]=0.25*(y_s[i][j]+y_s[i][j+1]+y_s[i-1][j]+y_s[i-1][j+1]);

	x_s_nw[i][j]=0.25*(x_s[i][j]+x_s[i][j-1]+x_s[i+1][j]+x_s[i+1][j-1]);
	y_s_nw[i][j]=0.25*(y_s[i][j]+y_s[i][j-1]+y_s[i+1][j]+y_s[i+1][j-1]);
	x_s_sw[i][j]=0.25*(x_s[i][j]+x_s[i][j-1]+x_s[i-1][j]+x_s[i-1][j-1]);
	y_s_sw[i][j]=0.25*(y_s[i][j]+y_s[i][j-1]+y_s[i-1][j]+y_s[i-1][j-1]);

	dt_s_e[i][j]=sqrt(pow(x_s_ne[i][j]-x_s_se[i][j],2) + pow(y_s_ne[i][j]-y_s_se[i][j],2));
	dt_s_w[i][j]=sqrt(pow(x_s_nw[i][j]-x_s_sw[i][j],2) + pow(y_s_nw[i][j]-y_s_sw[i][j],2));
	dt_s_n[i][j]=sqrt(pow(x_s_ne[i][j]-x_s_nw[i][j],2) + pow(y_s_ne[i][j]-y_s_nw[i][j],2));
	dt_s_s[i][j]=sqrt(pow(x_s_se[i][j]-x_s_sw[i][j],2) + pow(y_s_se[i][j]-y_s_sw[i][j],2));

	ds_s_e[i][j]=sqrt( pow(x_s[i][j+1]-x_s[i][j],2) + pow(y_s[i][j+1]-y_s[i][j],2) );
	ds_s_w[i][j]=sqrt( pow(x_s[i][j-1]-x_s[i][j],2) + pow(y_s[i][j-1]-y_s[i][j],2) );
	ds_s_n[i][j]=sqrt( pow(x_s[i+1][j]-x_s[i][j],2) + pow(y_s[i+1][j]-y_s[i][j],2) );
	ds_s_s[i][j]=sqrt( pow(x_s[i-1][j]-x_s[i][j],2) + pow(y_s[i-1][j]-y_s[i][j],2) );

	d1=sqrt( pow(x_s_ne[i][j]-x_s_sw[i][j],2) + pow(y_s_ne[i][j]-y_s_sw[i][j],2) );
	d2=sqrt( pow(x_s_nw[i][j]-x_s_se[i][j],2) + pow(y_s_nw[i][j]-y_s_se[i][j],2) );
	del_v_s[i][j]=0.25*sqrt( pow(2*d1*d2,2)-pow((dt_s_e[i][j]*dt_s_e[i][j]+dt_s_w[i][j]*dt_s_w[i][j]-dt_s_n[i][j]*dt_s_n[i][j]-dt_s_s[i][j]*dt_s_s[i][j]),2) );

	alpha_e[i][j]=asin((x_s_se[i][j]-x_s_ne[i][j])/dt_s_e[i][j]);
	beta_e[i][j]=asin((y_s[i][j+1]-y_s[i][j])/ds_s_e[i][j]);

	alpha_n[i][j]=PI/2.0-asin((y_s_nw[i][j]-y_s_ne[i][j])/dt_s_n[i][j]);
	beta_n[i][j]=acos((x_s[i+1][j]-x_s[i][j])/ds_s_n[i][j]);

	alpha_w[i][j]=PI+asin((x_s_sw[i][j]-x_s_nw[i][j])/dt_s_w[i][j]);
	beta_w[i][j]=PI-asin((y_s[i][j-1]-y_s[i][j])/ds_s_w[i][j]);

	if(i>1)
	{
		alpha_s[i][j]=PI+alpha_n[i-1][j];
		beta_s[i][j]=PI+beta_n[i-1][j];
	}

	if(i==1)
	{
		alpha_s[i][j]=PI+alpha_n[i][j];
		beta_s[i][j]=PI+beta_n[i][j];
	}

}

	
	/* Initial Conditions */
	for(j=0;j<Nx_s;j++)
	for(i=0;i<Ny_s;i++)
		{
			u_s[i][j]=0.0;
			v_s[i][j]=0.0;

			u_s_old[i][j]=0.0;
			v_s_old[i][j]=0.0;

			u_s_new[i][j]=0.0;
			v_s_new[i][j]=0.0;

			x_init[i][j]=x_s[i][j];
			y_init[i][j]=y_s[i][j];

			del_v_s[i][j] = del_eta_s*del_zeta_s;

		}

	FILE *fp2;
	fp2 = fopen("initial_solid.dat","w");
	fprintf(fp2,"Variables=\"X\",\"Y\",\"U\",\"V\",\"Alpha_N\",\"Beta_N\"\n");
	fprintf(fp2,"ZONE F=POINT\n");
	fprintf(fp2,"I=%d,J=%d\n",Nx_s,Ny_s);

	for(i=0;i<Ny_s;i++)
	for(j=0;j<Nx_s;j++)	
  	    fprintf( fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x_s[i][j],y_s[i][j],u_s[i][j],v_s[i][j],alpha_n[i][j],beta_n[i][j]);
	fclose(fp2);


}
