#include "DesignAndFlowData.h"

void initial_cond_fluid()
{

del_eta=D/(Ny-1.0);
del_zeta=L/(Nx-1.0);

printf("******************************\n");

for(i=0;i<Ny;i++)
for(j=0;j<Nx;j++)
{

  eta[i][j]=i*del_eta;
  zeta[i][j]=j*del_zeta;

  x[i][j]=zeta[i][j];
  y[i][j]=eta[i][j];
  
  if(i==0)
  x_oldo[j] = x[i][j];

  ds_e[i][j]=del_zeta;
  ds_w[i][j]=del_zeta;
  ds_n[i][j]=del_eta;
  ds_s[i][j]=del_eta;
  dt_e[i][j]=del_eta;
  dt_w[i][j]=del_eta;
  dt_n[i][j]=del_zeta;
  dt_s[i][j]=del_zeta;
 
}

//gridgen(0,Ny-1,0,Nx-1);
geometric_parameters();

/* initial Conditions */
for(i=0;i<Ny;i++)
for(j=0;j<Nx;j++)
{
 u[i][j]=0.0;
 v[i][j]=0.0;
 p[i][j]=0.0;
 u_b[i][j]=0.0;
 v_b[i][j]=0.0;
 p_old[i][j]=0.0;
 p_new[i][j]=0.0;

 if(j==0)
    del_v[i][0]=del_v[i][1];

 if(j==Nx-1)
    del_v[i][Nx-1]=del_v[i][Nx-2];

 if(i==0)
    del_v[0][j]=del_v[1][j];

 if(i==Ny-1)
    del_v[Ny-1][j]=del_v[Ny-2][j];
  
}

del_v[0][0]=del_v[0][1];
del_v[0][Nx-1]=del_v[0][Nx-2];

	/* writing to fjle */
	FILE *fp2;
	fp2 = fopen("init_fluid.dat","w");
	fprintf(fp2,"Variables=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
	fprintf(fp2,"ZONE F=POINT\n");
	fprintf(fp2,"I=%d,J=%d\n",Nx,Ny);

	for(i=0;i<Ny;i++)
	for(j=0;j<Nx;j++)	
  	    fprintf( fp2,"%lf\t%lf\t%lf\t%lf\t%lf\n",x[i][j],y[i][j],u[i][j],v[i][j],p[i][j]);
	fclose(fp2);

}
