/* FIPP and Fully-Implicit Coupling-based in-house FSI solver with Novel Hybrid AMM-MRR-SEMM Grid Genration 
   Methodology on a Curvilinear Structured Grid with Fully FVM-based discretisation */
   
/* Authors:
   Sumant R Morab, Dr. Atul Sharma, and Dr. Janani S Murallidharan,
   Department of Mechanical Engineering, Indian Institute of Technology Bombay, Mumbai, India */      
   
#include "DesignAndFlowData.h"

void main()
{

char* name 	= "re500_t";
char* nameWSS 	= "re500_WSS_t";
char* extension = ".dat";
char fileSpec[strlen(name)+strlen(extension)+5];
char fileSpec2[strlen(nameWSS)+strlen(extension)+5];

/* Define Flow Paramaters */
define_flow_para();

/* initial conditions -- fluid & solid */
initial_cond_fluid();
initial_cond_solid();

it	= 0;
time	= 0.0;
it2	= 0;
ct	= 0.1/delt; // interval for data storage

/* Enter the Time-Marching */
while(time<100.0)
{

time	= time+delt;
maxtol	= 1.0;
iter	= 0;

if(time>0)
   {   
      gridgen(0,Ny-1,0,Nx-1);        // novel fluid grid generation
      mesh_vel();		     // Calculate mesh velocity
      geometric_parameters();        // fluid geometric parameters update     
   }

bc_fluid();
maxerror = 0.0;

// Solve flow equations 
solve_fluid();

for(j=0;j<Nx;j++)
for(i=0;i<Ny;i++)
   {
     u_b[i][j]=u[i][j];
     v_b[i][j]=v[i][j];
   } 

// Solve solid equations
if(time>0)
solve_solid();

// Print the residuals
printf("time = %lf\t iterations = %d\n",time,iter);
it2++;

/* writing data to file */
if(it2%ct==0)
    {
	FILE *fp2;
	snprintf( fileSpec, sizeof( fileSpec ), "%s%d%s", name,it2/ct,extension );
	fp2 = fopen(fileSpec,"w");
	fprintf(fp2,"Variables=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
	fprintf(fp2,"ZONE F=POINT\n");
	fprintf(fp2,"I=%d,J=%d\n",Nx,Ny);

	for(i=0;i<Ny;i++)
	for(j=0;j<Nx;j++)	
  	    fprintf( fp2,"%lf\t%lf\t%lf\t%lf\t%lf\n",x[i][j],y[i][j],u[i][j],v[i][j],p_old[i][j]);
	fclose(fp2);
    }

} 

}
