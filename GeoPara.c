#include "DesignAndFlowData.h"

void geometric_parameters()
{

/* calculation of dt and ds and xp */
for(j=1;j<=Nx-2;j++)
for(i=1;i<=Ny-2;i++)
{
  x_ne[i][j]=0.25*(x[i][j]+x[i+1][j]+x[i][j+1]+x[i+1][j+1]);
  y_ne[i][j]=0.25*(y[i][j]+y[i+1][j]+y[i][j+1]+y[i+1][j+1]);

  x_se[i][j]=0.25*(x[i][j]+x[i-1][j]+x[i][j+1]+x[i-1][j+1]);
  y_se[i][j]=0.25*(y[i][j]+y[i-1][j]+y[i][j+1]+y[i-1][j+1]);

  x_nw[i][j]=0.25*(x[i][j]+x[i+1][j]+x[i][j-1]+x[i+1][j-1]);
  y_nw[i][j]=0.25*(y[i][j]+y[i+1][j]+y[i][j-1]+y[i+1][j-1]);

  x_sw[i][j]=0.25*(x[i][j]+x[i-1][j]+x[i][j-1]+x[i-1][j-1]);
  y_sw[i][j]=0.25*(y[i][j]+y[i-1][j]+y[i][j-1]+y[i-1][j-1]);

  dt_e[i][j]=sqrt(pow(x_ne[i][j]-x_se[i][j],2) + pow(y_ne[i][j]-y_se[i][j],2));
  dt_w[i][j]=sqrt(pow(x_nw[i][j]-x_sw[i][j],2) + pow(y_nw[i][j]-y_sw[i][j],2));
  dt_n[i][j]=sqrt(pow(x_ne[i][j]-x_nw[i][j],2) + pow(y_ne[i][j]-y_nw[i][j],2));
  dt_s[i][j]=sqrt(pow(x_se[i][j]-x_sw[i][j],2) + pow(y_se[i][j]-y_sw[i][j],2));

  ds_e[i][j]=sqrt( pow(x[i][j+1]-x[i][j],2) + pow(y[i][j+1]-y[i][j],2) );
  ds_w[i][j]=sqrt( pow(x[i][j-1]-x[i][j],2) + pow(y[i][j-1]-y[i][j],2) );
  ds_n[i][j]=sqrt( pow(x[i+1][j]-x[i][j],2) + pow(y[i+1][j]-y[i][j],2) );
  ds_s[i][j]=sqrt( pow(x[i-1][j]-x[i][j],2) + pow(y[i-1][j]-y[i][j],2) );

  d1=sqrt( pow(x_ne[i][j]-x_sw[i][j],2) + pow(y_ne[i][j]-y_sw[i][j],2) );
  d2=sqrt( pow(x_nw[i][j]-x_se[i][j],2) + pow(y_nw[i][j]-y_se[i][j],2) );
  del_v[i][j]=0.25*sqrt( pow(2*d1*d2,2)-pow((dt_e[i][j]*dt_e[i][j]+dt_w[i][j]*dt_w[i][j]-dt_n[i][j]*dt_n[i][j]-dt_s[i][j]*dt_s[i][j]),2) );

  dels_xe[i][j]=y_ne[i][j]-y_se[i][j];
  dels_ye[i][j]=x_se[i][j]-x_ne[i][j];
  dels_yn[i][j]=x_ne[i][j]-x_nw[i][j];
  dels_xn[i][j]=-(y_ne[i][j]-y_nw[i][j]);
  dels_yw[i][j]=-(x_nw[i][j]-x_sw[i][j]);
  dels_xw[i][j]=y_nw[i][j]-y_sw[i][j];
  dels_xs[i][j]=y_sw[i][j]-y_se[i][j];
  dels_ys[i][j]=x_se[i][j]-x_sw[i][j];

  c_th_tx=(x_ne[i][j]-x_se[i][j])/sqrt( pow(x_ne[i][j]-x_se[i][j],2)+pow(y_ne[i][j]-y_se[i][j],2) );
  c_th_ty=(y_ne[i][j]-y_se[i][j])/sqrt( pow(x_ne[i][j]-x_se[i][j],2)+pow(y_ne[i][j]-y_se[i][j],2) );
  c_th_sx=(x[i][j+1]-x[i][j])/sqrt( pow(x[i][j+1]-x[i][j],2)+pow(y[i][j+1]-y[i][j],2) );
  c_th_sy=(y[i][j+1]-y[i][j])/sqrt( pow(x[i][j+1]-x[i][j],2)+pow(y[i][j+1]-y[i][j],2) );

  dels_se[i][j]=(dels_xe[i][j]*c_th_ty-dels_ye[i][j]*c_th_tx)/(c_th_sx*c_th_ty-c_th_sy*c_th_tx);
  dels_te[i][j]=-(dels_xe[i][j]*c_th_sy-dels_ye[i][j]*c_th_sx)/(c_th_sx*c_th_ty-c_th_sy*c_th_tx);

  c_th_tx=(x_ne[i][j]-x_nw[i][j])/sqrt( pow(x_ne[i][j]-x_nw[i][j],2)+pow(y_ne[i][j]-y_nw[i][j],2) );
  c_th_ty=(y_ne[i][j]-y_nw[i][j])/sqrt( pow(x_ne[i][j]-x_nw[i][j],2)+pow(y_ne[i][j]-y_nw[i][j],2) );
  c_th_sx=(x[i+1][j]-x[i][j])/sqrt( pow(x[i+1][j]-x[i][j],2)+pow(y[i+1][j]-y[i][j],2) );
  c_th_sy=(y[i+1][j]-y[i][j])/sqrt( pow(x[i+1][j]-x[i][j],2)+pow(y[i+1][j]-y[i][j],2) );

  dels_sn[i][j]=(dels_xn[i][j]*c_th_ty-dels_yn[i][j]*c_th_tx)/(c_th_sx*c_th_ty-c_th_sy*c_th_tx);
  dels_tn[i][j]=-(dels_xn[i][j]*c_th_sy-dels_yn[i][j]*c_th_sx)/(c_th_sx*c_th_ty-c_th_sy*c_th_tx);

  if(j>1)
   {
    dels_sw[i][j]=dels_se[i][j-1];
    dels_tw[i][j]=dels_te[i][j-1];
   }

  if(i>1)
   {
    dels_ss[i][j]=dels_sn[i-1][j];
    dels_ts[i][j]=dels_tn[i-1][j];
   }

  if(j==1)
   {
    dels_sw[i][j]=dels_se[i][j];
    dels_tw[i][j]=dels_te[i][j];
   }

  if(i==1)
   {
    dels_ss[i][j]=dels_sn[i][j];
    dels_ts[i][j]=dels_tn[i][j];
   }

  dels_e[i][j]=sqrt( pow(dels_se[i][j],2) + pow(dels_te[i][j],2) );
  dels_w[i][j]=sqrt( pow(dels_sw[i][j],2) + pow(dels_tw[i][j],2) );
  dels_n[i][j]=sqrt( pow(dels_sn[i][j],2) + pow(dels_tn[i][j],2) );
  dels_s[i][j]=sqrt( pow(dels_ss[i][j],2) + pow(dels_ts[i][j],2) );

}

for(i=0;i<Ny;i++)
 {
    del_v[i][0]=del_v[i][1];
    del_v[i][Nx-1]=del_v[i][Nx-2];
 }
 
for(j=0;j<Nx;j++)
 {
    del_v[0][j]=del_v[1][j];
    del_v[Ny-1][j]=del_v[Ny-2][j];
 }

}
