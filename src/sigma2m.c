#include <R.h>
#include <Rmath.h>
#include <float.h>

void sigma2m(double *Nxk, double *m, int *n, int *nx, double *Yxsort, double *sigma2)
{
 long i, j, k, indY=0; 
 double aa, res, bb=0.0, cc=0.0;
   for(k=0; k < *nx; k++) {
    res=0.0;
     if(Nxk[k]==1){
     sigma2[k]=0;
     }
     else{
       for(i=0; i < Nxk[k]-1; i++) {
        cc=0.5*(R_pow(((double) i+1.0)/(Nxk[k]), (m[k]*2)-1))*(1.0-(((double) i+1.0)/(Nxk[k])))*(R_pow(Yxsort[indY+i+1]-Yxsort[indY+i],2));
         if(i>0){
          aa=0.0;
           for(j=0; j < i; j++){
            aa = aa + (R_pow(((double) j+1.0)/(Nxk[k]),m[k]))*(Yxsort[indY+j+1]-Yxsort[indY+j]);
            }                 
           bb=aa*(R_pow(((double) i+1.0)/(Nxk[k]),m[k]-1))*(1.0-((double) i+1.0)/(Nxk[k]))*(Yxsort[indY+i+1]-Yxsort[indY+i]);
           }
          else{
           bb=0.0;
           }
         res = res + cc + bb;
         }
      sigma2[k]=2.0*R_pow(m[k],2)**n*res/(Nxk[k]);
      indY=indY+Nxk[k];
      } 
    }
}
