//####################################################################
//ICA test 
//####################################################################
//2 May 2004
//Contact: kraskov@its.caltech.edu
//####################################################################
//uses mic_xnyn, mir_xnyn
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "miutils.h"

int main(int argc, char **argv) {

  FILE *fin;
  int i,d,d1,d2,N_,ddel;
  int N=4096;;
  int dim=3;
  int K=1;
  
  double **x;
  double **xx,**yy;
  double *psi,*min,*max,*scalxx;
  double t_d,t_d1,st_d,ct_d;
  int BOX1;
  double **mi;

  double s,me;

  int k;
  double *rmi;
  double **varmi;
  double angle;

  int tau=1;
  int edim=2;
  int method=1;
  double addnoise=-1;
  int nr=4;

  if (argc<5) {
    fprintf(stderr,"\nReliability Test using MI based on k-nearest neighbours statistics (rectangular)\n\n");
    fprintf(stderr,"Usage:\n%s <filename> <dim> <# points> <# neighbours> [tau] [edim] [method] [nrot] [addnoise]\n\n",argv[0]);
    fprintf(stderr,"Input:\n\t<filename>\ttext file with <dim> columns and <# points> rows\n");
    fprintf(stderr,"\t<dim>\t\tnumber of columns in file\n");
    fprintf(stderr,"\t<# points>\tnumber of rows (length of characteristic vector)\n");
    fprintf(stderr,"\t<# neighbours>\tnumber of the nearest neighbours for MI estimator\n");
    fprintf(stderr,"\t[tau]\t\tembedding delay, default 1\n");
    fprintf(stderr,"\t[edim]\t\tembedding dimension, default 2\n");
    fprintf(stderr,"\t[method]\teither cubic (0), or rectangular (1); default 1\n");
    fprintf(stderr,"\t[nrot]\t\tnumber of angles to check; default 4 (pi/8,pi/4,3pi/8,pi/2)\n");
    fprintf(stderr,"\t[addnoise]\tnoise amplitude; default 1e-8\n");
    fprintf(stderr,"\nOutput:\n");
    fprintf(stderr,"\nMI matrix and variance of MI matrix\n");
    fprintf(stderr,"\nContact: kraskov@its.caltech.edu\n");
    exit(-1);
  }

  dim=atoi(argv[2]);
  N=atoi(argv[3]);
  K=atoi(argv[4]);
  
  if (argc>=6) tau=atoi(argv[5]);
  if (argc>=7) edim=atoi(argv[6]);
  if (argc>=8) method=atoi(argv[7]);
  if (argc>=9) nr=atoi(argv[8]);
  if (argc>=10) addnoise=atof(argv[9]);
  if (argc>=11) {fprintf(stderr,"Too many input arguments\n");exit(-1);}

  yy=(double **)calloc(2,sizeof(double*));
  yy[0]=(double *)calloc(N,sizeof(double));
  yy[1]=(double *)calloc(N,sizeof(double));

  xx=(double **)calloc(edim*2,sizeof(double*));
  for (ddel=0;ddel<edim;ddel++) {
    xx[ddel]=(double *)calloc(N,sizeof(double));
    xx[ddel+edim]=(double *)calloc(N,sizeof(double));
  }  

  x=(double **)calloc(dim,sizeof(double*));
  
  for (d=0;d<dim;d++) {
    x[d]=(double *)calloc(N,sizeof(double));
  }

  //reading of the data
  fin=fopen(argv[1],"r");
  if (fin) 
    for (i=0;i<N;i++) {
      for (d=0;d<dim;d++) {
	fscanf(fin,"%lf",&(x[d][i]));
      }
    }
  else { fprintf(stderr,"File %s doesn't exist\n",argv[1]);exit(-1);}
  fclose(fin);
  

  // add noise
  if (addnoise) {
    srand((dim+edim+tau)*N*K*int(x[(dim)/2][N/10]));
    if (addnoise==-1) for (d=0;d<dim;d++) for (i=0;i<N;i++) x[d][i]+=(1.0*rand()/RAND_MAX)*1e-8;
    else for (d=0;d<dim;d++) for (i=0;i<N;i++) x[d][i]+=(1.0*rand()/RAND_MAX)*addnoise;
  }

  min=(double*)calloc(2*edim,sizeof(double));
  max=(double*)calloc(2*edim,sizeof(double));

  psi=(double*)calloc(N+1,sizeof(double));
  psi[1]=-(double).57721566490153;
  for (i=1;i<N;i++) psi[i+1]=psi[i]+1/(double)i; //cubic
  BOX1=N-5;
  scalxx=(double*)calloc(2*edim,sizeof(double));

  rmi=(double *)calloc(nr+1,sizeof(double));
 
  mi=(double **)calloc(dim,sizeof(double *));
  varmi=(double **)calloc(dim,sizeof(double *));
  for (d=0;d<dim;d++) {
    mi[d]=(double*)calloc(dim,sizeof(double));
    varmi[d]=(double*)calloc(dim,sizeof(double));
  }
  
  N_=N-(edim-1)*tau;
  
  for (d=0;d<dim;d++) {
    for (d2=d+1;d2<dim;d2++) {
      for (k=0;k<nr;k++) {
	angle=1.0*k/nr*M_PI/2;
	st_d=sin(angle);ct_d=cos(angle);
	for (i=0;i<N_;i++) {
	  t_d=x[d][i];t_d1=x[d2][i];
	  yy[0][i]=ct_d*t_d+st_d*t_d1;
	  yy[1][i]=ct_d*t_d1-st_d*t_d;
	}
	for (ddel=0;ddel<edim;ddel++) {
	  memcpy(xx[ddel],yy[0]+ddel*tau,N_*sizeof(double));
	  memcpy(xx[ddel+edim],yy[1]+ddel*tau,N_*sizeof(double));
	}
	
	//shift to positive values
	for (d1=0;d1<2*edim;d1++) {min[d1]=DBL_MAX/2;max[d1]=-DBL_MAX/2;}
	for (d1=0;d1<2*edim;d1++) {
	  for (i=0;i<N_;i++) {
	    if (xx[d1][i]<min[d1]) min[d1]=xx[d1][i]; 
	    if (xx[d1][i]>max[d1]) max[d1]=xx[d1][i];
	  }
	  for (i=0;i<N_;i++) xx[d1][i]=xx[d1][i]-min[d1];
	  scalxx[d1]=BOX1/(max[d1]-min[d1]);
	}
	switch (method) {
	case 0 :  mic_xnyn(xx,edim,edim,N_,K,psi,scalxx,&(rmi[k])); break;
	case 1 :  mir_xnyn(xx,edim,edim,N_,K,psi,scalxx,&(rmi[k])); break;
	}
	fprintf(stdout,"%f\t",rmi[k]);
      }
      fprintf(stdout,"\n");
      
      mi[d][d2]=rmi[0];
      mi[d2][d]=mi[d][d2];
      
      me=s=0; for (k=0;k<=nr;k++) me+=rmi[k];
      me/=(nr+1);  
      varmi[d][d2]=me-rmi[0];
      varmi[d2][d]=varmi[d][d2];
    }
  }
  
  fprintf(stdout,"\n");

  for (d=0;d<dim;d++) {
    for (d2=0;d2<dim;d2++) {
      if (d==d2) {fprintf(stdout," 0\t\t");continue;}
      if (mi[d][d2]>0) fprintf(stdout," ");
      fprintf(stdout,"%2.8f\t",mi[d][d2]);
    }
    fprintf(stdout,"\n");
  }
  
  fprintf(stdout,"\n");
    
  for (d=0;d<dim;d++) {
    for (d2=0;d2<dim;d2++) {
      if (d==d2) {fprintf(stdout," 0\t\t");continue;}
      if (varmi[d][d2]>0) fprintf(stdout," ");
      fprintf(stdout,"%2.8f\t",varmi[d][d2]);
    }
    fprintf(stdout,"\n");
  }

  free(yy[0]);free(yy[1]);free(yy);

  for (ddel=0;ddel<2*edim;ddel++) free(xx[ddel]);
  free(xx);
  for (d=0;d<dim;d++) free(x[d]);
  free(x);

  free(min);
  free(max);
  
  free(psi);
  free(scalxx);
  free(rmi);

  for (d=0;d<dim;d++) {
    free(mi[d]);free(varmi[d]);
  }
  free(mi);free(varmi);
}


