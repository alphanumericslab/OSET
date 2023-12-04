//####################################################################
//Mutual information estimator, rectangular version 
//####################################################################
//2 May 2004
//Contact: kraskov@its.caltech.edu
//####################################################################
//uses mir_xnyn
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "miutils.h"

int main(int argc, char **argv) {

  FILE *fin;
  int i;
  double **x;
  double *scal;
  double *min;
  double *max;
  double *psi;
  int K,N;
  int d;
  double mir;
  int dimx,dimy;


  int BOX1;

  double s,me;
  double addnoise=-1;

  if (argc<6) {
    fprintf(stderr,"\nMutual Infomation (MI) k-nearest neighbours statistics (rectangular)\n\n");
    fprintf(stderr,"Usage:\n%s <filename> <dimx> <dimy> <# points> <# neighbours> [addnoise]\n\n",argv[0]);
    fprintf(stderr,"Input:\n\t<filename>\ttext file with <dimx+dimy> columns and <# points> rows\n");
    fprintf(stderr,"\t<dimx>\t\tnumber of columns for X (dimension of X)\n");
    fprintf(stderr,"\t<dimy>\t\tnumber of columns for Y (dimension of Y)\n");
    fprintf(stderr,"\t<# points>\tnumber of rows (length of characteristic vector)\n");
    fprintf(stderr,"\t<# neighbours>\tnumber of the nearest neighbours for MI estimator\n");
    fprintf(stderr,"\t[addnoise]\tnoise amplitude; default 1e-8\n");
    fprintf(stderr,"\nOutput:\n");
    fprintf(stderr,"\nMI\n");
    fprintf(stderr,"\nContact: kraskov@its.caltech.edu\n");
    exit(-1);
  }



  dimx=atoi(argv[2]);
  dimy=atoi(argv[3]);
  N=atoi(argv[4]);
  K=atoi(argv[5]);
  if (argc==7) {addnoise=atof(argv[6]);}
  if (argc>=8) {fprintf(stderr,"Too many input arguments\n");exit(-1);}

  x=(double**)calloc(dimx+dimy,sizeof(double*));
  for (d=0;d<dimx+dimy;d++) x[d]=(double*)calloc(N,sizeof(double));
  scal=(double*)calloc(dimx+dimy,sizeof(double));
  min=(double*)calloc(dimx+dimy,sizeof(double));
  max=(double*)calloc(dimx+dimy,sizeof(double));
  for (d=0;d<dimx+dimy;d++) {min[d]=DBL_MAX/2;max[d]=-DBL_MAX/2;}
  //reading of the data
  fin=fopen(argv[1],"r");
  if (fin)
    for (i=0;i<N;i++) {
      for (d=0;d<dimx+dimy;d++) {
	fscanf(fin,"%lf",&(x[d][i]));
      }
    }
  else { fprintf(stderr,"File %s doesn't exist\n",argv[1]);exit(-1);}
  fclose(fin);
  // add noise
  if (addnoise) {
    srand((dimx+dimy)*N*K*int(x[(dimx+dimy)/2][N/10]));
    if (addnoise==-1) for (d=0;d<dimx+dimy;d++) for (i=0;i<N;i++) x[d][i]+=(1.0*rand()/RAND_MAX)*1e-8;
    else for (d=0;d<dimx+dimy;d++) for (i=0;i<N;i++) x[d][i]+=(1.0*rand()/RAND_MAX)*addnoise;
  }

  //normalization
  for (d=0;d<dimx+dimy;d++) {
    me=s=0; for (i=0;i<N;i++) me+=x[d][i];
    me/=N;  for (i=0;i<N;i++) s+=(x[d][i]-me)*(x[d][i]-me);
    s/=(N-1);s=sqrt(s);
    if (s==0) {;}
    for (i=0;i<N;i++) {
      x[d][i] = (x[d][i]-me)/s;
      if (x[d][i]<min[d]) min[d]=x[d][i]; 
      if (x[d][i]>max[d]) max[d]=x[d][i];
    }
    for (i=0;i<N;i++) x[d][i]=x[d][i]-min[d];
  }

  psi=(double*)calloc(N+1,sizeof(double)); 
  psi[1]=-(double).57721566490153;
  for (i=1;i<N;i++) psi[i+1]=psi[i]+1/(double)i;
  BOX1=N-5;
  for (d=0;d<dimx+dimy;d++) scal[d]=BOX1/(max[d]-min[d]); 

  mir_xnyn(x,dimx,dimy,N,K,psi,scal,&mir);
  fprintf(stdout,"%1.8f\n",mir);

  for (d=0;d<dimx+dimy;d++) free(x[d]); free(x);
  free(scal);
  free(min);free(max);
  free(psi);
}
