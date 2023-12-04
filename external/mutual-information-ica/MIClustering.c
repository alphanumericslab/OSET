//####################################################################
//Mutual Information Clustering, rectangular version 
//####################################################################
//2 May 2004
//Contact: kraskov@its.caltech.edu
//####################################################################
//uses mi2r, mir_xnyn, redr
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "miutils.h"

#define MIM 0
void maxx(double **mi, char *is, int dim, int *im, int *jm, double *mm);

int main(int argc, char **argv) {

  FILE *fin;
  int i,d,d2,im,jm,d1;
  int dx,dy,dmxy;
  int N;
  int dim;
  int K=10;
  double mm;
  char *is;
  int **noden;

  double **mis; // mi/(dx+dy) matrix
  double *clmis; // mi((X),(Y))/(dx+dy) of the combination

  double **mi; // mi matrix
  double *clmi; // mi((X),(Y)) of the combination

  int *cls[2];

  
  double **x;
  double **xx;
  double *psi,*min,*max,*scal,*scalxx;
  double t_d;
  double *red; //redundancy

  //read arguments
  double s,me;
  double addnoise=-1;
  int pdm=0;
  
  if (argc<5) {//first three arguments are necessary
    fprintf(stderr,"\nHierarchical Clustering Based on Mutual Infomation (MI)\n\n");
    fprintf(stderr,"Usage:\n%s <filename> <dim> <# points> <# neighbours> [addnoise] [P_MI]\n\n",argv[0]);
    fprintf(stderr,"Input:\n\t<filename>\ttext file with <dim> columns and <# points> rows\n");
    fprintf(stderr,"\t<dim>\t\tnumber of columns (number of characteristic vectors)\n");
    fprintf(stderr,"\t<# points>\tnumber of rows (length of characteristic vectors)\n");
    fprintf(stderr,"\t<# neighbours>\tnumber of the nearest neighbours for MI estimator\n");
    fprintf(stderr,"\t[addnoise]\tnoise amplitude; default 1e-8\n");
    fprintf(stderr,"\t[P_MI]\t\tprint_pairwise_MI; default do not print; to print say 1\n");
    fprintf(stderr,"\nOutput: (suitable for Matlab function dendrogram.m, just use [1 2 5])\n");
    fprintf(stderr,"\nnodeX nodeY MI(X,Y)/(dx+dy) MI(X,Y) RED(X,Y)\n");
    fprintf(stderr,"\tMI(X,Y)/(dx+dy)\tinter-cluster similarity measure for clusters\n\t\t\t\tnodeX and nodeY (used for topology)\n");
    fprintf(stderr,"\tMI(X,Y)\t\tmutual information between nodeX and nodeY\n");
    fprintf(stderr,"\tRED(X,Y)\tmutual information of the cluster (nodeXnodeY)\n\t\t\t\t(used for plotting a dendrogram)\n");
    fprintf(stderr,"\nContact: kraskov@its.caltech.edu\n");
    exit(-1);
  }
  dim=atoi(argv[2]);
  N=atoi(argv[3]);
  K=atoi(argv[4]);
  if (argc>=6) addnoise=atof(argv[5]);
  if (argc>=7) pdm=atoi(argv[6]);
  if (argc>=8) {fprintf(stderr,"Too many input arguments\n");exit(-1);}

  x=(double **)calloc(dim,sizeof(double*));
  xx=(double **)calloc(dim,sizeof(double*));
  for (d=0;d<dim;d++) {
    x[d]=(double *)calloc(N,sizeof(double));
    xx[d]=(double *)calloc(N,sizeof(double));
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
    srand(dim*N*K*int(x[dim/2][N/10])); //just random number
    if (addnoise==-1) for (d=0;d<dim;d++) for (i=0;i<N;i++) x[d][i]+=(1.0*rand()/RAND_MAX)*1e-8;
    else for (d=0;d<dim;d++) for (i=0;i<N;i++) x[d][i]+=(1.0*rand()/RAND_MAX)*addnoise;
  }

  min=(double*)calloc(dim,sizeof(double));
  max=(double*)calloc(dim,sizeof(double));
  for (d=0;d<dim;d++) {min[d]=DBL_MAX/2;max[d]=-DBL_MAX/2;}

  psi=(double*)calloc(N+1,sizeof(double)); 
  psi[1]=-(double).57721566490153;
  for (i=1;i<N;i++) psi[i+1]=psi[i]+1/(double)i; //cubic
  scal=(double*)calloc(dim,sizeof(double));
  scalxx=(double*)calloc(dim,sizeof(double));

  //normalization
  for (d=0;d<dim;d++) {
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

  is=(char *)calloc(2*dim-1,sizeof(char));

  noden=(int **)calloc(2*dim-1,sizeof(int *));
  mi=(double **)calloc(2*dim-1,sizeof(double *));
  mis=(double **)calloc(2*dim-1,sizeof(double *));
  for (d=0;d<2*dim-1;d++) {
    mi[d]=(double*)calloc(2*dim-1,sizeof(double));
    mis[d]=(double*)calloc(2*dim-1,sizeof(double));
    noden[d]=(int *)calloc(dim+1,sizeof(int));
  }
  cls[0]=(int*)calloc(dim-1,sizeof(int));
  cls[1]=(int*)calloc(dim-1,sizeof(int));
  clmi=(double*)calloc(dim-1,sizeof(double));
  clmis=(double*)calloc(dim-1,sizeof(double));
  red=(double*)calloc(dim-1,sizeof(double));
  for (d=0;d<dim-1;d++) {
    red[d]=clmi[d]=clmis[d]==0.0;
    cls[0][d]=cls[1][d]=0;
  }

  for (d=0;d<2*dim-1;d++) {
    for (d2=0;d2<2*dim-1;d2++) {
      mi[d][d2]=MIM;
      mis[d][d2]=MIM;
    }
  }
  
  for (d=0;d<dim;d++) {
    for (d2=d+1;d2<dim;d2++) {
      memcpy(xx[0],x[d],N*sizeof(double));
      memcpy(xx[1],x[d2],N*sizeof(double));
      scalxx[0]=scal[d];
      scalxx[1]=scal[d2];

      mi2r(xx,N,K,psi,scalxx,&(t_d));

      if (t_d<0) mi[d][d2]=MIM;
      else mi[d][d2]=t_d;
      mi[d2][d]=mi[d][d2];
      mis[d][d2]=mi[d][d2]/2;
      mis[d2][d]=mis[d][d2];
    }
  }

  if (pdm) {
    fprintf(stderr,"%d\n",dim);
    for (dx=0;dx<dim;dx++) {
      if (dx+1>=100) fprintf(stderr,"f%d",dx+1);
      else if (dx+1>=10) fprintf(stderr,"f0%d",dx+1);
      else fprintf(stderr,"f00%d",dx+1);
      
      for (dy=0;dy<dx;dy++) {
	fprintf(stderr,"\t%1.3f",mi[dy][dx]*(mi[dy][dx]>0));
      }
      fprintf(stderr,"\t0;\n");
    }
  }

  for (d=0;d<dim;d++) {
    is[d]=1;
    noden[d][0]=d;
    for (d2=1;d2<=dim;d2++) noden[d][d2]=-1;
  }
  for (d=dim;d<2*dim-1;d++) {
    is[d]=0;
    for (d2=0;d2<=dim;d2++) noden[d][d2]=-1;
  }
  
  
  for (d2=dim;d2<2*dim-1;d2++) {
    //find max over is==1
    maxx(mis,is,d2,&im,&jm,&mm); clmis[d2-dim]=mm;
    is[im]=is[jm]=0;
    is[d2]=1;
    cls[0][d2-dim]=im;cls[1][d2-dim]=jm;
    clmi[d2-dim]=mi[im][jm]; //save not divided mi
    for (d=0;d<=d2;d++) {
	mi[d][d2]=MIM;mi[d2][d]=MIM;
	mis[d][d2]=MIM;mis[d2][d]=MIM;
    }
    dx=0;
    while (noden[im][dx]!=-1) {noden[d2][dx]=noden[im][dx];dx++;}
    dy=0;
    while (noden[jm][dy]!=-1) {noden[d2][dx+dy]=noden[jm][dy];dy++;}
    dmxy = (dx>dy) ? dx : dy ;
    
    if ( (im<dim) && (jm<dim) ) red[d2-dim]=mi[im][jm]; 
    else {
      for (d=0;d<dx+dy;d++) {
	memcpy(xx[d],x[noden[d2][d]],N*sizeof(double));
	scalxx[d]=scal[noden[d2][d]];
      }
      
      redr(xx,dx+dy,N,K,psi,scalxx,&t_d);
      red[d2-dim]=t_d;
    }
    
    dx+=dy; //size of new node
    for (d=0;d<d2;d++) {
      if (is[d]) {
	dy=0;
	while (noden[d][dy]!=-1) dy++;
	for (d1=0;d1<dx;d1++) {
	  memcpy(xx[d1],x[noden[d2][d1]],N*sizeof(double));
	  scalxx[d1]=scal[noden[d2][d1]];
	}
	for (d1=0;d1<dy;d1++) {
	  memcpy(xx[d1+dx],x[noden[d][d1]],N*sizeof(double));
	  scalxx[d1+dx]=scal[noden[d][d1]];
	}
	
	mir_xnyn(xx,dx,dy,N,K,psi,scalxx,&(t_d)); 
	
	if (t_d<0) mi[d][d2]=MIM+1e-7;
	else mi[d][d2]=t_d;
	dmxy = (dx>dy) ? dx : dy;
	mis[d][d2]=mi[d][d2]/(dx+dy);
      }
    }
  }
  for (d=0;d<dim-1;d++) {
    fprintf(stdout,"%d\t%d\t%f\t%f\t%f\n",
	    cls[0][d]+1,cls[1][d]+1,clmis[d],clmi[d],red[d]);
  }

  for (d=0;d<dim;d++) {
    free(x[d]);free(xx[d]);
  }
  free(x);free(xx);
  free(min);free(max);
  free(psi);free(scal);free(scalxx);
  free(is);

  for (d=0;d<2*dim-1;d++) {
    free(noden[d]);free(mi[d]);free(mis[d]);
  }
  free(noden);free(mi);free(mis);
  
  free(cls[0]);free(cls[1]);
  free(clmi);free(clmis);
  free(red);
}

void maxx(double **mi, char *is, int dim, int *im, int *jm, double *mm) {
  int d,dd;
  *mm=MIM;
  for (d=0;d<dim;d++) for (dd=d+1;dd<dim;dd++) 
    if (is[d]*is[dd]) if (mi[d][dd]>*mm) { *im=d;*jm=dd;*mm=mi[*im][*jm]; }
}

