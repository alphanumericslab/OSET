//####################################################################
//MILCA with delay embedding
//####################################################################
//2 May 2004
//Contact: kraskov@its.caltech.edu
//####################################################################
//uses mic_xnyn mir_xnyn
//####################################################################
//Nov 18, 2005
//fixed a bug reported by Yasuhiro Ishikawa,
//the memory for variables max and min was allocated incorrectly
//that is why it crashed if edim was smaller than dim/2

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "miutils.h"

#define M 7
#define NSTACK 50
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr

void four1(double *data, unsigned long nn, int isign);
void realft(double *data, unsigned long n, int isign);


int main(int argc, char **argv) {

  FILE *fin;
  int i,d,d1,d2, d_,d2_,ddel;
  int N=4096,N_;
  int dim=3;
  int K=1;
  
  double **x;
  double **xx,**yy;
  double *psi,*min,*max,*scalxx;
  double t_d,t_d1;
  int BOX1;

  double s,me;


  int k;
  double *rmi;
  double angle;
  double st_d, ct_d;

  double **Rall,**R,**Rt;
  int p, pMax;
  int count;

  double mimin;
  int angI;

  double addnoise=-1;
  int method=1;
  int tau=1;
  int edim=2;
  int ns=1;
  int nr=128;
  
  if (argc<7) {
    fprintf(stderr,"\nMutual Information Least dependent Component Analysis (MILCA) with delay embedding\n\n");
    fprintf(stderr,"Usage:\n%s <filename> <dim> <# points> <# neighbours> <tau> <edim> [method] [nrot] [nharm] [addnoise]\n\n",argv[0]);
    fprintf(stderr,"Input:\n\t<filename>\ttext file with <dim> columns and <# points> rows\n");
    fprintf(stderr,"\t<dim>\t\tnumber of columns in file\n");
    fprintf(stderr,"\t<# points>\tnumber of rows (length of characteristic vector)\n");
    fprintf(stderr,"\t<# neighbours>\tnumber of the nearest neighbours for MI estimator\n");
    fprintf(stderr,"\t<tau>\t\tembedding delay; default 1\n");
    fprintf(stderr,"\t<edim>\t\tembedding dimension; default 2\n");
    fprintf(stderr,"\t[method]\teither cubic (0), or rectangular (1); default 1\n");
    fprintf(stderr,"\t[nrot]\t\tnumber of angles to minimize over; default 128\n");
    fprintf(stderr,"\t[nharm]\t\tnumber of harmonics to approximate MI(angle) dependence; default 1\n");
    fprintf(stderr,"\t[addnoise]\tnoise amplitude; default 1e-8\n");
    fprintf(stderr,"\nOutput:\n");
    fprintf(stderr,"\nde-mixing matrix\n");
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
  if (argc>=10) ns=atoi(argv[9]);
  if (argc>=11) addnoise=atof(argv[10]);
  if (argc>=12) {fprintf(stderr,"Too many input arguments\n");exit(-1);}

  pMax=dim-1;

  rmi=(double *)calloc(nr,sizeof(double));

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

  Rall=(double **)calloc(dim,sizeof(double*));
  R=(double **)calloc(dim,sizeof(double*));
  for (d=0;d<dim;d++) {
    Rall[d]=(double *)calloc(dim,sizeof(double));
    R[d]=(double *)calloc(dim,sizeof(double));
    for (d2=0;d2<dim;d2++) Rall[d][d2]=(d==d2);
  }
  Rt=(double **)calloc(2,sizeof(double*));
  Rt[0]=(double *)calloc(dim,sizeof(double));
  Rt[1]=(double *)calloc(dim,sizeof(double));
  
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

  min=(double*)calloc(dim,sizeof(double));
  max=(double*)calloc(dim,sizeof(double));
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
  free(min); free(max);
  
  psi=(double*)calloc(N+1,sizeof(double)); 
  psi[1]=-(double).57721566490153;
  for (i=1;i<N;i++) psi[i+1]=psi[i]+1/(double)i; //cubic
  BOX1=N-5;

  scalxx=(double*)calloc(2*edim,sizeof(double));
  min=(double*)calloc(2*edim,sizeof(double));
  max=(double*)calloc(2*edim,sizeof(double));


  count=0;

  N_=N-(edim-1)*tau;

  for (p=0;p<pMax;p++) {
    for (d=0;d<dim;d++) {
      for (d2=0;d2<dim;d2++) {
	if ( (d<d2) || ((dim>2) && (d!=d2)) ) {
	  for (k=0;k<nr;k++) {
	    angle=1.0*k/nr*M_PI/2;
	    st_d=sin(angle);ct_d=cos(angle);
	    for (i=0;i<N;i++) {
	      t_d=x[d][i];t_d1=x[d2][i];
	      yy[0][i]=ct_d*t_d+st_d*t_d1;
	      yy[1][i]=ct_d*t_d1-st_d*t_d;
	    }
	    for (ddel=0;ddel<edim;ddel++) {
	      memcpy(xx[ddel],yy[0]+ddel*tau,N_*sizeof(double));
	      memcpy(xx[ddel+edim],yy[1]+ddel*tau,N_*sizeof(double));
	    }
	    
	    for (d1=0;d1<2*edim;d1++) {min[d1]=DBL_MAX/2;max[d1]=-DBL_MAX/2;}

	    for (d1=0;d1<2*edim;d1++) {
	      for (i=0;i<N_;i++) {
		if (xx[d1][i]<min[d1]) min[d1]=xx[d1][i]; 
		if (xx[d1][i]>max[d1]) max[d1]=xx[d1][i];
	      }
	      for (i=0;i<N_;i++) xx[d1][i]=xx[d1][i]-min[d1];
	    }
	    
	    for (d1=0;d1<2*edim;d1++) scalxx[d1]=BOX1/(max[d1]-min[d1]); 
      
	    switch (method) {
	    case 0 :  mic_xnyn(xx,edim,edim,N_,K,psi,scalxx,&(rmi[k])); break;
	    case 1 :  mir_xnyn(xx,edim,edim,N_,K,psi,scalxx,&(rmi[k])); break;
	    }
	  }
	  realft(rmi-1,nr,1);
	  for (k=2+ns*2;k<nr;k++) rmi[k]=0;

	  realft(rmi-1,nr,-1);
	  mimin=rmi[0];angI=0;
	  for (k=1;k<nr;k++) {
	    if (rmi[k]<mimin) {
	      mimin=rmi[k];angI=k;
	    }
	  }
	  angle=1.0*angI/nr*M_PI/2;
	  
	  if ( (angI>=1) && (angI<=nr-1) ) {//minimun angle rotation are too small, and not considered
	    for (d_=0;d_<dim;d_++) for (d2_=0;d2_<dim;d2_++) R[d_][d2_]=(d_==d2_);
	    
	    st_d=sin(angle);ct_d=cos(angle);
	    
	    R[d][d]=ct_d;
	    R[d][d2]=st_d;
	    R[d2][d]=-st_d;
	    R[d2][d2]=ct_d;
	    
	    for (i=0;i<N;i++) {
	      t_d=x[d][i];t_d1=x[d2][i];
	      x[d][i]=ct_d*t_d+st_d*t_d1;
	      x[d2][i]=ct_d*t_d1-st_d*t_d;
	    }
	    for (d_=0;d_<dim;d_++) {
	      Rt[0][d_]=0;Rt[1][d_]=0;
	      for (d2_=0;d2_<dim;d2_++) {
		Rt[0][d_]+=R[d][d2_]*Rall[d2_][d_];
		Rt[1][d_]+=R[d2][d2_]*Rall[d2_][d_];
	      }
	    }
	    memcpy(Rall[d],Rt[0],dim*sizeof(double));
	    memcpy(Rall[d2],Rt[1],dim*sizeof(double));
	    count++;
	  }
	}
      }
    }     

    if (dim==2) break;
    if (count==0) break;
    count=0;
  }

  for (d=0;d<dim;d++) {
    for (d2=0;d2<dim;d2++) {
      fprintf(stdout,"%f\t",Rall[d][d2]);
    }
    fprintf(stdout,"\n");
  }

  for (d=0;d<dim;d++) {
    free(x[d]);
    free(Rall[d]);
    free(R[d]);
  }

  free(Rt[0]);free(Rt[1]);
  free(Rall);free(R);
  free(Rt);
  free(yy[0]);free(yy[1]);free(yy);
  free(xx[0]);free(xx[1]);free(xx);
  free(x);
  free(min);free(max);
  free(psi);free(scalxx);
  free(rmi);
}

void four1(double *data, unsigned long nn, int isign)
  // The length of data should be 2*nn, because data is complex array
  // nn - number of complex points
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software -0)0+3$j3D.. */
void realft(double *data, unsigned long n, int isign) { 
  //Calculates the Fourier transform of a set of n real-valued data
  //points. Replaces this data (which is stored in array data[1..n]) by
  //the positive frequency half of its complex Fourier transform. The
  //real-valued  rst and last components of the complex transform are
  //returned as elements data[1] and data[2], respectively. n must be a
  //power of 2. This routine also calculates the inverse transform of a
  //complex data array if it is the transform of real data. (Result in
  //this case must be multiplied by 2/n.)
  unsigned long i,i1,i2,i3,i4,np3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  
  theta=M_PI/(double)(n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software -0)0+3$j3D.. */
