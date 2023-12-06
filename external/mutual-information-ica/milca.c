//####################################################################
//MILCA
//####################################################################
//2 May 2004
//Contact: kraskov@its.caltech.edu
//####################################################################
//uses mi2c mi2r
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
  int i,d,d1,d2, d_,d2_;
  int N=4096;;
  int dim=3;
  int K=1;

  double **x;
  double **xx;
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
  int ns=1;

  double mimin;
  int angI;

  int method=1;
  int nr=128;
  double addnoise=-1;
  
  if (argc<5) {
    fprintf(stderr,"\nMutual Information Least dependent Component Analysis (MILCA)\n\n");
    fprintf(stderr,"Usage:\n%s <filename> <dim> <# points> <# neighbours> [method] [nrot] [nharm] [addnoise]\n\n",argv[0]);
    fprintf(stderr,"Input:\n\t<filename>\ttext file with <dim> columns and <# points> rows\n");
    fprintf(stderr,"\t<dim>\t\tnumber of columns in file\n");
    fprintf(stderr,"\t<# points>\tnumber of rows (length of characteristic vector)\n");
    fprintf(stderr,"\t<# neighbours>\tnumber of the nearest neighbours for MI estimator\n");
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
  if (argc>=6) method=atoi(argv[5]);
  if (argc>=7) nr=atoi(argv[6]);
  if (argc>=8) ns=atoi(argv[7]);
  if (argc>=9) addnoise=atof(argv[8]);
  if (argc>=10) {fprintf(stderr,"Too many input arguments\n");exit(-1);}


  pMax=dim-1;

  rmi=(double *)calloc(nr+1,sizeof(double));

  xx=(double **)calloc(2,sizeof(double*));
  xx[0]=(double *)calloc(N,sizeof(double));
  xx[1]=(double *)calloc(N,sizeof(double));

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
    srand((dim)*N*K*int(x[(dim)/2][N/10]));
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

  psi=(double*)calloc(N+1,sizeof(double));
  psi[1]=-(double).57721566490153;
  for (i=1;i<N;i++) psi[i+1]=psi[i]+1/(double)i; //cubic
  BOX1=N-5;
  scalxx=(double*)calloc(2,sizeof(double));

  count=0;
  for (p=0;p<pMax;p++) {
    for (d=0;d<dim;d++) {
      for (d2=0;d2<dim;d2++) {
	if ( (d<d2) || ((dim>2) && (d!=d2)) ) {
	  for (k=0;k<nr;k++) {
	    angle=1.0*k/nr*M_PI/2;
	    st_d=sin(angle);ct_d=cos(angle);
	    for (i=0;i<N;i++) {
	      t_d=x[d][i];t_d1=x[d2][i];
	      xx[0][i]=ct_d*t_d+st_d*t_d1;
	      xx[1][i]=ct_d*t_d1-st_d*t_d;
	    }
	    for (d1=0;d1<2;d1++) {min[d1]=DBL_MAX/2;max[d1]=-DBL_MAX/2;}
	    for (d1=0;d1<2;d1++) {
	      for (i=0;i<N;i++) {
		if (xx[d1][i]<min[d1]) min[d1]=xx[d1][i];
		if (xx[d1][i]>max[d1]) max[d1]=xx[d1][i];
	      }
	      for (i=0;i<N;i++) xx[d1][i]=xx[d1][i]-min[d1];
	    }
	    for (d1=0;d1<2;d1++) scalxx[d1]=BOX1/(max[d1]-min[d1]);
	    switch (method) {
	    case 0 :  mi2c(xx,N,K,psi,scalxx,&(rmi[k])); break;
	    case 1 :  mi2r(xx,N,K,psi,scalxx,&(rmi[k])); break;
	    }
	  }
	  realft(rmi-1,nr,1);
//	  rmi[0]=0;rmi[1]=0;
	  for (k=2+ns*2;k<nr;k++) rmi[k]=0;
	  realft(rmi-1,nr,-1);
	  mimin=rmi[0];angI=0;
	  for (k=1;k<nr;k++) {
	    if (rmi[k]<mimin) {
	      mimin=rmi[k];angI=k;
	    }
	  }

	  angle=1.0*angI/nr*M_PI/2;

	  //fprintf(stdout,"i=%d, j=%d \t p=%d, angle %1.6f\n",d,d2,p,angle);
	  
	  if ( (angI>=1) && (angI<=nr-1) ) {//minumin angle rotation are too small, and not considered
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
  //output
  for (d=0;d<dim;d++) {
    for (d2=0;d2<dim;d2++) {
      fprintf(stdout,"%f\t",Rall[d][d2]);
    }
    fprintf(stdout,"\n");
  }

  free(rmi);
  free(xx[0]);free(xx[1]);free(xx);
  for (d=0;d<dim;d++) {
    free(Rall[d]);free(R[d]);free(x[d]);
  }
  free(Rall);free(R);free(x);
  free(Rt[0]);free(Rt[1]);free(Rt);
  free(min);
  free(max);
  free(psi);
  free(scalxx);
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
