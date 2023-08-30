#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void make_box1(double *x, int N, double scal, int bs, 
	       int *box, int *lis, int *mxi) {
  int i, ix;
  for (i=0;i<=bs;i++) box[i]=-1;
  for (i=0;i<=bs;i++) mxi[i]=0;
  for(i=0;i<N;i++) {ix=(int)(x[i]*scal);lis[i]=box[ix]; box[ix]=i; mxi[ix]++;}
  for (i=1;i<=bs;i++) mxi[i]+=mxi[i-1];
}
void make_box2(double **x, int dim, int N, int comp1, int comp2, int bs, int inveps, 
	       int **box, int *lis) {
  int d,i,ix,iy,ixy;
  int ib=bs-1;
  double **xx;

  xx=(double **)calloc(dim, sizeof(double*)); 
  for (d=0;d<dim;d++) {
    xx[d]=(double *)calloc(N, sizeof(double)); 
    memcpy(xx[d],x[d],N*sizeof(double));
  }
  for (ix=0;ix<bs;ix++) for (iy=0;iy<bs;iy++) box[ix][iy] = -1;
  
  for (i=0;i<N;i++) {
    ix=(int)(x[comp1][i]*inveps)&ib;iy=(int)(x[comp2][i]*inveps)&ib;
    lis[i]=box[ix][iy];box[ix][iy]=i;
  }
  i=-1;
  for (ix=0;ix<bs;ix++) for (iy=0;iy<bs;iy++) {
    ixy=box[ix][iy];
    while(ixy>=0) {
      i++;
      for (d=0;d<dim;d++) x[d][i]=xx[d][ixy];
      ixy=lis[ixy];
    }
    box[ix][iy]=-1;
  }
  for (i=0;i<N;i++) {
    ix=(int)(x[comp1][i]*inveps)&ib;iy=(int)(x[comp2][i]*inveps)&ib;
    lis[i]=box[ix][iy];box[ix][iy]=i;
  }
  for (d=0;d<dim;d++) free(xx[d]); free(xx);
}
void make_box2ind(double **x, int dim, int N, int comp1, int comp2, int bs, int inveps, 
		  int *ind,  int **box, int *lis) {
  int d,i,ix,iy,ixy;
  int ib=bs-1;
  double **xx;

  xx=(double **)calloc(dim, sizeof(double*)); 
  for (d=0;d<dim;d++) {
    xx[d]=(double *)calloc(N, sizeof(double)); 
    memcpy(xx[d],x[d],N*sizeof(double));
  }
  for (ix=0;ix<bs;ix++) for (iy=0;iy<bs;iy++) box[ix][iy] = -1;
  
  for (i=0;i<N;i++) {
    ix=(int)(x[comp1][i]*inveps)&ib;iy=(int)(x[comp2][i]*inveps)&ib;
    lis[i]=box[ix][iy];box[ix][iy]=i;
  }
  i=-1;
  for (ix=0;ix<bs;ix++) for (iy=0;iy<bs;iy++) {
    ixy=box[ix][iy];
    while(ixy>=0) {
      i++;
      for (d=0;d<dim;d++) x[d][i]=xx[d][ixy];
      ind[ixy]=i;ixy=lis[ixy];
    }
    box[ix][iy]=-1;
  }
  for (i=0;i<N;i++) {
    ix=(int)(x[comp1][i]*inveps)&ib;iy=(int)(x[comp2][i]*inveps)&ib;
    lis[i]=box[ix][iy];box[ix][iy]=i;
  }
  for (d=0;d<dim;d++) free(xx[d]); free(xx);
}
int neiE1(double *x, int i, double scal, int bs, double eps,  int *box, int *lis, int *mxi) {
  double dd;
  int mi,mp,mm,nx=0;
  double xc=x[i];

  mp=int((xc+eps)*scal); if(mp>bs) mp=bs;
  mm=int((xc-eps)*scal); if(mm<0)  mm=0;
  mi=box[mp];
  while(mi>=0) {
    dd=x[mi]-xc;if(fabs(dd)<=eps) nx++;
    mi=lis[mi];
  }
  if(mm>=mp) return nx-1;
  mi=box[mm];
  while(mi>=0) {
    dd=xc-x[mi];if(fabs(dd)<=eps) nx++;
    mi=lis[mi];
  }
  nx+=mxi[mp-1]-mxi[mm];
  return nx-1;
}
int neiE(double **x, int i, int comp1, int comp2, int dim, int bs, double epsgrid, double eps, int **box, int *lis) {
  int ix,iy,ix1,iy1,ix2,jj,step,d;
  int el,nx,ib=bs-1;
  double dd,dy;
  double *xx;

  xx=(double *)calloc(dim,sizeof(double));
  for (d=0;d<dim;d++) xx[d]=x[d][i];

  ix=(int)(xx[comp1]/epsgrid)&ib; iy=(int)(xx[comp2]/epsgrid)&ib;
  jj=0; nx=0;
  while (eps>epsgrid*(jj-1)) {
    step = (jj) ? 2*jj : 1;
    for (ix1=ix-jj;ix1<=ix+jj;ix1++) {
      ix2=ix1&ib;
      for (iy1=iy-jj;iy1<=iy+jj;iy1+=step) {
	el=box[ix2][iy1&ib];
	while (el != -1) { 
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps ) break;
	  if (dd<=eps) nx++;
	  el=lis[el];
	}
      }
    }
    for (ix1=ix-jj;ix1<=ix+jj;ix1+=step) {
      ix2=ix1&ib;
      for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
	el=box[ix2][iy1&ib];
	while (el != -1) { 
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps ) break;
	  if (dd<=eps) nx++;
	  el=lis[el];
	}
      }
    }
    jj++;
    if (jj==(bs/2)) break;
  }
  if ( jj==(bs/2) ) { //half of the layer
    for (ix1=ix-jj;ix1<ix+jj;ix1++) {
      ix2=ix1&ib; iy1=iy-jj;
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	dd=fabs(xx[0]-x[0][el]); 
	for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps ) break;
	if (dd<=eps) nx++;
	el=lis[el];
      }
    }
    ix1=ix-jj; ix2=ix1&ib;
    for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	dd=fabs(xx[0]-x[0][el]); 
	for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps ) break;
	if (dd<=eps) nx++;
	el=lis[el];
      }
    }
  }
  free(xx);
  return nx-1;
}

void neiEK(double **x, int i, int comp1, int comp2, int dim, int K,
	   int bs, double epsgrid, double *eps, int **box, int *lis,
	   int *nx) {
  int ix,iy,ix1,iy1,ix2,jj,step,d,ik;
  int el,ib=bs-1;
  double dd,dy;
  double *xx;

  xx=(double *)calloc(dim,sizeof(double));
  for (d=0;d<dim;d++) xx[d]=x[d][i];

  ix=(int)(xx[comp1]/epsgrid)&ib; iy=(int)(xx[comp2]/epsgrid)&ib;
  jj=0; 
  for (ik=0;ik<K;ik++) nx[ik]=0;
  while (eps[K-1]>epsgrid*(jj-1)) {
    step = (jj) ? 2*jj : 1;
    for (ix1=ix-jj;ix1<=ix+jj;ix1++) {
      ix2=ix1&ib;
      for (iy1=iy-jj;iy1<=iy+jj;iy1+=step) {
	el=box[ix2][iy1&ib];
	while (el != -1) { 
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps[K-1] ) break;; 
	  for (ik=0;ik<K;ik++) {
	    if (dd<=eps[ik]) nx[ik]++;
	  }
	  el=lis[el];
	}
      }
    }
    for (ix1=ix-jj;ix1<=ix+jj;ix1+=step) {
      ix2=ix1&ib;
      for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
	el=box[ix2][iy1&ib];
	while (el != -1) { 
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps[K-1] ) break;
	  for (ik=0;ik<K;ik++) {
	    if (dd<=eps[ik]) nx[ik]++;
	  }
	  el=lis[el];
	}
      }
    }
    jj++;
    if (jj==(bs/2)) break;
  }
  if ( jj==(bs/2) ) { //half of the layer
    for (ix1=ix-jj;ix1<ix+jj;ix1++) {
      ix2=ix1&ib; iy1=iy-jj;
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	dd=fabs(xx[0]-x[0][el]); 
	for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps[K-1] ) break;
	for (ik=0;ik<K;ik++) {
	  if (dd<=eps[ik]) nx[ik]++;
	}
	el=lis[el];
      }
    }
    ix1=ix-jj; ix2=ix1&ib;
    for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	dd=fabs(xx[0]-x[0][el]); 
	for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) if ( (dd=dy) > eps[K-1] ) break;
	for (ik=0;ik<K;ik++) {
	  if (dd<=eps[ik]) nx[ik]++;
	}
	el=lis[el];
      }
    }
  }
  free(xx);
  for (ik=0;ik<K;ik++) nx[ik]-=1;
}


void neiK(double **x, int dim, int comp1, int comp2, int i, 
	  int bs, double epsgrid, int K, int **box, int *lis,
	  int *nn) {
  double *dn,*xx;
  int k,ix,iy,ix1,iy1,ix2,jj,step,ib=bs-1;
  int el;
  double dd,dy;
  int d;

  dn=(double*)calloc(K+1,sizeof(double));
  xx=(double*)calloc(dim,sizeof(double));
  for(k=0;k<dim;k++) xx[k]=x[k][i];
  dn[0]=0;
  for(k=1;k<=K;k++) dn[k]=1e30;
  ix=(int)(xx[comp1]/epsgrid)&ib; iy=(int)(xx[comp2]/epsgrid)&ib;
  jj=0; 
  while (dn[K]>epsgrid*(jj-1)) {
    step = (jj) ? 2*jj : 1;
    for (ix1=ix-jj;ix1<=ix+jj;ix1++) {
      ix2=ix1&ib;
      for (iy1=iy-jj;iy1<=iy+jj;iy1+=step) {
	el=box[ix2][iy1&ib];
	while (el != -1) { if (el!=i) {
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) dd=dy;
	  if(dd<dn[K]) {
	    k=K; while(dd<dn[k]) { if (k<K) { dn[k+1]=dn[k]; nn[k+1]=nn[k]; } k--; }
	    dn[k+1]=dd; nn[k+1]=el;
	  }
	} el=lis[el];
	}
      }
    }
    for (ix1=ix-jj;ix1<=ix+jj;ix1+=step) {
      ix2=ix1&ib;
      for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
	el=box[ix2][iy1&ib];
	while (el != -1) { if (el!=i) {
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) dd=dy;
	  if ( (dy=fabs(xx[1]-x[1][el]))>dd ) dd=dy;
	  if(dd<dn[K]) {
	    k=K; while(dd<dn[k]) { if (k<K) { dn[k+1]=dn[k]; nn[k+1]=nn[k]; } k--; }
	    dn[k+1]=dd; nn[k+1]=el;
	  } 
	} el=lis[el];
	}
      }
    }
    jj++;
    if (jj==bs/2) break;
  }

  if ( jj==(bs/2) ) { //half of the layer
    for (ix1=ix-jj;ix1<ix+jj;ix1++) {
      ix2=ix1&ib; iy1=iy-jj;
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	if (el!=i) {
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) dd=dy;
	  if(dd<dn[K]) {
	    k=K; while(dd<dn[k]) { if (k<K) { dn[k+1]=dn[k]; nn[k+1]=nn[k]; } k--; }
	    dn[k+1]=dd; nn[k+1]=el;
	  }
	}
	el=lis[el];
      }
    }
    ix1=ix-jj; ix2=ix1&ib;
    for (iy1=iy-jj+1;iy1<=iy+jj-1;iy1++) {
      el=box[ix2][iy1&ib];
      while (el != -1) { 
	if (el!=i) {
	  dd=fabs(xx[0]-x[0][el]); 
	  for (d=1;d<dim;d++) if ( (dy=fabs(xx[d]-x[d][el]))>dd ) dd=dy;
	  if(dd<dn[K]) {
	    k=K; while(dd<dn[k]) { if (k<K) { dn[k+1]=dn[k]; nn[k+1]=nn[k]; } k--; }
	    dn[k+1]=dd; nn[k+1]=el;
	  }
	} 
	el=lis[el];
      }
    }
  }
  free(dn);free(xx);
}


/*********/
void mi2(double **x, int N, int K,
	 double *psi, 
	 double *scal,
	 double *mic, double *mir) {
  
  int i,k,*n1,*n2;
  double *xc,dx;
  double *eps,Eps;
  int maxdim;
  double dxy1,dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  int dim=2;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 


  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n1=(int*)calloc(dim,sizeof(int));
  n2=(int*)calloc(dim,sizeof(int));

  make_box2(x,dim,N,0,1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  dxy1=dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    Eps=0;maxdim=-1;
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for(k=1;k<=K;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[d] ) eps[d]=dx; }
      if (eps[d]>Eps) {Eps=eps[d];maxdim=d;}
    }	
    for (d=0;d<dim;d++) {
      n2[d]=neiE1(x[d],i,scal[d],BOX1,eps[d],box1[d],lis1[d],mxi[d]);
      if (d==maxdim) { n1[d]=n2[d]; dxy1+=psi[n1[d]]; }
      else { 
	n1[d]=neiE1(x[d],i,scal[d],BOX1,Eps,box1[d],lis1[d],mxi[d]); 
	dxy1+=psi[n1[d]+1];
      }
      dxy2+=psi[n2[d]];
    }
    //    fprintf(stdout,"%f %f %d %d %d %d %d %f %d\n",xc[0],xc[1],(eps[0]<=eps[1]),n1[0],n1[1],n2[0],n2[1],Eps,nn[K]);
  }
  dxy1/=N;*mic=psi[N]+psi[K]-dxy1;
  dxy2/=N;*mir=psi[N]+phi[K]-dxy2;
  
  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  free(eps);free(n1);free(n2);
  free(phi);
}

void mi2c(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mic) {
  
  int i,k,*n1;
  double *xc,dx;
  double *eps,Eps;
  int maxdim;
  double dxy1;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  int dim=2;

  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n1=(int*)calloc(dim,sizeof(int));

  make_box2(x,dim,N,0,1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  dxy1=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    Eps=0;maxdim=-1;
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for(k=1;k<=K;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[d] ) eps[d]=dx; }
      if (eps[d]>Eps) {Eps=eps[d];maxdim=d;}
    }	
    for (d=0;d<dim;d++) {
      n1[d]=neiE1(x[d],i,scal[d],BOX1,Eps,box1[d],lis1[d],mxi[d]); 
      if (d==maxdim) { 
	dxy1+=psi[n1[d]];
      } else { 
	dxy1+=psi[n1[d]+1];
      }
    }
    //    fprintf(stdout,"%f %f %d %d %d %f %d\n",xc[0],xc[1],(eps[0]<=eps[1]),n1[0],n1[1],Eps,nn[K]);
  }
  dxy1/=N;*mic=psi[N]+psi[K]-dxy1;
  
  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  free(eps);free(n1);
}

void mi2r(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mir) {
  
  int i,k,*n2;
  double *xc,dx;
  double *eps;
  double dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  int dim=2;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 


  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n2=(int*)calloc(dim,sizeof(int));

  make_box2(x,dim,N,0,1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for(k=1;k<=K;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[d] ) eps[d]=dx; }
    }	
    for (d=0;d<dim;d++) {
      n2[d]=neiE1(x[d],i,scal[d],BOX1,eps[d],box1[d],lis1[d],mxi[d]);
      dxy2+=psi[n2[d]];
    }
    //    fprintf(stdout,"%f %f %d %d %d %d\n",xc[0],xc[1],(eps[0]<=eps[1]),n2[0],n2[1],nn[K]);
  }
  dxy2/=N;*mir=psi[N]+phi[K]-dxy2;
  
  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  free(eps);free(n2);
  free(phi);
}

void red(double **x, int dim, int N, int K, 
	 double *psi,
	 double *scal,
	 double *mic, double *mir) {
  int i,k,*n1,*n2;
  double *xc,dx;
  double *eps,Eps;
  int maxdim;
  double dxy1,dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 


  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n1=(int*)calloc(dim,sizeof(int));
  n2=(int*)calloc(dim,sizeof(int));

  make_box2(x,dim,N,0,dim-1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  dxy1=dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    //    for (k=1;k<=K;k++) {
    //      fprintf(stdout,"k%d %f %f %f %f\n",k,x[0][i],x[0][nn[k]],x[1][i],x[1][nn[k]]);
    //    }
    
    Eps=0;maxdim=-1;
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for(k=1;k<=K;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[d] ) eps[d]=dx; }
      if (eps[d]>Eps) {Eps=eps[d];maxdim=d;}
    }	
    for (d=0;d<dim;d++) {
      n2[d]=neiE1(x[d],i,scal[d],BOX1,eps[d],box1[d],lis1[d],mxi[d]);
      if (d==maxdim) { n1[d]=n2[d]; dxy1+=psi[n1[d]]; }
      else { 
	n1[d]=neiE1(x[d],i,scal[d],BOX1,Eps,box1[d],lis1[d],mxi[d]); 
	dxy1+=psi[n1[d]+1];
      }
      dxy2+=psi[n2[d]];
    }
    //    fprintf(stdout,"%f %f %e %e %e %d %d %d %d %f %d\n",xc[0],xc[1],eps[0],eps[1],(eps[0]-eps[1]),n1[0],n1[1],n2[0],n2[1],Eps,nn[K]);
  }
  dxy1/=N;*mic=(dim-1)*psi[N]+psi[K]-dxy1;
  dxy2/=N;*mir=(dim-1)*psi[N]+phi[K]-dxy2;
  
  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  free(eps);free(n1);free(n2);
  free(phi);
}


void redc(double **x, int dim, int N, int K, 
	  double *psi,
	  double *scal,
	  double *mic) {
  
  int i,k,*n1;
  double *xc,dx;
  double *eps,Eps;
  int maxdim;
  double dxy1;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n1=(int*)calloc(dim,sizeof(int));

  make_box2(x,dim,N,0,dim-1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  dxy1=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    Eps=0;maxdim=-1;
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for(k=1;k<=K;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[d] ) eps[d]=dx; }
      if (eps[d]>Eps) {Eps=eps[d];maxdim=d;}
    }	
    for (d=0;d<dim;d++) {
      n1[d]=neiE1(x[d],i,scal[d],BOX1,Eps,box1[d],lis1[d],mxi[d]); 
      if (d==maxdim) { 
	dxy1+=psi[n1[d]];
      } else { 
	dxy1+=psi[n1[d]+1];
      }
    }
    //    fprintf(stdout,"%f %f %d %d %d %f %d\n",xc[0],xc[1],(eps[0]<=eps[1]),n1[0],n1[1],Eps,nn[K]);
  }
  dxy1/=N;*mic=(dim-1)*psi[N]+psi[K]-dxy1;
  
  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  free(eps);free(n1);
}

  


void redr(double **x, int dim, int N, int K, 
	  double *psi,
	  double *scal,
	  double *mir) {
  int i,k,*n2;
  double *xc,dx;
  double *eps;
  double dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 


  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n2=(int*)calloc(dim,sizeof(int));

  make_box2(x,dim,N,0,dim-1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for(k=1;k<=K;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[d] ) eps[d]=dx; }
    }	
    for (d=0;d<dim;d++) {
      n2[d]=neiE1(x[d],i,scal[d],BOX1,eps[d],box1[d],lis1[d],mxi[d]);
      dxy2+=psi[n2[d]];
    }
    //    fprintf(stdout,"%f %f %d %d %d %d\n",xc[0],xc[1],(eps[0]<=eps[1]),n2[0],n2[1],nn[K]);
  }
  dxy2/=N;*mir=(dim-1)*psi[N]+phi[K]-dxy2;
  
  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  free(eps);free(n2);
  free(phi);
}
  

void redr_embed(double **x, int dim, int edim, int tau, int N, int K, 
		double *psi,
		double *mir) {
  int i,k,*n2;
  double *xc,dx;
  double *eps;
  double dxy2;
  int *nn;
  int d,de;
  int N_;

  double **xx;
  double ***xxx;
  
  int BOX;
  int **box,*lis; // two dimensional boxes
  int ***boxxx; // twodimensional boxes for each dim
  int **lisxx; // lists for each dim
  double epsilon;
  int inveps;

  int *ind;
  int **indxx;

  double *phi;

  //  if (edim<2) Just call usual redundancy 

  N_=N-(edim-1)*tau;
  xx=(double **)calloc(dim*edim,sizeof(double*));
  // reservation, and embedding
  for (d=0; d<dim*edim; d++) {
    xx[d]=(double*)calloc(N_,sizeof(double));
    memcpy(xx[d],x[d/edim]+tau*(d%edim),N_*sizeof(double));
  }
  xxx=(double ***)calloc(dim,sizeof(double**));
  // just copy, make bo will destroy an order
  for (d=0; d<dim; d++) {
    xxx[d]=(double**)calloc(edim,sizeof(double*));
    for (de=0; de<edim; de++) {
      xxx[d][de]=(double*)calloc(N_,sizeof(double));
      memcpy(xxx[d][de],xx[d*edim+de],N_*sizeof(double));
    }
  }

  ind=(int*)calloc(N_,sizeof(int));
  phi=(double*)calloc(K+1,sizeof(double));

  for (i=1;i<=K;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 

  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim*edim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  boxxx=(int***)calloc(dim,sizeof(int**)); 
  lisxx=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    boxxx[d]=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxxx[d][i]=(int*)calloc(BOX,sizeof(int));
    lisxx[d]=(int*)calloc(N_,sizeof(int)); 
  }
  indxx=(int**)calloc(dim,sizeof(int*));
  for (d=0;d<dim;d++) {
    indxx[d]=(int*)calloc(N_,sizeof(int)); 
  }


  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n2=(int*)calloc(dim,sizeof(int));

  //for searching neighbours in product (dim*edim) space
  make_box2ind(xx,dim*edim,N_,0,edim,BOX,inveps,ind,box,lis); 
  //second component is arbitratry, lets take the first component of second vector 
  // it should be the same for neiK procedure

  for (d=0;d<dim;d++) 
    make_box2ind(xxx[d],edim,N_,0,edim-1,BOX,inveps,indxx[d],boxxx[d],lisxx[d]);
  
  dxy2=0.0;
  for (i=0;i<N_;i++) {
    for (d=0;d<dim*edim;d++) xc[d]=xx[d][ind[i]];
    neiK(xx,dim*edim,0,edim,ind[i],BOX,epsilon,K,box,lis,nn);
    
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for (de=0;de<edim;de++) 
	for(k=1;k<=K;k++) 
	  if( (dx=fabs(xc[d*edim+de]-xx[d*edim+de][nn[k]]))>eps[d] ) eps[d]=dx;
    }	
    for (d=0;d<dim;d++) {
      n2[d]=neiE(xxx[d],indxx[d][i],0,edim-1,edim,BOX,epsilon,eps[d],boxxx[d],lisxx[d]);
      dxy2+=psi[n2[d]];
    }
    //    fprintf(stdout,"%f %f %d %d %d %d\n",xc[0],xc[1],(eps[0]<=eps[1]),n2[0],n2[1],nn[K]);
  }
  dxy2/=N_;*mir=(dim-1)*psi[N_]+phi[K]-dxy2;

  for (d=0; d<dim*edim; d++) free(xx[d]);
  free(xx);

  for (d=0; d<dim; d++) {
    for (de=0; de<edim; de++) free(xxx[d][de]);
    for (i=0;i<BOX;i++) free(boxxx[d][i]);
    free(boxxx[d]);free(lisxx[d]);free(xxx[d]);free(indxx[d]);
  }
  free(xxx);free(boxxx);free(lisxx);free(indxx);
  free(ind);free(phi);free(nn);
  free(xc);

  for (i=0;i<BOX;i++) free(box[i]);
  free(box);free(lis);
  free(eps);free(n2);
}
  

void mi_xnyn(double **x, int dimx, int dimy, int N, int K, 
	     double *psi, 
	     double *scal,
	     double *mic, double *mir) {
  int i,k,nx1,ny1,nx2,ny2;
  double *xc,dy,dx;
  double epsx,epsy,Eps;
  double dxy1,dxy2;
  double **xx,**yy;;
  double scalx, scaly;
  int *nn;
  int d;

  int BOX,BOX1;
  int **box,**boxy,**boxx,*lis; // two dimensional boxes
  int *lisy,*lisx; // lists for two dimensions
  int *boxx1, *boxy1; // onedimensional boxes
  int *lisy1,*lisx1; // lists for one dimensions
  int *mxi, *myi; //accumulative lists of points in oned boxes
  int *ind,*indx,*indy; //indeces of original data (the data resorted during box creating)
  double epsilon;
  int inveps;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-1/(double(i));   // 

  nn=(int*)calloc(K+1,sizeof(int));

  xc=(double*)calloc(dimx+dimy,sizeof(double));

  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;
  BOX1=N-5;

  if (dimx>1) {
    xx=(double**)calloc(dimx,sizeof(double*));
    for (d=0;d<dimx;d++) xx[d]=(double*)calloc(N,sizeof(double));
    boxx=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxx[i]=(int*)calloc(BOX,sizeof(int));
    lisx=(int*)calloc(N,sizeof(int));
  } else { boxx1=(int*)calloc(BOX1+1,sizeof(int)); 
  lisx1=(int*)calloc(N,sizeof(int)); mxi=(int*)calloc(BOX1+1,sizeof(int)); }
  if (dimy>1) {
    yy=(double**)calloc(dimy,sizeof(double*));
    for (d=0;d<dimy;d++) yy[d]=(double*)calloc(N,sizeof(double));
    boxy=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxy[i]=(int*)calloc(BOX,sizeof(int));
    lisy=(int*)calloc(N,sizeof(int));
  } else { boxy1=(int*)calloc(BOX1+1,sizeof(int)); 
  lisy1=(int*)calloc(N,sizeof(int)); myi=(int*)calloc(BOX1+1,sizeof(int)); }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  ind=(int*)calloc(N,sizeof(int));
  indx=(int*)calloc(N,sizeof(int));
  indy=(int*)calloc(N,sizeof(int));

  //save x if it would be reordered
  if (dimx>1) for (d=0;d<dimx;d++) memcpy(xx[d],x[d],N*sizeof(double)); 
  if (dimy>1) for (d=0;d<dimy;d++) memcpy(yy[d],x[d+dimx],N*sizeof(double));

  make_box2ind(x,dimx+dimy,N,0,dimx,BOX,inveps,ind,box,lis); 
  //for searching neighbours in product space

  if (dimx==1) {scalx=scal[0]; make_box1(x[0],N,scalx,BOX1,boxx1,lisx1,mxi);}
  else make_box2ind(xx,dimx,N,0,dimx-1,BOX,inveps,indx,boxx,lisx); 
  if (dimy==1) {scaly=scal[dimx]; make_box1(x[dimx],N,scaly,BOX1,boxy1,lisy1,myi); }
  else make_box2ind(yy,dimy,N,0,dimy-1,BOX,inveps,indy,boxy,lisy); 

  dxy1=dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dimx+dimy;d++) xc[d]=x[d][ind[i]];

    neiK(x,dimx+dimy,0,dimx,ind[i],BOX,epsilon,K,box,lis,nn);
    epsx=0; for (d=0;d<dimx;d++) 
      for(k=1;k<=K;k++) if( (dx=fabs(xc[d]-x[d][nn[k]]))>epsx ) epsx=dx;
    epsy=0; for (d=dimx;d<dimx+dimy;d++) 
      for(k=1;k<=K;k++) if( (dy=fabs(xc[d]-x[d][nn[k]]))>epsy ) epsy=dy;
    if (dimx>1) nx2=neiE(xx,indx[i],0,dimx-1,dimx,BOX,epsilon,epsx,boxx,lisx);
    else  nx2=neiE1(x[0],ind[i],scalx,BOX1,epsx,boxx1,lisx1,mxi);
    if (dimy>1) ny2=neiE(yy,indy[i],0,dimy-1,dimy,BOX,epsilon,epsy,boxy,lisy);
    else ny2=neiE1(x[dimx],ind[i],scaly,BOX1,epsy,boxy1,lisy1,myi);
    
    if (epsx>epsy) {
      Eps=epsx;nx1=nx2;
      if (dimy>1) ny1=neiE(yy,indy[i],0,dimy-1,dimy,BOX,epsilon,Eps,boxy,lisy);
      else ny1=neiE1(x[dimx],ind[i],scaly,BOX1,Eps,boxy1,lisy1,myi);
      dxy1+=psi[nx1]+psi[ny1+1];
    } else {
      Eps=epsy;ny1=ny2;
      if (dimx>1) nx1=neiE(xx,indx[i],0,dimx-1,dimx,BOX,epsilon,Eps,boxx,lisx);
      else nx1=neiE1(x[0],ind[i],scalx,BOX1,Eps,boxx1,lisx1,mxi);
      dxy1+=psi[nx1+1]+psi[ny1];
    }
    dxy2+=psi[nx2]+psi[ny2];
  }
  dxy1/=N;*mic=psi[N]+psi[K]-dxy1;
  dxy2/=N;*mir=psi[N]+phi[K]-dxy2;

  free(xc);free(nn);free(lis);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(ind);free(indx);free(indy);
  if (dimx==1) {free(mxi);free(boxx1);free(lisx1);} 
  else { for (i=0;i<BOX;i++) free(boxx[i]); free(boxx); free(lisx); for (d=0;d<dimx;d++) free(xx[d]); free(xx); }
  if (dimy==1) {free(myi);free(boxy1);free(lisy1);} 
  else { for (i=0;i<BOX;i++) free(boxy[i]); free(boxy); free(lisy); for (d=0;d<dimy;d++) free(yy[d]); free(yy); }
  free(phi);
}

void mic_xnyn(double **x, int dimx, int dimy, int N, int K, 
	      double *psi, 
	      double *scal,
	      double *mic) {
  int i,k,nx1,ny1;
  double *xc,dy,dx;
  double epsx,epsy,Eps;
  double dxy1;
  double **xx,**yy;;
  double scalx, scaly;
  int *nn;
  int d;

  int BOX,BOX1;
  int **box,**boxy,**boxx,*lis; // two dimensional boxes
  int *lisy,*lisx; // lists for two dimensions
  int *boxx1, *boxy1; // onedimensional boxes
  int *lisy1,*lisx1; // lists for one dimensions
  int *mxi, *myi; //accumulative lists of points in oned boxes
  int *ind,*indx,*indy; //indeces of original data (the data resorted during box creating)
  double epsilon;
  int inveps;

  nn=(int*)calloc(K+1,sizeof(int));

  xc=(double*)calloc(dimx+dimy,sizeof(double));

  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;
  BOX1=N-5;

  if (dimx>1) {
    xx=(double**)calloc(dimx,sizeof(double*));
    for (d=0;d<dimx;d++) xx[d]=(double*)calloc(N,sizeof(double));
    boxx=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxx[i]=(int*)calloc(BOX,sizeof(int));
    lisx=(int*)calloc(N,sizeof(int));
  } else { boxx1=(int*)calloc(BOX1+1,sizeof(int)); 
  lisx1=(int*)calloc(N,sizeof(int)); mxi=(int*)calloc(BOX1+1,sizeof(int)); }
  if (dimy>1) {
    yy=(double**)calloc(dimy,sizeof(double*));
    for (d=0;d<dimy;d++) yy[d]=(double*)calloc(N,sizeof(double));
    boxy=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxy[i]=(int*)calloc(BOX,sizeof(int));
    lisy=(int*)calloc(N,sizeof(int));
  } else { boxy1=(int*)calloc(BOX1+1,sizeof(int)); 
  lisy1=(int*)calloc(N,sizeof(int)); myi=(int*)calloc(BOX1+1,sizeof(int)); }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  ind=(int*)calloc(N,sizeof(int));
  indx=(int*)calloc(N,sizeof(int));
  indy=(int*)calloc(N,sizeof(int));

  //save x if it would be reordered
  if (dimx>1) for (d=0;d<dimx;d++) memcpy(xx[d],x[d],N*sizeof(double)); 
  if (dimy>1) for (d=0;d<dimy;d++) memcpy(yy[d],x[d+dimx],N*sizeof(double));

  make_box2ind(x,dimx+dimy,N,0,dimx,BOX,inveps,ind,box,lis); 
  //for searching neighbours in product space

  if (dimx==1) {scalx=scal[0]; make_box1(x[0],N,scalx,BOX1,boxx1,lisx1,mxi);}
  else make_box2ind(xx,dimx,N,0,dimx-1,BOX,inveps,indx,boxx,lisx); 
  if (dimy==1) {scaly=scal[dimx]; make_box1(x[dimx],N,scaly,BOX1,boxy1,lisy1,myi); }
  else make_box2ind(yy,dimy,N,0,dimy-1,BOX,inveps,indy,boxy,lisy); 

  dxy1=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dimx+dimy;d++) xc[d]=x[d][ind[i]];

    neiK(x,dimx+dimy,0,dimx,ind[i],BOX,epsilon,K,box,lis,nn);
    epsx=0; for (d=0;d<dimx;d++) 
      for(k=1;k<=K;k++) if( (dx=fabs(xc[d]-x[d][nn[k]]))>epsx ) epsx=dx;
    epsy=0; for (d=dimx;d<dimx+dimy;d++) 
      for(k=1;k<=K;k++) if( (dy=fabs(xc[d]-x[d][nn[k]]))>epsy ) epsy=dy;

    if (epsx>epsy) { Eps=epsx;} else {Eps=epsy;}

    if (dimy>1) ny1=neiE(yy,indy[i],0,dimy-1,dimy,BOX,epsilon,Eps,boxy,lisy);
    else ny1=neiE1(x[dimx],ind[i],scaly,BOX1,Eps,boxy1,lisy1,myi);
    if (dimx>1) nx1=neiE(xx,indx[i],0,dimx-1,dimx,BOX,epsilon,Eps,boxx,lisx);
    else nx1=neiE1(x[0],ind[i],scalx,BOX1,Eps,boxx1,lisx1,mxi);

    if (epsx>epsy) {
      dxy1+=psi[nx1]+psi[ny1+1];
    } else {
      dxy1+=psi[nx1+1]+psi[ny1];
    }
  }
  dxy1/=N;*mic=psi[N]+psi[K]-dxy1;

  free(xc);free(nn);free(lis);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(ind);free(indx);free(indy);
  if (dimx==1) {free(mxi);free(boxx1);free(lisx1);} 
  else { for (i=0;i<BOX;i++) free(boxx[i]); free(boxx); free(lisx); for (d=0;d<dimx;d++) free(xx[d]); free(xx); }
  if (dimy==1) {free(myi);free(boxy1);free(lisy1);} 
  else { for (i=0;i<BOX;i++) free(boxy[i]); free(boxy); free(lisy); for (d=0;d<dimy;d++) free(yy[d]); free(yy); }
}


void mir_xnyn(double **x, int dimx, int dimy, int N, int K, 
	      double *psi, 
	      double *scal,
	      double *mir) {
  int i,k,nx2,ny2;
  double *xc,dy,dx;
  double epsx,epsy;
  double dxy2;
  double **xx,**yy;;
  double scalx, scaly;
  int *nn;
  int d;

  int BOX,BOX1;
  int **box,**boxy,**boxx,*lis; // two dimensional boxes
  int *lisy,*lisx; // lists for two dimensions
  int *boxx1, *boxy1; // onedimensional boxes
  int *lisy1,*lisx1; // lists for one dimensions
  int *mxi, *myi; //accumulative lists of points in oned boxes
  int *ind,*indx,*indy; //indeces of original data (the data resorted during box creating)
  double epsilon;
  int inveps;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-1/(double(i));   // 

  nn=(int*)calloc(K+1,sizeof(int));

  xc=(double*)calloc(dimx+dimy,sizeof(double));

  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;
  BOX1=N-5;

  if (dimx>1) {
    xx=(double**)calloc(dimx,sizeof(double*));
    for (d=0;d<dimx;d++) xx[d]=(double*)calloc(N,sizeof(double));
    boxx=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxx[i]=(int*)calloc(BOX,sizeof(int));
    lisx=(int*)calloc(N,sizeof(int));
  } else { boxx1=(int*)calloc(BOX1+1,sizeof(int)); 
  lisx1=(int*)calloc(N,sizeof(int)); mxi=(int*)calloc(BOX1+1,sizeof(int)); }
  if (dimy>1) {
    yy=(double**)calloc(dimy,sizeof(double*));
    for (d=0;d<dimy;d++) yy[d]=(double*)calloc(N,sizeof(double));
    boxy=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxy[i]=(int*)calloc(BOX,sizeof(int));
    lisy=(int*)calloc(N,sizeof(int));
  } else { boxy1=(int*)calloc(BOX1+1,sizeof(int)); 
  lisy1=(int*)calloc(N,sizeof(int)); myi=(int*)calloc(BOX1+1,sizeof(int)); }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  ind=(int*)calloc(N,sizeof(int));
  indx=(int*)calloc(N,sizeof(int));
  indy=(int*)calloc(N,sizeof(int));

  //save x if it would be reordered
  if (dimx>1) for (d=0;d<dimx;d++) memcpy(xx[d],x[d],N*sizeof(double)); 
  if (dimy>1) for (d=0;d<dimy;d++) memcpy(yy[d],x[d+dimx],N*sizeof(double));

  make_box2ind(x,dimx+dimy,N,0,dimx,BOX,inveps,ind,box,lis); 
  //for searching neighbours in product space

  if (dimx==1) {scalx=scal[0]; make_box1(x[0],N,scalx,BOX1,boxx1,lisx1,mxi);}
  else make_box2ind(xx,dimx,N,0,dimx-1,BOX,inveps,indx,boxx,lisx); 
  if (dimy==1) {scaly=scal[dimx]; make_box1(x[dimx],N,scaly,BOX1,boxy1,lisy1,myi); }
  else make_box2ind(yy,dimy,N,0,dimy-1,BOX,inveps,indy,boxy,lisy); 

  dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dimx+dimy;d++) xc[d]=x[d][ind[i]];

    neiK(x,dimx+dimy,0,dimx,ind[i],BOX,epsilon,K,box,lis,nn);
    epsx=0; for (d=0;d<dimx;d++) 
      for(k=1;k<=K;k++) if( (dx=fabs(xc[d]-x[d][nn[k]]))>epsx ) epsx=dx;
    epsy=0; for (d=dimx;d<dimx+dimy;d++) 
      for(k=1;k<=K;k++) if( (dy=fabs(xc[d]-x[d][nn[k]]))>epsy ) epsy=dy;
    if (dimx>1) nx2=neiE(xx,indx[i],0,dimx-1,dimx,BOX,epsilon,epsx,boxx,lisx);
    else  nx2=neiE1(x[0],ind[i],scalx,BOX1,epsx,boxx1,lisx1,mxi);
    if (dimy>1) ny2=neiE(yy,indy[i],0,dimy-1,dimy,BOX,epsilon,epsy,boxy,lisy);
    else ny2=neiE1(x[dimx],ind[i],scaly,BOX1,epsy,boxy1,lisy1,myi);

    dxy2+=psi[nx2]+psi[ny2];
  }
  dxy2/=N;*mir=psi[N]+phi[K]-dxy2;

  free(xc);free(nn);free(lis);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(ind);free(indx);free(indy);
  if (dimx==1) {free(mxi);free(boxx1);free(lisx1);} 
  else { for (i=0;i<BOX;i++) free(boxx[i]); free(boxx); free(lisx); for (d=0;d<dimx;d++) free(xx[d]); free(xx); }
  if (dimy==1) {free(myi);free(boxy1);free(lisy1);} 
  else { for (i=0;i<BOX;i++) free(boxy[i]); free(boxy); free(lisy); for (d=0;d<dimy;d++) free(yy[d]); free(yy); }
  free(phi);
}

void redK(double **x, int dim, int N, int K, 
	  double *psi,
	  double *scal,
	  double *mi_cr) {
  
  int i,k,ik,**n1,**n2;
  double *xc,dx;
  double **eps,*Eps;
  int *maxdim;
  double *dxy1,*dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;
  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 

  nn=(int*)calloc(K+1,sizeof(int));

  xc=(double*)calloc(dim+1,sizeof(double));

  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;
  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double**)calloc(K,sizeof(double*));
  Eps=(double*)calloc(K,sizeof(double));
  n1=(int**)calloc(K,sizeof(int*));
  n2=(int**)calloc(K,sizeof(int*));
  maxdim=(int*)calloc(K,sizeof(int));

  for (ik=0;ik<K;ik++) {
    eps[ik]=(double*)calloc(dim,sizeof(double));
    n1[ik]=(int*)calloc(dim,sizeof(int));
    n2[ik]=(int*)calloc(dim,sizeof(int));
    mi_cr[ik*2]=0;
    mi_cr[ik*2+1]=0;
  }

  dxy1=(double*)calloc(K,sizeof(double));
  dxy2=(double*)calloc(K,sizeof(double));

  make_box2(x,dim,N,0,dim-1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);
  
  for (ik=0;ik<K;ik++) {dxy1[ik]=dxy2[ik]=0.0;}
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
   
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    for (ik=0;ik<K;ik++) {
      Eps[ik]=0;maxdim[ik]=-1;
      for (d=0;d<dim;d++) {
	eps[ik][d]=0;
	for(k=1;k<=ik+1;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[ik][d] ) eps[ik][d]=dx; }
	if (eps[ik][d]>Eps[ik]) {Eps[ik]=eps[ik][d];maxdim[ik]=d;}
      }	
    }
    
    for (ik=0;ik<K;ik++) {
      for (d=0;d<dim;d++) {
	n2[ik][d]=neiE1(x[d],i,scal[d],BOX1,eps[ik][d],box1[d],lis1[d],mxi[d]);
	if (d==maxdim[ik]) { n1[ik][d]=n2[ik][d]; dxy1[ik]+=psi[n1[ik][d]]; }
	else { 
	  n1[ik][d]=neiE1(x[d],i,scal[d],BOX1,Eps[ik],box1[d],lis1[d],mxi[d]); 
	  dxy1[ik]+=psi[n1[ik][d]+1];
	}
	dxy2[ik]+=psi[n2[ik][d]];
      }
    }
  }
  for (ik=0;ik<K;ik++) {
    dxy1[ik]/=N;mi_cr[ik*2]=(dim-1)*psi[N]+psi[ik+1]-dxy1[ik];
    dxy2[ik]/=N;mi_cr[ik*2+1]=(dim-1)*psi[N]+phi[ik+1]-dxy2[ik];
  }

  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  
  for (ik=0;ik<K;ik++) {
    free(eps[ik]);
    free(n1[ik]);
    free(n2[ik]);
  }
  free(eps);free(n1);free(n2);
  free(Eps);free(maxdim);

  free(dxy1); free(dxy2);
  free(phi);
}


void mi2K(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mi_cr) {
  
  int i,k,ik,**n1,**n2;
  double *xc,dx;
  double **eps,*Eps;
  int *maxdim;
  double *dxy1,*dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  int dim=2;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 


  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double**)calloc(K,sizeof(double*));
  Eps=(double*)calloc(K,sizeof(double));
  n1=(int**)calloc(K,sizeof(int*));
  n2=(int**)calloc(K,sizeof(int*));
  maxdim=(int*)calloc(K,sizeof(int));

  for (ik=0;ik<K;ik++) {
    eps[ik]=(double*)calloc(dim,sizeof(double));
    n1[ik]=(int*)calloc(dim,sizeof(int));
    n2[ik]=(int*)calloc(dim,sizeof(int));
    mi_cr[ik*2]=0;
    mi_cr[ik*2+1]=0;
  }

  dxy1=(double*)calloc(K,sizeof(double));
  dxy2=(double*)calloc(K,sizeof(double));

  make_box2(x,dim,N,0,1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  for (ik=0;ik<K;ik++) {dxy1[ik]=dxy2[ik]=0.0;}
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    for (ik=0;ik<K;ik++) {
      Eps[ik]=0;maxdim[ik]=-1;
      for (d=0;d<dim;d++) {
	eps[ik][d]=0;
	for(k=1;k<=ik+1;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[ik][d] ) eps[ik][d]=dx; }
	if (eps[ik][d]>Eps[ik]) {Eps[ik]=eps[ik][d];maxdim[ik]=d;}
      }	
    }
    for (ik=0;ik<K;ik++) {
      for (d=0;d<dim;d++) {
	n2[ik][d]=neiE1(x[d],i,scal[d],BOX1,eps[ik][d],box1[d],lis1[d],mxi[d]);
	if (d==maxdim[ik]) { n1[ik][d]=n2[ik][d]; dxy1[ik]+=psi[n1[ik][d]]; }
	else { 
	  n1[ik][d]=neiE1(x[d],i,scal[d],BOX1,Eps[ik],box1[d],lis1[d],mxi[d]); 
	  dxy1[ik]+=psi[n1[ik][d]+1];
	}
	dxy2[ik]+=psi[n2[ik][d]];
      }
    }
  }
  for (ik=0;ik<K;ik++) {
    dxy1[ik]/=N;mi_cr[ik*2]=psi[N]+psi[ik+1]-dxy1[ik];
    dxy2[ik]/=N;mi_cr[ik*2+1]=psi[N]+phi[ik+1]-dxy2[ik];
  }

  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  
  for (ik=0;ik<K;ik++) {
    free(eps[ik]);
    free(n1[ik]);
    free(n2[ik]);
  }
  free(eps);free(n1);free(n2);
  free(Eps);free(maxdim);

  free(dxy1); free(dxy2);
  free(phi);
}

void mi_xnynK(double **x, int dimx, int dimy, int N, int K, 
	      double *psi, 
	      double *scal,
	      double *mi_cr) {

  int i,ik,k,*nx1,*ny1,*nx2,*ny2;
  double *xc,dy,dx;
  double *epsx,*epsy,*Eps;
  double *dxy1,*dxy2;
  double **xx,**yy;;
  double scalx, scaly;
  int *nn;
  int d;

  int BOX,BOX1;
  int **box,**boxy,**boxx,*lis; // two dimensional boxes
  int *lisy,*lisx; // lists for two dimensions
  int *boxx1, *boxy1; // onedimensional boxes
  int *lisy1,*lisx1; // lists for one dimensions
  int *mxi, *myi; //accumulative lists of points in oned boxes
  int *ind,*indx,*indy; //indeces of original data (the data resorted during box creating)
  double epsilon;
  int inveps;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-1/(double(i));   // 

  nn=(int*)calloc(K+1,sizeof(int));

  xc=(double*)calloc(dimx+dimy,sizeof(double));
  dxy1=(double*)calloc(K,sizeof(double));
  dxy2=(double*)calloc(K,sizeof(double));
  nx1=(int*)calloc(K,sizeof(int));
  nx2=(int*)calloc(K,sizeof(int));
  ny1=(int*)calloc(K,sizeof(int));
  ny2=(int*)calloc(K,sizeof(int));
  epsx=(double*)calloc(K,sizeof(double));
  epsy=(double*)calloc(K,sizeof(double));
  Eps=(double*)calloc(K,sizeof(double));

  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;
  BOX1=N-5;

  if (dimx>1) {
    xx=(double**)calloc(dimx,sizeof(double*));
    for (d=0;d<dimx;d++) xx[d]=(double*)calloc(N,sizeof(double));
    boxx=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxx[i]=(int*)calloc(BOX,sizeof(int));
    lisx=(int*)calloc(N,sizeof(int));
  } else { boxx1=(int*)calloc(BOX1+1,sizeof(int)); 
  lisx1=(int*)calloc(N,sizeof(int)); mxi=(int*)calloc(BOX1+1,sizeof(int)); }
  if (dimy>1) {
    yy=(double**)calloc(dimy,sizeof(double*));
    for (d=0;d<dimy;d++) yy[d]=(double*)calloc(N,sizeof(double));
    boxy=(int**)calloc(BOX,sizeof(int*)); 
    for (i=0;i<BOX;i++) boxy[i]=(int*)calloc(BOX,sizeof(int));
    lisy=(int*)calloc(N,sizeof(int));
  } else { boxy1=(int*)calloc(BOX1+1,sizeof(int)); 
  lisy1=(int*)calloc(N,sizeof(int)); myi=(int*)calloc(BOX1+1,sizeof(int)); }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  ind=(int*)calloc(N,sizeof(int));
  indx=(int*)calloc(N,sizeof(int));
  indy=(int*)calloc(N,sizeof(int));

  //save x if it would be reordered
  if (dimx>1) for (d=0;d<dimx;d++) memcpy(xx[d],x[d],N*sizeof(double)); 
  if (dimy>1) for (d=0;d<dimy;d++) memcpy(yy[d],x[d+dimx],N*sizeof(double));

  make_box2ind(x,dimx+dimy,N,0,dimx,BOX,inveps,ind,box,lis); 
  //for searching neighbours in product space

  if (dimx==1) {scalx=scal[0]; make_box1(x[0],N,scalx,BOX1,boxx1,lisx1,mxi);}
  else make_box2ind(xx,dimx,N,0,dimx-1,BOX,inveps,indx,boxx,lisx); 
  if (dimy==1) {scaly=scal[dimx]; make_box1(x[dimx],N,scaly,BOX1,boxy1,lisy1,myi); }
  else make_box2ind(yy,dimy,N,0,dimy-1,BOX,inveps,indy,boxy,lisy); 

  for (ik=0;ik<K;ik++) {dxy1[ik]=dxy2[ik]=0.0;}
  for (i=0;i<N;i++) {
    for (d=0;d<dimx+dimy;d++) xc[d]=x[d][ind[i]];

    neiK(x,dimx+dimy,0,dimx,ind[i],BOX,epsilon,K,box,lis,nn);
    for (ik=0;ik<K;ik++) {
      epsx[ik]=0; for (d=0;d<dimx;d++) 
	for(k=1;k<=ik+1;k++) if( (dx=fabs(xc[d]-x[d][nn[k]]))>epsx[ik] ) epsx[ik]=dx;
      epsy[ik]=0; for (d=dimx;d<dimx+dimy;d++) 
	for(k=1;k<=ik+1;k++) if( (dy=fabs(xc[d]-x[d][nn[k]]))>epsy[ik] ) epsy[ik]=dy;
      if (dimx>1) nx2[ik]=neiE(xx,indx[i],0,dimx-1,dimx,BOX,epsilon,epsx[ik],boxx,lisx);
      else  nx2[ik]=neiE1(x[0],ind[i],scalx,BOX1,epsx[ik],boxx1,lisx1,mxi);
      if (dimy>1) ny2[ik]=neiE(yy,indy[i],0,dimy-1,dimy,BOX,epsilon,epsy[ik],boxy,lisy);
      else ny2[ik]=neiE1(x[dimx],ind[i],scaly,BOX1,epsy[ik],boxy1,lisy1,myi);
    
      if (epsx[ik]>epsy[ik]) {
	Eps[ik]=epsx[ik];nx1[ik]=nx2[ik];
	if (dimy>1) ny1[ik]=neiE(yy,indy[i],0,dimy-1,dimy,BOX,epsilon,Eps[ik],boxy,lisy);
	else ny1[ik]=neiE1(x[dimx],ind[i],scaly,BOX1,Eps[ik],boxy1,lisy1,myi);
	dxy1[ik]+=psi[nx1[ik]]+psi[ny1[ik]+1];
      } else {
	Eps[ik]=epsy[ik];ny1[ik]=ny2[ik];
	if (dimx>1) nx1[ik]=neiE(xx,indx[i],0,dimx-1,dimx,BOX,epsilon,Eps[ik],boxx,lisx);
	else nx1[ik]=neiE1(x[0],ind[i],scalx,BOX1,Eps[ik],boxx1,lisx1,mxi);
	dxy1[ik]+=psi[nx1[ik]+1]+psi[ny1[ik]];
      }
      dxy2[ik]+=psi[nx2[ik]]+psi[ny2[ik]];
    }
  }
  for (ik=0;ik<K;ik++) {
    dxy1[ik]/=N;mi_cr[2*ik]=psi[N]+psi[ik+1]-dxy1[ik];
    dxy2[ik]/=N;mi_cr[2*ik+1]=psi[N]+phi[ik+1]-dxy2[ik];
  }

  free(xc);free(nn);free(lis);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(ind);free(indx);free(indy);
  if (dimx==1) {free(mxi);free(boxx1);free(lisx1);} 
  else { for (i=0;i<BOX;i++) free(boxx[i]); free(boxx); free(lisx); for (d=0;d<dimx;d++) free(xx[d]); free(xx); }
  if (dimy==1) {free(myi);free(boxy1);free(lisy1);} 
  else { for (i=0;i<BOX;i++) free(boxy[i]); free(boxy); free(lisy); for (d=0;d<dimy;d++) free(yy[d]); free(yy); }
  free(phi);
  free(epsx);free(epsy);free(Eps);
  free(nx1);free(nx2);free(ny1);free(ny2);
  free(dxy1);free(dxy2);
}

void mi_xnynKembed(double **x, int dim, int N, int K, 
		   double **xx, double **yy, 
		   int **boxx, int *lisx, int *indx,
		   int **boxy, int *lisy, int *indy,
		   double *psi, 
		   double *mi_cr) {

  int i,ik,k,*nx1,*ny1,*nx2,*ny2;
  double *xc,dy,dx;
  double *epsx,*epsy,*Eps;
  char *maxdim;
  double *dxy1,*dxy2;
  int *nn;
  int d;

  int BOX;
  int **box,*lis; // two dimensional boxes
  int *ind; //indeces of original data (the data resorted during box creating)
  double epsilon;
  int inveps;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-1/(double(i));   // 

  nn=(int*)calloc(K+1,sizeof(int));

  xc=(double*)calloc(2*dim,sizeof(double));
  dxy1=(double*)calloc(K,sizeof(double));
  dxy2=(double*)calloc(K,sizeof(double));
  nx1=(int*)calloc(K,sizeof(int));
  nx2=(int*)calloc(K,sizeof(int));
  ny1=(int*)calloc(K,sizeof(int));
  ny2=(int*)calloc(K,sizeof(int));
  epsx=(double*)calloc(K,sizeof(double));
  epsy=(double*)calloc(K,sizeof(double));
  Eps=(double*)calloc(K,sizeof(double));
  maxdim=(char*)calloc(K,sizeof(char));

  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  ind=(int*)calloc(N,sizeof(int));


  make_box2ind(x,2*dim,N,0,dim,BOX,inveps,ind,box,lis); 

  for (ik=0;ik<K;ik++) {dxy1[ik]=dxy2[ik]=0.0;}
  for (i=0;i<N;i++) {
    for (d=0;d<2*dim;d++) xc[d]=x[d][ind[i]];

    neiK(x,2*dim,0,dim,ind[i],BOX,epsilon,K,box,lis,nn);
    for (ik=0;ik<K;ik++) {
      epsx[ik]=0; for (d=0;d<dim;d++) 
	for(k=1;k<=ik+1;k++) if( (dx=fabs(xc[d]-x[d][nn[k]]))>epsx[ik] ) epsx[ik]=dx;
      epsy[ik]=0; for (d=dim;d<2*dim;d++) 
	for(k=1;k<=ik+1;k++) if( (dy=fabs(xc[d]-x[d][nn[k]]))>epsy[ik] ) epsy[ik]=dy;
    }
    neiEK(xx,indx[i],0,dim-1,dim,K,BOX,epsilon,epsx,boxx,lisx,nx2);
    neiEK(yy,indy[i],0,dim-1,dim,K,BOX,epsilon,epsy,boxy,lisy,ny2);
    for (ik=0;ik<K;ik++) {
      maxdim[ik]=0;
      if (epsx[ik]>epsy[ik]) { Eps[ik]=epsx[ik]; } else { Eps[ik]=epsy[ik]; maxdim[ik]=1; }
    }
    neiEK(yy,indy[i],0,dim-1,dim,K,BOX,epsilon,Eps,boxy,lisy,ny1);
    neiEK(xx,indx[i],0,dim-1,dim,K,BOX,epsilon,Eps,boxx,lisx,nx1);

    for (ik=0;ik<K;ik++) {
      if (maxdim[ik]) dxy1[ik]+=psi[nx1[ik]+1]+psi[ny1[ik]]; else dxy1[ik]+=psi[nx1[ik]]+psi[ny1[ik]+1];
      dxy2[ik]+=psi[nx2[ik]]+psi[ny2[ik]];
    }
   
  }
  for (ik=0;ik<K;ik++) {
    dxy1[ik]/=N;mi_cr[2*ik]=psi[N]+psi[ik+1]-dxy1[ik];
    dxy2[ik]/=N;mi_cr[2*ik+1]=psi[N]+phi[ik+1]-dxy2[ik];
  }

  free(xc);free(nn);free(lis);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(ind);
  free(phi);
  free(epsx);free(epsy);free(Eps);free(maxdim);
  free(nx1);free(nx2);free(ny1);free(ny2);
  free(dxy1);free(dxy2);
}



void mi2h(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mic, double *mir, double *hc, double *hr) {
  
  int i,k,*n1,*n2;
  double *xc,dx;
  double *eps,Eps;
  double *h1,*h2,hj1,hj2;
  int maxdim;
  double dxy1,dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  int dim=2;

  double *phi;
  double t_d1,t_d2;

  phi=(double*)calloc(N+1,sizeof(double));
  for (i=1;i<=N;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 


  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n1=(int*)calloc(dim,sizeof(int));
  n2=(int*)calloc(dim,sizeof(int));

  h1=(double*)calloc(dim,sizeof(double));
  h2=(double*)calloc(dim,sizeof(double));
  for (d=0;d<dim;d++) h1[d]=h2[d]=0;
  hj1=hj2=0;

  make_box2(x,dim,N,0,1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  dxy1=dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    Eps=0;maxdim=-1;
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for(k=1;k<=K;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[d] ) eps[d]=dx; }
      if (eps[d]>Eps) {Eps=eps[d];maxdim=d;}
    }	
    for (d=0;d<dim;d++) {
      n2[d]=neiE1(x[d],i,scal[d],BOX1,eps[d],box1[d],lis1[d],mxi[d]);
      if (d==maxdim) { 
	n1[d]=n2[d]; dxy1+=psi[n1[d]];
	h1[d]+=log(Eps)-psi[n1[d]];
      }
      else { 
	n1[d]=neiE1(x[d],i,scal[d],BOX1,Eps,box1[d],lis1[d],mxi[d]); 
	dxy1+=psi[n1[d]+1];
	h1[d]+=log(Eps)-psi[n1[d]+1];
      }
      dxy2+=psi[n2[d]];

      h2[d]+=log(eps[d])-phi[n2[d]];
      hj1+=log(Eps);
      hj2+=log(eps[d]);
    }
    //    fprintf(stdout,"%f %f %d %d %d %d %d %f %d\n",xc[0],xc[1],(eps[0]<=eps[1]),n1[0],n1[1],n2[0],n2[1],Eps,nn[K]);
    
  }
  dxy1/=N;*mic=psi[N]+psi[K]-dxy1;
  dxy2/=N;*mir=psi[N]+phi[K]-dxy2;
  hj1/=N; *hc=psi[N]-psi[K]+hj1;//dim-1
  hj2/=N; *hr=psi[N]-phi[K]+hj2;//dim-1
  t_d1=+1e30;
  t_d2=+1e30;
  fprintf(stdout,  "hj_c\t%1.3f\t\t\t\thj_r\t%1.3f\n",*hc,*hr);
  for (d=0;d<dim;d++) {
    h1[d]/=N;h1[d]+=psi[N];
    h2[d]/=N;h2[d]+=psi[N];
    if (h1[d]<t_d1) t_d1=h1[d];
    if (h2[d]<t_d2) t_d2=h2[d];
    if ((*hc-h1[d])<t_d1) t_d1=*hc-h1[d];
    if ((*hr-h2[d])<t_d2) t_d2=*hr-h2[d];
    fprintf(stdout,"h[%d]_c\t%1.3f\thcond[%d]_c\t%1.3f\th[%d]_r\t%1.3f\thcond[%d]_r\t%1.3f\n",
	    d,h1[d],d,*hc-h1[d],d,h2[d],d,*hr-h2[d]);
  }
  *hc-=2*t_d1;
  *hr-=2*t_d2;
  
  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  free(eps);free(n1);free(n2);
  free(phi);
  free(h1);free(h2);
}










void mi_d(double **x, int dim, int N, int K, float *mi_cr,
	  double *psi, double *phi, double minx, double maxx, double miny, double maxy) {
  
  int i,k,ik,**n1,**n2;
  double *xc,dx;
  double **eps,*Eps;
  int *maxdim;
  double *dxy1,*dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  double scal[2];

  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX*K<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;
  scal[0]=BOX1/(maxx-minx);
  scal[1]=BOX1/(maxy-miny);

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double**)calloc(K,sizeof(double*));
  Eps=(double*)calloc(K,sizeof(double));
  n1=(int**)calloc(K,sizeof(int*));
  n2=(int**)calloc(K,sizeof(int*));
  maxdim=(int*)calloc(K,sizeof(int));

  for (ik=0;ik<K;ik++) {
    eps[ik]=(double*)calloc(dim,sizeof(double));
    n1[ik]=(int*)calloc(dim,sizeof(int));
    n2[ik]=(int*)calloc(dim,sizeof(int));
    mi_cr[ik*2]=0;
    mi_cr[ik*2+1]=0;
  }

  dxy1=(double*)calloc(K,sizeof(double));
  dxy2=(double*)calloc(K,sizeof(double));

  make_box2(x,dim,N,0,1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) {
    make_box1(x[d],N,scal[d&1],BOX1,box1[d],lis1[d],mxi[d]);      
  }
  
  for (ik=0;ik<K;ik++) {dxy1[ik]=dxy2[ik]=0.0;}
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,1,i,BOX,epsilon,K,box,lis,nn);
    
    for (ik=0;ik<K;ik++) {
      Eps[ik]=0;maxdim[ik]=-1;
      for (d=0;d<dim;d++) {
	eps[ik][d]=0;
	for(k=1;k<=ik+1;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[ik][d] ) eps[ik][d]=dx; }
	if (eps[ik][d]>Eps[ik]) {Eps[ik]=eps[ik][d];maxdim[ik]=d;}
      }	
    }
    for (ik=0;ik<K;ik++) {
      for (d=0;d<dim;d++) {
	n2[ik][d]=neiE1(x[d],i,scal[d&1],BOX1,eps[ik][d],box1[d],lis1[d],mxi[d]);
	if (d==maxdim[ik]) { n1[ik][d]=n2[ik][d]; dxy1[ik]+=psi[n1[ik][d]]; }
	else { 
	  n1[ik][d]=neiE1(x[d],i,scal[d&1],BOX1,Eps[ik],box1[d],lis1[d],mxi[d]); 
	  dxy1[ik]+=psi[n1[ik][d]+1];
	}
	dxy2[ik]+=psi[n2[ik][d]];
      }
    }
  }
  for (ik=0;ik<K;ik++) {
    dxy1[ik]/=N;mi_cr[ik*2]=float((dim-1)*psi[N]+psi[ik+1]-dxy1[ik]);
    dxy2[ik]/=N;mi_cr[ik*2+1]=float((dim-1)*psi[N]+phi[ik+1]-dxy2[ik]);
  }

  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  
  for (ik=0;ik<K;ik++) {
    free(eps[ik]);
    free(n1[ik]);
    free(n2[ik]);
  }
  free(eps);free(n1);free(n2);
  free(Eps);free(maxdim);

  free(dxy1); free(dxy2);
}




void mi2r_(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mir) {
  
  int i,k,*n2;
  double *xc,dx;
  double *eps;
  double dxy2;
  int *nn;
  int d;
  
  int BOX,BOX1;
  int **box,*lis; // two dimensional boxes
  int **box1; // onedimensional boxes
  int **lis1; // lists for one dimensions
  int **mxi; //accumulative lists of points in oned boxes
  double epsilon;
  int inveps;

  int dim=2;

  double *phi;

  phi=(double*)calloc(K+1,sizeof(double));
  for (i=1;i<=K;i++) phi[i]=psi[i]-(dim-1)/(double(i));   // 


  nn=(int*)calloc(K+1,sizeof(int));
  xc=(double*)calloc(dim+1,sizeof(double));
  BOX=1; while (0.5*BOX*BOX<N) BOX*=2;
  epsilon=4.0/BOX;
  inveps=BOX/4;

  BOX1=N-5;

  box1=(int**)calloc(dim,sizeof(int*)); 
  lis1=(int**)calloc(dim,sizeof(int*)); 
  mxi=(int**)calloc(dim,sizeof(int*)); 
  for (d=0;d<dim;d++) {
    box1[d]=(int*)calloc(BOX1+1,sizeof(int)); 
    lis1[d]=(int*)calloc(N,sizeof(int)); 
    mxi[d]=(int*)calloc(BOX1+1,sizeof(int));
  }
 
  box=(int**)calloc(BOX,sizeof(int*));
  for (i=0;i<BOX;i++) box[i]=(int*)calloc(BOX,sizeof(int));
  lis=(int*)calloc(N,sizeof(int));

  eps=(double*)calloc(dim,sizeof(double));
  n2=(int*)calloc(dim,sizeof(int));

  make_box2(x,dim,N,0,1,BOX,inveps,box,lis); //for searching neighbours in prodict space
  for (d=0;d<dim;d++) make_box1(x[d],N,scal[d],BOX1,box1[d],lis1[d],mxi[d]);      
  
  dxy2=0.0;
  for (i=0;i<N;i++) {
    for (d=0;d<dim;d++) xc[d]=x[d][i];
    neiK(x,dim,0,dim-1,i,BOX,epsilon,K,box,lis,nn);
    
    for (d=0;d<dim;d++) {
      eps[d]=0;
      for(k=1;k<=K;k++) {if( (dx=fabs(xc[d]-x[d][nn[k]]))>eps[d] ) eps[d]=dx; }
    }	
    for (d=0;d<dim;d++) {
      n2[d]=neiE1(x[d],i,scal[d],BOX1,eps[d],box1[d],lis1[d],mxi[d]);
      dxy2+=psi[n2[d]];
    }
    //    fprintf(stdout,"%f %f %d %d %d %d\n",xc[0],xc[1],(eps[0]<=eps[1]),n2[0],n2[1],nn[K]);
  }
  dxy2/=N;*mir=psi[N]+phi[K]-dxy2;
  
  free(xc);free(nn);
  for (i=0;i<BOX;i++) free(box[i]); free(box);
  free(lis);
  for (d=0;d<dim;d++) {
    free(box1[d]);free(lis1[d]);free(mxi[d]);
  }
  free(box1);free(lis1);free(mxi);
  free(eps);free(n2);
  free(phi);
}
