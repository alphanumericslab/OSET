void make_box1(double *x, int N, double scal, int bs, 
	       int *box, int *lis, int *mxi);
/* make one-dimensional box
input
   x - input sequence
   N - length of the input sequence
   scal - renormalization factor, scal=#boxes/(max(x)-min(x))
   bs - #boxes
output
   box - each element of the array is ordinary number of the last point in the box
   lis - each element of the array is ordinary number of the previous points in the box or -1
   mxi - accumulative number of the points in the box
 */
void make_box2(double **x, int dim, int N, int comp1, int comp2, int bs, int inveps, 
	       int **box, int *lis);
/* uses comp1 and comp2 component for two-dimensional grid */
/* make two-dimensional box !!! changing the order of the points
input
   x - two-dimensional input sequence
   dim - dimension of x
   N - length of the input sequence
   comp1 - first component for two dimensional grid
   comp2 - second component for two dimensional grid
   bs - #boxes
   inveps - box size
output
   box - each element of the array is ordinary number of the last point in the box
   lis - each element of the array is ordinary number of the previous points in the box or -1
 */
void make_box2ind(double **x, int dim, int N, int comp1, int comp2, int bs, int inveps, 
		  int *ind,  int **box, int *lis);
/* make two-dimensional box !!! changing the order of the points, but saving the order in ind
input
   x - dim-dimensional input sequence
   dim - dimension of input sequence 
   N - length of the input sequence
   comp1 - first component for creating grid
   comp2 - first component for creating grid
   bs - #boxes
   inveps - box size
output
   ind - index of original data 
   box - each element of the array is ordinary number of the last point in the box
   lis - each element of the array is ordinary number of the previous points in the box or -1
 */
int neiE1(double *x, int i, double scal, int bs, double eps,  int *box, int *lis, int *mxi);
/* searching for neighbors of point i in eps-neighborhood in one dimension
input
   x - one-dimensional input sequence
   i - current point
   scal - renormalization factor, scal=#boxes/(max(x)-min(x))
   bs - #boxes
   eps - neighborhood
   box - each element of the array is ordinary number of the last point in the box
   lis - each element of the array is ordinary number of the previous points in the box or -1
   mxi - accumulative number of the points in the box
output
   number of neighbors
 */
int neiE(double **x, int i, int comp1, int comp2, int dim, int bs, double epsgrid, double eps, int **box, int *lis);
/* searching for neighbors of point i in eps-neighborhood in dim dimension
input
   x - dim-dimensional input sequence
   i - current point
   dim -dimension on input sequence 
   comp1 - 
   comp2 -
   bs - #boxes
   epsgrid - size of the grid
   eps - neighborhood
   box - each element of the array is ordinary number of the last point in the box
   lis - each element of the array is ordinary number of the previous points in the box or -1
output
   number of neighbors
 */
void neiEK(double **x, int i, int comp1, int comp2, int dim, int K,
	   int bs, double epsgrid, double *eps, int **box, int *lis,
	   int *nx); 
/* searching for neighbors of point i in eps-neighborhood in dim dimension for K
input
   x - dim-dimensional input sequence
   i - current point
   dim -dimension on input sequence 
   K - max number of neighbors
   comp1 - 
   comp2 -
   bs - #boxes
   epsgrid - size of the grid
   eps - neighborhood
   box - each element of the array is ordinary number of the last point in the box
   lis - each element of the array is ordinary number of the previous points in the box or -1
output
   nx - number of neighbors
 */


void neiK(double **x, int dim, int comp1, int comp2, int i, 
	  int bs, double epsgrid, int K, int **box, int *lis,
	  int *nn);
/* searching for K neighbors of point i in dim dimension
input
   x - dim-dimensional input sequence
   dim - dimension of x
   comp1 - first component for two dimensional grid
   comp2 - second component for two dimensional grid
   i - current point
   bs - #boxes
   epsgrid - size of the grid
   box - each element of the array is ordinary number of the last point in the box
   lis - each element of the array is ordinary number of the previous points in the box or -1
output
   nn - indices of K neighbors
 */

  //void neiK(double **x, int i, int bs, double epsgrid, int K, int **box, int *lis,
  //	  int *nn);
/* searching for K neighbours of point i in two dimension
input
   x - 2-dimensional input sequence
   i - current point
   bs - #boxes
   epsgrid - size of the grid
   box - each element of the array is ordinary number of the last point in the box
   lis - each element of the array is ordinary number of the previous points in the box or -1
output
   nn - indices of K neighbors
 */
void mi2(double **x, int N, int K,
	 double *psi, 
	 double *scal,
	 double *mic, double *mir);
void mi2c(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mic);
void mi2r(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mir);
void mi2r_(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mir);

void red(double **x, int dim, int N, int K, 
	 double *psi,
	 double *scal,
	 double *mic, double *mir);
void redc(double **x, int dim, int N, int K, 
	  double *psi,
	  double *scal,
	  double *mic);
void redr(double **x, int dim, int N, int K, 
	  double *psi,
	  double *scal,
	  double *mir);
void redr_embed(double **x, int dim, int edim, int tau, int N, int K, 
		double *psi,
		double *mir);
void mi_xnyn(double **x, int dimx, int dimy, int N, int K, 
	     double *psi, 
	     double *scal,
	     double *mic, double *mir);
/*
calculating of mutual information between vectors, only for one K

input
   x -  (dimx+dimy)-dimensional input
   dimx - dimension of vector x
   dimy - dimension of vector y
   N - length of vectors
   K - max number of neighbors
   psi - digamma function
   scal - scale for one dimensional box (BOX1/(max-min))
output
   mic - cubic method
   mir - rectange method
*/
void mic_xnyn(double **x, int dimx, int dimy, int N, int K, 
	      double *psi, 
	      double *scal,
	      double *mic);
void mir_xnyn(double **x, int dimx, int dimy, int N, int K, 
	      double *psi, 
	      double *scal,
	      double *mir);



void redK(double **x, int dim, int N, int K, 
	  double *psi, 
	  double *scal,
	  double *mi_cr);
/*
calculating of redundancy 

input
   x -  dim-dimensional input
   dim - dimension
   N - length of vectors
   K - max number of neighbors
   psi - digamma function
   scal - scale for one dimensional box (BOX1/(max-min))
output
   mi_cr - MI using cubic and rectangle method for different K
      array contains the values of mutual information [cub_k1,rec_k1,cub_k2,rec_k2,...,cub_kK,rec_kK]
*/

void mi2K(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mi_cr);
/*
calculating of mutual information

input
   x -  two-dimensional input
   N - length of vectors
   K - max number of neighbors
   psi - digamma function
   scal - scale for one dimensional box (BOX1/(max-min))
output
   mi_cr - MI using cubic and rectangle method for different K
      array contains the values of mutual information [cub_k1,rec_k1,cub_k2,rec_k2,...,cub_kK,rec_kK]
*/

void mi_xnynK(double **x, int dimx, int dimy, int N, int K, 
	      double *psi, 
	      double *scal,
	      double *mi_cr);
void mi_xnynKembed(double **x, int dim, int N, int K, 
		   double **xx, double **yy, 
		   int **boxx, int *lisx, int *indx,
		   int **boxy, int *lisy, int *indy,
		   double *psi, 
		   double *mi_cr);


void mi2h(double **x, int N, int K,
	  double *psi, 
	  double *scal,
	  double *mic, double *mir, double *hc, double *hr);


void mi_embed(double **x, int dim, int N, int K, float *mi_cr,
	  double *psi, double *phi, double minx, double maxx, double miny, double maxy);


