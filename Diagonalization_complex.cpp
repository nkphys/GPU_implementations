#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <cuda.h>
#include <iomanip>
#include <complex.h>
#include "magma_v2.h"
#include "magma_lapack.h"

using namespace std;
int main (int argc ,char ** argv ) {
magma_init (); // initialize Magma
magma_queue_t queue = NULL ;

magma_int_t dev =0;
magma_queue_create(dev ,&queue);
double gpu_time , cpu_time ;
magma_int_t n=3000 , n2=n*n;
magmaDoubleComplex *a, *r; // a, r - nxn matrices on the host
magmaDoubleComplex *d_r ; // nxn matrix on the device
magmaDoubleComplex *h_work ; // workspace
double *rwork;
magma_int_t lwork ; // h_work size
magma_int_t *iwork ; // workspace
magma_int_t liwork ; // iwork size
magma_int_t lrwork; //rwork size
double *w1 , *w2; // w1 ,w2 - vectors of eigenvalues
double error , work[1]; // used in difference computations

magma_int_t ione = 1 , info ;
double mione = -1.0;
magma_int_t incr = 1;
magma_int_t ISEED[4] = {0 ,0 ,0 ,1}; // seed
magma_dmalloc_cpu(&w1 ,n); // host memory for real
magma_dmalloc_cpu(&w2 ,n); // eigenvalues
magma_zmalloc_cpu(&a,n2 ); // host memory for a
magma_zmalloc_cpu(&r,n2 ); // host memory for r
magma_zmalloc(&d_r ,n2 ); // device memory for d_r

magmaDoubleComplex aux_work[1];
magma_int_t aux_iwork[1];
double aux_rwork[1];

magma_zheevd_gpu(MagmaVec , MagmaLower ,n,d_r ,n,w1 ,r,n, aux_work ,-1, aux_rwork, -1,  aux_iwork , -1 ,& info );
lrwork= (magma_int_t)aux_rwork[0];
lwork = (magma_int_t)(MAGMA_Z_REAL(aux_work[0]));

//cout<<(magma_int_t)(MAGMA_Z_REAL(aux_work[0]))<<"   "<<(magma_int_t)(MAGMA_Z_IMAG(aux_work[0]))<<endl;

liwork = aux_iwork[0];
iwork =( magma_int_t *) malloc ( liwork * sizeof ( magma_int_t ));
magma_zmalloc_cpu(&h_work , lwork ); // memory for workspace
magma_dmalloc_cpu(&rwork , lrwork ); // memory for workspace

cout<<liwork<< "(int), "<<lwork<<"(complex), "<<lrwork<<"(double) "<<endl;
// Randomize the matrix a and copy a -> r
lapackf77_zlarnv(&ione ,ISEED ,&n2 ,a);
for(int i=0;i<n;i++){
a[i+i*n]=MAGMA_Z_MAKE(0.0,0.0);
}


// 1d NN-tight binding model
/*
for(int i=0;i<n*n;i++){
a[i]=MAGMA_Z_MAKE(0.0,0.0);
}
for(int i=1;i<n-1;i++){
a[i+(i-1)*n]=MAGMA_Z_MAKE(1.0,0.0);
a[i+(i+1)*n]=MAGMA_Z_MAKE(1.0,0.0);
}
a[0+(n-1)*n]=MAGMA_Z_MAKE(1.0,0.0);a[0+(1)*n]=MAGMA_Z_MAKE(1.0,0.0);
a[n-1+(0)*n]=MAGMA_Z_MAKE(1.0,0.0);a[n-1+(n-2)*n]=MAGMA_Z_MAKE(1.0,0.0);
*/

/*
complex<double> temp_z;
printf("a: \n");
  for (int i = 0; i < 4; i++){
for(int j=0;j<4;j++){
//printf("%2.1g + i%2.1g ", a[i+j*4]);
temp_z.real(MAGMA_Z_REAL(a[i+j*4]));
temp_z.imag(MAGMA_Z_IMAG(a[i+j*4]));
cout<<temp_z<<" ";
}
printf("\n");
}
printf("\n");
*/



lapackf77_zlacpy(MagmaFullStr ,&n ,&n,a ,&n,r ,&n);
magma_zsetmatrix( n, n, a, n, d_r ,n, queue ); // copy a -> d_r
// compute the eigenvalues and eigenvectors for a symmetric ,
// real nxn matrix ; Magma version
 gpu_time = magma_sync_wtime ( NULL );

magma_zheevd_gpu(MagmaVec,MagmaLower,n,d_r,n,w1,r,n,h_work,lwork,rwork, lrwork, iwork,liwork,&info);
gpu_time = magma_sync_wtime ( NULL ) - gpu_time ;
printf (" zheevd_gpu gpu time : %7.5f sec .\n",gpu_time ); // Magma time

cout<<setprecision(10)<<endl;
for(int i=0;i<8;i++){
cout<<w1[i]<<"  ";
}
cout<<"......";
for(int i=0;i<8;i++){
cout<<w1[n+i-8]<<"  ";
}
cout<<endl;



// Lapack version
cpu_time = magma_wtime ();
lapackf77_zheevd ("V","L" ,&n,a ,&n,w2 ,h_work ,&lwork , rwork, &lrwork,  iwork ,&liwork ,&info );
cpu_time = magma_wtime () - cpu_time ;
printf (" zheevd cpu time : %7.5f sec .\n",cpu_time ); // Lapack
// difference in eigenvalues

for(int i=0;i<8;i++){
cout<<w2[i]<<"  ";
}
cout<<"......";
for(int i=0;i<8;i++){
cout<<w2[n+i-8]<<"  ";
}
cout<<endl;



blasf77_daxpy ( &n, &mione , w1 , &incr , w2 , &incr );
error = lapackf77_dlange ( "M", &n, &ione , w2 , &n, work );

printf (" difference in eigenvalues : %e\n",error );
free (w1 ); // free host memory
free (w2 ); // free host memory
free (a); // free host memory
free (r); // free host memory
free ( h_work ); // free host memory
magma_free (d_r ); // free device memory
magma_queue_destroy ( queue ); // destroy queue
magma_finalize (); // finalize Magma
return EXIT_SUCCESS ;
}
