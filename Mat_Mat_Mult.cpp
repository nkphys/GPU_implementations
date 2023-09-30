#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "magma_v2.h"
#include "magma_lapack.h"

int main ( int argc , char ** argv ){
magma_init (); // initialize Magma
magma_queue_t queue = NULL ;
magma_int_t dev =0;
magma_queue_create (dev ,& queue );
real_Double_t gpu_time ;
magma_int_t m = 1280;//92; // a - mxk matrix
magma_int_t n = 1280;//96; // b - kxn matrix
magma_int_t k = 1280;//48; // c - mxn matrix
magma_int_t mk=m*k; // size of a
magma_int_t kn=k*n; // size of b
magma_int_t mn=m*n; // size of c
float *a; // a- mxk matrix
float *b; // b- kxn matrix
float *c; // c- mxn matrix
float alpha = MAGMA_S_MAKE ( 1.0 , 0.0 ); // alpha =1
float beta = MAGMA_S_MAKE ( 0.0 , 0.0 ); // beta =1
magma_int_t ione = 1; // random uniform distr . in (0 ,1)
magma_int_t ISEED [4] = { 0 ,1 ,2 ,3 }; // seed
cudaMallocManaged (&a,mk* sizeof ( float )); // unified mem. for a
cudaMallocManaged (&b,kn* sizeof ( float )); // unified mem. for b
cudaMallocManaged (&c,mn* sizeof ( float )); // unified mem. for c

// generate random matrices a, b, c;
lapackf77_slarnv (& ione ,ISEED ,&mk ,a); // randomize a
lapackf77_slarnv (& ione ,ISEED ,&kn ,b); // randomize b
lapackf77_slarnv (& ione ,ISEED ,&mn ,c); // randomize c
// matrix - matrix multiplication : c = al*a*b + bet *c
// a -mxk matrix , b -kxn matrix , c -mxn matrix ;
// al ,bet - scalars


gpu_time = magma_sync_wtime ( NULL );
magma_sgemm(MagmaNoTrans,MagmaNoTrans,m,n,k,alpha,a,m,b,k,beta,c,m,queue);

/*
for(int ra=0;ra<m;ra++){
for(int cb=0;cb<n;cb++){
c[ra + m*cb]=0.0;
for(int ca=0;ca<k;ca++){
c[ra + m*cb] += a[ra + m*ca]*b[ca + k*cb];
}
}
}
*/

gpu_time = magma_sync_wtime ( NULL ) - gpu_time ;
printf ("magma_sgemm time : %7.5f sec .\n",gpu_time );
printf ("after magma_sgemm :\n");
printf ("c:\n");
for(int i=0;i <4;i ++){
for(int j=0;j <4;j++) printf (" %10.4f,",c[i*m+j]);
printf (" ...\n");}
printf (" ...............................................\n");
magma_free (a); // free memory
magma_free (b); // free memory
magma_free (c); // free memory
magma_queue_destroy ( queue ); // destroy queue
magma_finalize (); // finalize Magma
return 0;
}

