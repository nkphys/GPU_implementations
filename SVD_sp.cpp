#include <stdio.h>
#include <cuda.h>
#include "magma_v2.h"
#include "magma_lapack.h"
#define min(a,b)  (((a)<(b))?(a):(b))
int main ( int argc , char ** argv )
{
magma_init (); // initialize Magma
real_Double_t gpu_time , cpu_time ;

magma_int_t m=9600 , n =9600 , n2=m*n, min_mn = min(m,n);
float *a, *r; // a,r - mxn matrices
float *u, *vt;// u - mxm matrix , vt - nxn matrix on the host
float *s1 , *s2; // vectors of singular values
magma_int_t info ;
magma_int_t ione = 1;
float work[1] , error = 1.; // used in difference computations
float mone = -1.0 , * h_work ; // h_work - workspace
magma_int_t lwork ; // workspace size
magma_int_t ISEED[4] = {0 ,0 ,0 ,1}; 

magma_smalloc_cpu(&a,n2 ); // host memory for a
magma_smalloc_cpu(&vt ,n*n); // host memory for vt
magma_smalloc_cpu(&u,m*m); // host memory for u
magma_smalloc_cpu(&s1 , min_mn ); // host memory for s1

magma_smalloc_cpu(&s2 , min_mn ); // host memory for s2
magma_smalloc_pinned(&r,n2 ); // host memory for r
magma_int_t nb = magma_get_sgesvd_nb(m,n); // optim . block size
lwork = min_mn*min_mn +2*min_mn +2*min_mn *nb;
magma_smalloc_pinned(&h_work , lwork ); // host mem . for h_work


// Randomize the matrix a
lapackf77_slarnv(&ione ,ISEED ,&n2 ,a);
lapackf77_slacpy(MagmaFullStr ,&m ,&n,a ,&m,r ,&m); //a- >r
// MAGMA
 gpu_time = magma_wtime();

magma_sgesvd(MagmaNoVec,MagmaNoVec,m,n,r,m,s1,u,m,vt,n,h_work,lwork,&info );
gpu_time = magma_wtime () - gpu_time ;
printf (" sgesvd gpu time : %7.5f\n", gpu_time ); // Magma time


cpu_time = magma_wtime ();
lapackf77_sgesvd("N","N" ,&m ,&n,a ,&m,s2 ,u ,&m,vt ,&n,h_work ,&lwork ,&info);
cpu_time = magma_wtime () - cpu_time ;
printf (" sgesvd cpu time : %7.5f\n", cpu_time ); // Lapack time


error = lapackf77_slange("f" ,&min_mn ,&ione ,s1 ,&min_mn , work );
blasf77_saxpy(&min_mn ,&mone ,s1 ,&ione ,s2 ,&ione );
error = lapackf77_slange("f" ,&min_mn ,&ione ,s2 ,&min_mn , work );

// error ;
 printf(" difference : %e\n", error); // difference in singul .


free(a); // free host memory
free(vt); // free host memory
free(s1); // free host memory
free(s2); // free host memory
free(u); // free host memory
magma_free_pinned(h_work); // free host memory
magma_free_pinned(r); // free host memory
magma_finalize( ); // finalize Magma
return EXIT_SUCCESS ;
}

