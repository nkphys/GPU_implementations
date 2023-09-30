#include<stdio.h>
#include<cuda.h>
#include"magma_v2.h"
#include"magma_lapack.h"
int main (int argc ,char ** argv ) {
magma_init (); // initialize Magma
magma_queue_t queue = NULL ;

magma_int_t dev =0;
magma_queue_create(dev ,&queue);
double gpu_time , cpu_time ;
magma_int_t n=6000 , n2=n*n;
double *a, *r; // a, r - nxn matrices on the host
double *d_r ; // nxn matrix on the device
double * h_work ; // workspace
magma_int_t lwork ; // h_work size
magma_int_t * iwork ; // workspace
magma_int_t liwork ; // iwork size
double *w1 , *w2; // w1 ,w2 - vectors of eigenvalues
double error , work[1]; // used in difference computations
magma_int_t ione = 1 , info ;
double mione = -1.0;
magma_int_t incr = 1;
magma_int_t ISEED[4] = {0 ,0 ,0 ,1}; // seed
magma_dmalloc_cpu(&w1 ,n); // host memory for real
magma_dmalloc_cpu(&w2 ,n); // eigenvalues
magma_dmalloc_cpu(&a,n2 ); // host memory for a
magma_dmalloc_cpu(&r,n2 ); // host memory for r
magma_dmalloc(&d_r ,n2 ); // device memory for d_r

double aux_work[1];
magma_int_t aux_iwork[1];
magma_dsyevd_gpu(MagmaVec , MagmaLower ,n,d_r ,n,w1 ,r,n, aux_work ,-1 , aux_iwork , -1 ,& info );
lwork = (magma_int_t)aux_work[0];
liwork = aux_iwork[0];
iwork =( magma_int_t *) malloc ( liwork * sizeof ( magma_int_t ));
magma_dmalloc_cpu(&h_work , lwork ); // memory for workspace
// Randomize the matrix a and copy a -> r
lapackf77_dlarnv(&ione ,ISEED ,&n2 ,a);
lapackf77_dlacpy(MagmaFullStr ,&n ,&n,a ,&n,r ,&n);
magma_dsetmatrix( n, n, a, n, d_r ,n, queue ); // copy a -> d_r
// compute the eigenvalues and eigenvectors for a symmetric ,
// real nxn matrix ; Magma version
 gpu_time = magma_sync_wtime ( NULL );

magma_dsyevd_gpu(MagmaVec,MagmaLower,n,d_r,n,w1,r,n,h_work,lwork,iwork,liwork,&info);
gpu_time = magma_sync_wtime ( NULL ) - gpu_time ;
printf (" dsyevd_gpu gpu time : %7.5f sec .\n",gpu_time ); // Magma time
// Lapack version
cpu_time = magma_wtime ();
lapackf77_dsyevd ("V","L" ,&n,a ,&n,w2 ,h_work ,& lwork ,iwork ,&liwork ,& info );
cpu_time = magma_wtime () - cpu_time ;
printf (" dsyevd cpu time : %7.5f sec .\n",cpu_time ); // Lapack
// difference in eigenvalues

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
