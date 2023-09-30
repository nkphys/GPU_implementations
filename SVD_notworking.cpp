#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "magma_v2.h"
#include "magma_lapack.h"
#define min(a,b)   (((a)<(b))?(a):(b))
int main ( int argc , char ** argv )
{
magma_init (); // initialize Magma
real_Double_t gpu_time , cpu_time ;

magma_int_t m=320 , n =320 , n2=m*n, min_mn = min(m,n);
float *a, *r; // a,r - mxn matrices
float *u, *vt; //u - mxm matrix , vt - nxn matrix
float *s1 , *s2; // vectors of singular values
magma_int_t info ;
magma_int_t ione = 1;
float work[1] , error = 1.; // used in difference computations
float mone = -1.0 , * h_work ; // h_work - workspace
magma_int_t lwork ; // workspace size
magma_int_t ISEED[4] = {0 ,0 ,0 ,1}; // seed


//printf ("here \n");

cudaMallocManaged(&a , n2*sizeof(float)); // unif . memory for a
cudaMallocManaged(&vt , n*n*sizeof(float)); // uni. memory for vt
cudaMallocManaged(&u , m*m*sizeof(float)); // unif . memory for u
cudaMallocManaged(&s1 , min_mn*sizeof(float)); // uni. mem. for s1
cudaMallocManaged(&s2 , min_mn*sizeof(float)); // uni. mem. for s2
cudaMallocManaged(&r , n2*sizeof(float)); // unif . memory for r
magma_int_t nb = magma_get_sgesvd_nb(m,n); // optim . block size
lwork = min_mn*min_mn +2*min_mn +2*min_mn*nb;
cudaMallocManaged(&h_work , lwork*sizeof(float)); //m.f. h_work

//printf ("here \n");

// Randomize the matrix a
lapackf77_slarnv(&ione ,ISEED ,&n2 ,a);

printf ("a: ");
for(int j=0;j <20;j++){
printf (" %6.4f,",a[j]);
}
printf (".........,");
for(int j=0;j<20;j++){
printf (" %6.4f,",a[n2-20+j]);
}
printf ("\n");


lapackf77_slacpy(MagmaFullStr ,&m ,&n,a ,&m,r ,&m); //a- >r


printf ("r: ");
for(int j=0;j <20;j++){
printf (" %6.4f,",r[j]);
}
printf (".........,");
for(int j=0;j<20;j++){
printf (" %6.4f,",r[n2-20+j]);
}

printf ("\n");



// MAGMA
 gpu_time = magma_wtime ();
// compute the singular value decomposition of r ( copy of a)
// and optionally the left and right singular vectors :
// r = u* sigma *vt; the diagonal elements of sigma (s1 array )
// are the singular values of a in descending order
// the first min (m,n) columns of u contain the left sing . vec 
// the first min (m,n) columns of vt contain the right sing .vec .

printf ("here1 \n");

magma_sgesvd(MagmaNoVec,MagmaNoVec,m,n,r,m,s1,u,m,vt,n,h_work,lwork,&info);

printf ("here2 \n");

gpu_time = magma_wtime () - gpu_time ;
printf (" sgesvd gpu time : %7.5f \n", gpu_time ); // Magma time


// LAPACK
cpu_time = magma_wtime ();
lapackf77_sgesvd ("N","N" ,&m ,&n,a ,&m,s2 ,u ,&m, vt ,&n, h_work , &lwork ,&info);


cpu_time = magma_wtime () - cpu_time ;
printf (" sgesvd cpu time : %7.5f \n", cpu_time ); // Lapack time
// difference
error = lapackf77_slange ("f" ,&min_mn ,& ione ,s1 ,& min_mn , work );
blasf77_saxpy (& min_mn ,& mone ,s1 ,& ione ,s2 ,& ione );
error = lapackf77_slange ("f" ,&min_mn ,& ione ,s2 ,& min_mn , work );
// error ;
printf (" difference : %e\n", error ); // difference in singul .
// values

magma_free (a); // free memory
magma_free (vt ); // free memory
magma_free (s1 ); // free memory
magma_free (s2 ); // free memory
magma_free (u); // free memory
magma_free ( h_work ); // free memory
magma_free (r); // free memory
magma_finalize ( ); // finalize Magma
return EXIT_SUCCESS ;
}

