#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include "magma_v2.h"
int main ( int argc , char ** argv ){
magma_init (); // initialize Magma

magma_queue_t queue = NULL ;
magma_int_t dev =0;
magma_queue_create (dev ,& queue );
magma_int_t m = 32; // length of a,b
float *a; // a- m- vector
float *b; // b- m- vector
cudaMallocManaged (&a,m* sizeof ( float )); // unif . memory for a
cudaMallocManaged (&b,m* sizeof ( float )); // unif . memory for b

for(int j=0;j<m;j++) a[j]= sin (( float )j);
for(int j=0;j<m;j++) b[j]= cos (( float )j);

printf ("a: ");
for(int j=0;j <4;j++){
printf (" %6.4f,",a[j]); 
}
printf (" ...\n");

printf ("b: ");
for(int j=0;j <4;j++){ 
printf (" %6.4f,",b[j]);}
printf (" ...\n");

magma_sswap( m, a, 1, b, 1, queue );
cudaDeviceSynchronize ();
printf ("after magma_sswap :\n");

printf ("a: ");
for(int j=0;j <4;j++){
printf (" %6.4f,",a[j]);} 
printf (" ...\n");

printf ("b: ");
for(int j=0;j <4;j++){
printf (" %6.4f,",b[j]);}
printf (" ...\n");

magma_free (a); // free memory
magma_free (b); // free memory
magma_queue_destroy ( queue );
magma_finalize ();
return 0;
}
