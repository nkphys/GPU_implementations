//#include <cmath>
//#include <cstdio>
//#include <cstdlib>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "magma_v2.h"

int main(int argc, char **argv) {
  magma_init(); // initialize Magma
  magma_queue_t queue = NULL;
  magma_int_t dev = 0;
  magma_queue_create(dev, &queue);
  magma_int_t m = 32; // length of a
  float *a;           // a - m- vector on the host
  float *b;           // b - m- vector on the host
  float *a2;
  float *b2;
  float *d_a; // d_a - m- vector a on the device
  float *d_b; // d_b - m- vector a on the device

  magma_int_t err;
  // allocate the vectors on the host
  err = magma_smalloc_cpu(&a, m); // host mem. for a
  err = magma_smalloc_cpu(&b, m); // host mem. for b

  err = magma_smalloc_cpu(&a2, m); // host mem. for a
  err = magma_smalloc_cpu(&b2, m); // host mem. for b

  // allocate the vector on the device
  err = magma_smalloc(&d_a, m); // device memory for a
  err = magma_smalloc(&d_b, m); // device memory for b
  // a={ sin (0) , sin (1) ,... , sin (m -1)}
  for (int j = 0; j < m; j++)
    a[j] = sin((float)j);
  // b={ cos (0) , cos (1) ,... , cos (m -1)}
  for (int j = 0; j < m; j++)
    b[j] = cos((float)j);

  printf("a: ");
  for (int j = 0; j < 4; j++)
    printf(" %6.4f,", a[j]);
  printf(" ...\n");

  printf("b: ");
  for (int j = 0; j < 4; j++)
    printf(" %6.4f,", b[j]);
  printf(" ...\n");

  // copy data from host to device
  magma_setvector(m, sizeof(float), a, 1, d_a, 1, queue); // copy a -> d_a
  magma_setvector(m, sizeof(float), b, 1, d_b, 1, queue); // copy b -> d_b

  printf("d_a from gpu: ");
  magma_sprint_gpu(m, 1, d_a, m, queue);
  printf("d_b from gpu: ");
  magma_sprint_gpu(m, 1, d_b, m, queue);

  // swap the vectors
  magma_sswap(m, d_a, 1, d_b, 1, queue);
  magma_sgetvector(m, d_a, 1, a2, 1, queue); // copy d_a -> a2
  magma_sgetvector(m, d_b, 1, b2, 1, queue); // copy d_b -> b2

  printf(" after magma_sswap :\n");

  printf("a2: ");
  for (int j = 0; j < 4; j++)
    printf(" %6.4f,", a2[j]);
  printf(" ...\n");

  printf("b2: ");
  for (int j = 0; j < 4; j++)
    printf(" %6.4f,", b2[j]);
  printf(" ...\n");

  free(a); // free host memory
  free(b); // free host memory
  free(a2);
  free(b2);
  magma_free(d_a); // free device memory
  magma_free(d_b); // free device memory
  magma_queue_destroy(queue);
  magma_finalize();
  return 0;
}
