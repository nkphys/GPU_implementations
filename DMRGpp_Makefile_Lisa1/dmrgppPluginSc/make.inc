#CC=gcc  -std=c99
#CXX=g++ -std=c++98
CXX=g++ -std=c++11
CC=g++ -std=c++11
CFLAGS=-DUSE_CUMEMADVISE -DUSE_FLOAT -Wall  -O3 -fopenmp -DUSE_MAGMA -I/spack/var/spack/environments/magma/.spack-env/view/include -I/spack/var/spack/environments/magma/.spack-env/view/include
#CFLAGS=-Wall  -g3   -fopenmp  -UUSE_MAGMA -I/spack/var/spack/environments/magma/.spack-env/view/include -I/spack/var/spack/environments/magma/.spack-env/view/include
CXXFLAGS=$(CFLAGS)
LIBS= -L/spack/var/spack/environments/magma/.spack-env/view/lib -lmagma_sparse -lmagma -lopenblas -L/spack/var/spack/environments/magma/.spack-env/view/lib64 -lcublas -lcudart -lcusparse -Wl,-rpath -Wl,/spack/var/spack/environments/magma/.spack-env/view/lib -Wl,-rpath -Wl,/spack/var/spack/environments/magma/.spack-env/view/lib64
    


#-L/spack/opt/spack/linux-centos7-broadwell/gcc-11.4.0/magma-2.7.2-k7svlaqky5xzl7ffz2os7o4ty3qs26u6/lib/  \
#    -L/spack/var/spack/environments/magma/.spack-env/view/lib -lmagma  \
#    -L/spack/var/spack/environments/magma/.spack-env/view/lib64   -lcudart -lcublas -lcusparse \
#    -L/spack/opt/spack/linux-centos7-broadwell/gcc-11.4.0/openblas-0.3.23-df7sz7ivpjt6phry3olucndpb2y5h6cm/lib/ -lopenblas -lm
