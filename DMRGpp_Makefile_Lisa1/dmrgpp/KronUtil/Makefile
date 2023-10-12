# DO NOT EDIT!!! Changes will be lost. Use the PsiTag system to configure instead
# This Makefile was written by configure.pl
# KronUtil


LDFLAGS += -L../../PsimagLite/lib -L../../../PsimagLite/lib -lpsimaglite

CXX = mpic++ -frecord-gcc-switches



# We're using ansi C++
CPPFLAGS += -pedantic -std=c++11

# Enable warnings and treat warnings as errors
CPPFLAGS += -Wall -Wendif-labels

LDFLAGS += -LKronUtil

CPPFLAGS += -I/opt/hdf5/code/1.12.1/install/include
LDFLAGS += -L/opt/hdf5/code/1.12.1/install/lib
LDFLAGS += -lhdf5_hl -lhdf5

# If on MacOs please say LDFLAGS += -framework Accelerate
LDFLAGS += #-llapack

LDFLAGS += -L/spack/var/spack/environments/magma/.spack-env/view/lib -lmagma_sparse -lmagma -lopenblas -L/spack/var/spack/environments/magma/.spack-env/view/lib64 -lcublas -lcudart -lcusparse -Wl,-rpath -Wl,/spack/var/spack/environments/magma/.spack-env/view/lib -Wl,-rpath -Wl,/spack/var/spack/environments/magma/.spack-env/view/lib64
CPPFLAGS += -DUSE_MAGMA -I/spack/var/spack/environments/magma/.spack-env/view/include -I/spack/var/spack/environments/magma/.spack-env/view/include


#-L/spack/opt/spack/linux-centos7-broadwell/gcc-11.4.0/openblas-0.3.23-df7sz7ivpjt6phry3olucndpb2y5h6cm/lib/ -lopenblas

CPPFLAGS += -DUSE_FLOAT -DUSE_BOOST -I/home/nitink/DMRGpp_latest/boost_1_81_0

STRIP_COMMAND = strip

# Here add -lpthread if threading is needed and also 
# set -DUSE_PTHREADS below
LDFLAGS += -lpthread

# Enable pthreads
CPPFLAGS += -DUSE_PTHREADS

CPPFLAGS += -O3

CPPFLAGS += -DNDEBUG


CPPFLAGS += -DPLUGIN_SC \
-I ../../../../dmrgppPluginSc/src \
-I ../../../dmrgppPluginSc/src \
-I ../../../../dmrgppPluginSc/include \
-I ../../../dmrgppPluginSc/include
CPPFLAGS +=    -fopenmp

#This adds linkage for the PLUGIN_SC libraries
LDFLAGS += ../../../dmrgppPluginSc/src/libdmrgppPluginSc.a
#LDFLAGS +=  -lgomp
LDFLAGS +=  -fopenmp



CPPFLAGS += -I../../../PsimagLite -I../../../PsimagLite/src -I../Engine
CPPFLAGS += 
LDFLAGS += 
all:   libkronutil.a
# test1 test2

KronUtil.o: KronUtil.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c KronUtil.cpp

util.o: util.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c util.cpp

utilComplex.o: utilComplex.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c utilComplex.cpp

csc_nnz.o: csc_nnz.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c csc_nnz.cpp


../../../PsimagLite/lib/libpsimaglite.a:
	$(MAKE) -f Makefile -C ../../../PsimagLite/lib/

Makefile.dep: KronUtil.cpp util.cpp utilComplex.cpp csc_nnz.cpp   
	$(CXX) $(CPPFLAGS) -MM KronUtil.cpp util.cpp utilComplex.cpp csc_nnz.cpp   > Makefile.dep

clean::
	rm -f core*  *.o *.a *.dep  

include Makefile.dep

libkronutil.a:  KronUtil.o  util.o  utilComplex.o  csc_nnz.o  Makefile
	ar rc libkronutil.a  KronUtil.o  util.o  utilComplex.o  csc_nnz.o 
	ranlib libkronutil.a

#test1: libkronutil.a test1.o  KronUtil.o  util.o  utilComplex.o  csc_nnz.o 
#	$(CXX) $(CFLAGS) -o test1 test1.o libkronutil.a $(LDFLAGS)

#test2: libkronutil.a test2.o
#	$(CXX) $(CFLAGS) -o test2 test2.o libkronutil.a $(LDFLAGS)

