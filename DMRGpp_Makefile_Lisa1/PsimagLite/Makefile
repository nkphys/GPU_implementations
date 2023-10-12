# DO NOT EDIT!!! Changes will be lost. Use the PsiTag system to configure instead
# This Makefile was written by configure.pl
# PsimagLite/lib


CXX = mpic++ -frecord-gcc-switches



# We're using ansi C++
CPPFLAGS += -pedantic -std=c++11

# Enable warnings and treat warnings as errors
CPPFLAGS += -Wall -Wendif-labels

CPPFLAGS += -I/opt/hdf5/code/1.12.1/install/include
LDFLAGS += -L/opt/hdf5/code/1.12.1/install/lib
LDFLAGS += -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5


# If on MacOs please say LDFLAGS += -framework Accelerate
#LDFLAGS += -llapack
#LDFLAGS += -lblas

LDFLAGS += -L/spack/var/spack/environments/magma/.spack-env/view/lib -lmagma_sparse -lmagma -lopenblas -L/spack/var/spack/environments/magma/.spack-env/view/lib64 -lcublas -lcudart -lcusparse -Wl,-rpath -Wl,/spack/var/spack/environments/magma/.spack-env/view/lib -Wl,-rpath -Wl,/spack/var/spack/environments/magma/.spack-env/view/lib64

CPPFLAGS += -DUSE_MAGMA -I/spack/var/spack/environments/magma/.spack-env/view/include -I/spack/var/spack/environments/magma/.spack-env/view/include
CPPFLAGS += -DUSE_FLOAT -DUSE_BOOST -I/home/nitink/DMRGpp_latest/boost_1_81_0



CPPFLAGS += -DNDEBUG

CPPFLAGS += -O3

# Here add -lpthread if threading is needed and also 
# set -DUSE_PTHREADS below
LDFLAGS += -lpthread

# Enable pthreads
CPPFLAGS += -DUSE_PTHREADS



all: libpsimaglite.a
libpsimaglite.a: Makefile MersenneTwister.o Matrix.o Mpi.o Concurrency.o ProgressIndicator.o MemResolv.o PsimagLite.o PsiBase64.o SpecialFunctions.o Io/TypeToH5.o TridiagonalMatrix.o PredicateSimple.o  AinurSpirit.o AinurConvert.o
	ar rcs libpsimaglite.a MersenneTwister.o Matrix.o Mpi.o Concurrency.o ProgressIndicator.o MemResolv.o PsimagLite.o PsiBase64.o SpecialFunctions.o TypeToH5.o TridiagonalMatrix.o PredicateSimple.o  AinurSpirit.o AinurConvert.o
MersenneTwister.o: ../src/MersenneTwister.cpp ../src/MersenneTwister.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/MersenneTwister.cpp
Matrix.o: ../src/Matrix.cpp ../src/Matrix.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/Matrix.cpp
Mpi.o: ../src/Mpi.cpp ../src/Mpi.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/Mpi.cpp
Concurrency.o: ../src/Concurrency.cpp ../src/Concurrency.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/Concurrency.cpp
ProgressIndicator.o: ../src/ProgressIndicator.cpp ../src/ProgressIndicator.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/ProgressIndicator.cpp
MemResolv.o: ../src/MemResolv.cpp ../src/MemResolv.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/MemResolv.cpp
PsimagLite.o: ../src/PsimagLite.cpp ../src/PsimagLite.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/PsimagLite.cpp
PsiBase64.o: ../src/PsiBase64.cpp ../src/PsiBase64.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/PsiBase64.cpp
SpecialFunctions.o: ../src/SpecialFunctions.cpp ../src/SpecialFunctions.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/SpecialFunctions.cpp
Io/TypeToH5.o: ../src/Io/TypeToH5.cpp ../src/Io/TypeToH5.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/Io/TypeToH5.cpp
TridiagonalMatrix.o: ../src/TridiagonalMatrix.cpp ../src/TridiagonalMatrix.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/TridiagonalMatrix.cpp
PredicateSimple.o: ../src/PredicateSimple.cpp ../src/PredicateSimple.h Makefile
	$(CXX) $(CPPFLAGS) -c ../src/PredicateSimple.cpp
AinurSpirit.o: ../src/Ainur/AinurSpirit.cpp ../src/Ainur/AinurSpirit.h Makefile ../src/Ainur/AinurSpirit1.cpp
	$(CXX) $(CPPFLAGS) -c ../src/Ainur/AinurSpirit.cpp

AinurConvert.o: ../src/Ainur/AinurConvert.cpp ../src/Ainur/AinurConvert1.cpp ../src/Ainur/AinurConvert.hh Makefile
	$(CXX) $(CPPFLAGS) -c ../src/Ainur/AinurConvert.cpp

Makefile.dep: ../src/MersenneTwister.cpp ../src/Matrix.cpp ../src/Mpi.cpp ../src/Concurrency.cpp ../src/ProgressIndicator.cpp ../src/MemResolv.cpp ../src/PsimagLite.cpp ../src/PsiBase64.cpp ../src/SpecialFunctions.cpp ../src/Io/TypeToH5.cpp ../src/TridiagonalMatrix.cpp ../src/PredicateSimple.cpp  ../src/Ainur/AinurSpirit.cpp ../src/Ainur/AinurConvert.cpp ../src/Ainur/AinurConvert1.cpp
	$(CXX) $(CPPFLAGS) -MM  ../src/MersenneTwister.cpp ../src/Matrix.cpp ../src/Mpi.cpp ../src/Concurrency.cpp ../src/ProgressIndicator.cpp ../src/MemResolv.cpp ../src/PsimagLite.cpp ../src/PsiBase64.cpp ../src/SpecialFunctions.cpp ../src/Io/TypeToH5.cpp ../src/TridiagonalMatrix.cpp ../src/PredicateSimple.cpp   > Makefile.dep

clean: Makefile.dep
	rm -f core* *.o *.dep *.a

include Makefile.dep

