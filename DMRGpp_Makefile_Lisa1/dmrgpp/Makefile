# DO NOT EDIT!!! Changes will be lost. Use the PsiTag system to configure instead
# This Makefile was written by configure.pl
# DMRG++


LDFLAGS += -L../../PsimagLite/lib -L../../../PsimagLite/lib -lpsimaglite

CXX = mpic++ -frecord-gcc-switches



# We're using ansi C++
CPPFLAGS += -pedantic -std=c++11

# Enable warnings and treat warnings as errors
CPPFLAGS += -Wall -Wendif-labels

LDFLAGS += -LKronUtil

CPPFLAGS += -I/opt/hdf5/code/1.12.1/install/include
LDFLAGS += -L/opt/hdf5/code/1.12.1/install/lib
LDFLAGS += -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5

# If on MacOs please say LDFLAGS += -framework Accelerate
LDFLAGS += #-llapack

LDFLAGS += #-lblas

CPPFLAGS += -DUSE_FLOAT -DUSE_BOOST -I/home/nitink/DMRGpp_latest/boost_1_81_0

LDFLAGS += -L/spack/var/spack/environments/magma/.spack-env/view/lib -lmagma_sparse -lmagma -lopenblas -L/spack/var/spack/environments/magma/.spack-env/view/lib64 -lcublas -lcudart -lcusparse -Wl,-rpath -Wl,/spack/var/spack/environments/magma/.spack-env/view/lib -Wl,-rpath -Wl,/spack/var/spack/environments/magma/.spack-env/view/lib64
CPPFLAGS += -DUSE_MAGMA -I/spack/var/spack/environments/magma/.spack-env/view/include -I/spack/var/spack/environments/magma/.spack-env/view/include


#CPPFLAGS += -DUSE_MAGMA -I/spack/var/spack/environments/magma/.spack-env/view/include -I/spack/var/spack/environments/magma/.spack-env/view/include
#LDFLAGS += -L/spack/var/spack/environments/magma/.spack-env/view/lib -lmagma -L/spack/var/spack/environments/magma/.spack-env/view/lib64 -lcudart -lcublas -lcusparse -L/spack/opt/spack/linux-centos7-broadwell/gcc-11.4.0/openblas-0.3.23-df7sz7ivpjt6phry3olucndpb2y5h6cm/lib/ -lopenblas -lm


STRIP_COMMAND = strip

# Here add -lpthread if threading is needed and also 
# set -DUSE_PTHREADS below
LDFLAGS += -lpthread

# Enable pthreads
CPPFLAGS += -DUSE_PTHREADS

CPPFLAGS += -O3

CPPFLAGS += -DNDEBUG


CPPFLAGS += -DPLUGIN_SC \
-I ../../../dmrgppPluginSc/src \
-I ../../dmrgppPluginSc/src \
-I ../../../dmrgppPluginSc/include \
-I ../../dmrgppPluginSc/include
CPPFLAGS +=    -fopenmp

#This adds linkage for the PLUGIN_SC libraries
LDFLAGS += ../../dmrgppPluginSc/src/libdmrgppPluginSc.a
#LDFLAGS +=  -lgomp
LDFLAGS +=  -fopenmp



CPPFLAGS += -I ../../PsimagLite -I ../../PsimagLite/src -I Engine
CPPFLAGS += 
LDFLAGS += 
all: GitRevision.h observe toolboxdmrg manyOmegas procOmegas dmrg  operator

Provenance.o: Provenance.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c Provenance.cpp

Su2Related.o: Su2Related.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c Su2Related.cpp

ProgramGlobals.o: ProgramGlobals.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c ProgramGlobals.cpp

Utils.o: Utils.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c Utils.cpp

Qn.o: Qn.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c Qn.cpp

observe.o: observe.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c observe.cpp

observe: observe.o ProgramGlobals.o Provenance.o Utils.o Su2Related.o Qn.o  ObserveDriver0.o ObserveDriver1.o ObserveDriver2.o  
	$(CXX) -o  observe observe.o ProgramGlobals.o Provenance.o Utils.o Su2Related.o Qn.o  ObserveDriver0.o ObserveDriver1.o ObserveDriver2.o  $(LDFLAGS)  $(CPPFLAGS)
	$(STRIP_COMMAND) observe

toolboxdmrg.o: toolboxdmrg.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c toolboxdmrg.cpp

toolboxdmrg: toolboxdmrg.o ProgramGlobals.o Provenance.o Utils.o Qn.o 
	$(CXX) -o  toolboxdmrg toolboxdmrg.o ProgramGlobals.o Provenance.o Utils.o Qn.o $(LDFLAGS)  $(CPPFLAGS)
	$(STRIP_COMMAND) toolboxdmrg

ObserveDriver0.o: ObserveDriver0.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c ObserveDriver0.cpp

ObserveDriver1.o: ObserveDriver1.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c ObserveDriver1.cpp

ObserveDriver2.o: ObserveDriver2.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c ObserveDriver2.cpp

manyOmegas.o: manyOmegas.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c manyOmegas.cpp

manyOmegas: manyOmegas.o libdmrgpp.a KronUtil/libkronutil.a 
	$(CXX) -o  manyOmegas manyOmegas.o libdmrgpp.a KronUtil/libkronutil.a $(LDFLAGS)  $(CPPFLAGS)
	$(STRIP_COMMAND) manyOmegas

procOmegas.o: procOmegas.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c procOmegas.cpp

procOmegas: procOmegas.o libdmrgpp.a KronUtil/libkronutil.a 
	$(CXX) -o  procOmegas procOmegas.o libdmrgpp.a KronUtil/libkronutil.a $(LDFLAGS)  $(CPPFLAGS)
	$(STRIP_COMMAND) procOmegas

DmrgDriver0.o: DmrgDriver0.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver0.cpp

DmrgDriver1.o: DmrgDriver1.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver1.cpp

DmrgDriver2.o: DmrgDriver2.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver2.cpp

DmrgDriver3.o: DmrgDriver3.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver3.cpp

DmrgDriver4.o: DmrgDriver4.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver4.cpp

DmrgDriver5.o: DmrgDriver5.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver5.cpp

DmrgDriver6.o: DmrgDriver6.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver6.cpp

DmrgDriver7.o: DmrgDriver7.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver7.cpp

DmrgDriver8.o: DmrgDriver8.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver8.cpp

DmrgDriver9.o: DmrgDriver9.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver9.cpp

DmrgDriver10.o: DmrgDriver10.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver10.cpp

DmrgDriver11.o: DmrgDriver11.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c DmrgDriver11.cpp

dmrg.o: dmrg.cpp  Makefile  
	$(CXX) $(CPPFLAGS) -c dmrg.cpp

dmrg: dmrg.o Provenance.o Utils.o Qn.o  ProgramGlobals.o Su2Related.o DmrgDriver0.o  DmrgDriver1.o  DmrgDriver2.o  DmrgDriver3.o  DmrgDriver4.o  DmrgDriver5.o  DmrgDriver6.o  DmrgDriver7.o  DmrgDriver8.o  DmrgDriver9.o  DmrgDriver10.o  DmrgDriver11.o  KronUtil/libkronutil.a 
	$(CXX) -o  dmrg dmrg.o Provenance.o Utils.o Qn.o  ProgramGlobals.o Su2Related.o DmrgDriver0.o  DmrgDriver1.o  DmrgDriver2.o  DmrgDriver3.o  DmrgDriver4.o  DmrgDriver5.o  DmrgDriver6.o  DmrgDriver7.o  DmrgDriver8.o  DmrgDriver9.o  DmrgDriver10.o  DmrgDriver11.o  KronUtil/libkronutil.a $(LDFLAGS)  $(CPPFLAGS)
	$(STRIP_COMMAND) dmrg


 ../../PsimagLite/lib/libpsimaglite.a:
	$(MAKE) -f Makefile -C  ../../PsimagLite/lib/

Makefile.dep: Provenance.cpp Su2Related.cpp ProgramGlobals.cpp Utils.cpp Qn.cpp observe.cpp toolboxdmrg.cpp ObserveDriver0.cpp ObserveDriver1.cpp ObserveDriver2.cpp manyOmegas.cpp procOmegas.cpp DmrgDriver0.cpp DmrgDriver1.cpp DmrgDriver2.cpp DmrgDriver3.cpp DmrgDriver4.cpp DmrgDriver5.cpp DmrgDriver6.cpp DmrgDriver7.cpp DmrgDriver8.cpp DmrgDriver9.cpp DmrgDriver10.cpp DmrgDriver11.cpp dmrg.cpp   
	$(CXX) $(CPPFLAGS) -MM Provenance.cpp Su2Related.cpp ProgramGlobals.cpp Utils.cpp Qn.cpp observe.cpp toolboxdmrg.cpp ObserveDriver0.cpp ObserveDriver1.cpp ObserveDriver2.cpp manyOmegas.cpp procOmegas.cpp DmrgDriver0.cpp DmrgDriver1.cpp DmrgDriver2.cpp DmrgDriver3.cpp DmrgDriver4.cpp DmrgDriver5.cpp DmrgDriver6.cpp DmrgDriver7.cpp DmrgDriver8.cpp DmrgDriver9.cpp DmrgDriver10.cpp DmrgDriver11.cpp dmrg.cpp   > Makefile.dep

clean::
	rm -f core* observe toolboxdmrg manyOmegas procOmegas dmrg  *.o *.a *.dep  

include Makefile.dep

.PHONY: GitRevision.h

GitRevision.h:
	./createGitRevision.pl GitRevision.h

operator: dmrg
	cp dmrg operator

testQn: testQn.o Qn.o
	$(CXX) Qn.o testQn.o $(LDFLAGS) -o testQn

KronUtil/libkronutil.a:
	$(MAKE) -C KronUtil

libdmrgpp.a: Makefile Provenance.o ProgramGlobals.o Qn.o Utils.o
	ar rcs libdmrgpp.a Provenance.o ProgramGlobals.o Qn.o Utils.o

../doc/manual.pdf: ../doc/manual.tex
	cd ../doc; pdflatex manual.tex; pdflatex manual.tex; pdflatex manual.tex

../doc/manual.tex: ../doc/manual.ptex
	cd ../doc; find ../src -iname "*.h" -or -iname "*.cpp" | ../../PsimagLite/scripts/doc.pl manual.ptex

clean::
	$(MAKE) -C KronUtil clean
