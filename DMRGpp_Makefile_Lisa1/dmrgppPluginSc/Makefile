include make.inc

OBJS=\
        Xcopy_.o \
	Xgemm_.o \
	Xlacpy_.o \
	makeFloat.o \
	dmrg_Xgemm_vbatch.o \
	dmrg_Xgetmatrix.o \
	dmrg_Xsetmatrix.o \
	dmrg_Xgetvector.o \
	dmrg_Xsetvector.o \
	get_total_memory.o \
	\
	apply_Htarget_vbatch.o \
	cal_kron_flops.o \
	dmrg_malloc.o \
	dmrg_vbatch.o \
	gen_patches_comb.o \
	setup_matrix.o \
	setup_vbatch.o \
	setup_sparse_batch.o \
	unsetup_vbatch.o \
	unsetup_sparse_batch.o \
	setup_nC.o \
        apply_Htarget_sparse.o \
        estimate_work.o \
        BatchedGemm.o


libdmrgppPluginSc.a: $(OBJS)
	ar cr libdmrgppPluginSc.a $(OBJS)

test_vbatch: test_vbatch.o libdmrgppPluginSc.a
	$(CC) $(CFLAGS) -o test_vbatch  test_vbatch.o \
		libdmrgppPluginSc.a $(LIBS)

BatchedGemm.o: BatchedGemm.h BatchedGemm.cpp
	$(CXX) $(CFLAGS) -c BatchedGemm.cpp

include make.dep

clean:
	touch libdmrgppPluginSc.a $(OBJS)
	rm libdmrgppPluginSc.a $(OBJS)
	touch test_vbatch test_vbatch.o
	rm test_vbatch test_vbatch.o
