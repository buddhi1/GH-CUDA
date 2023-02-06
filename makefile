GCC = g++

LIBB = 
LIBCUDA = -L/usr/local/cuda/lib64
LIBRA = 
CFLAGS = -m64 -O2 -Wall -c -std=c++11 -I/usr/local/include
DEBUG =

filters: gppolyclip_mDIV512.o ghcuda.o 
	$(GCC) -O2 $(LFLAGS) -o program ghcuda.o gppolyclip_mDIV512.o $(LIBB) $(LIBRA) $(LIBCUDA) -lcudart

count: gppolyclip_mDIV512_withcount.o ghcuda.o 
	$(GCC) -O2 $(LFLAGS) -o program ghcuda.o gppolyclip_mDIV512_withcount.o $(LIBB) $(LIBRA) $(LIBCUDA) -lcudart

cmfcf: gppolyclip_mDIV512-CMF_CF.o ghcuda.o 
	$(GCC) -O2 $(LFLAGS) -o program ghcuda.o gppolyclip_mDIV512-CMF_CF.o $(LIBB) $(LIBRA) $(LIBCUDA) -lcudart

lsmf: gppolyclip_mDIV512-LSMF.o ghcuda.o 
	$(GCC) -O2 $(LFLAGS) -o program ghcuda.o gppolyclip_mDIV512-LSMF.o $(LIBB) $(LIBRA) $(LIBCUDA) -lcudart

nofilters: gppolyclip_mDIV512-noFilters.o ghcuda.o 
	$(GCC) -O2 $(LFLAGS) -o program ghcuda.o gppolyclip_mDIV512-noFilters.o $(LIBB) $(LIBRA) $(LIBCUDA) -lcudart



ghcuda.o: src/ghcuda.cpp
	$(GCC) $(DEBUG) $(CFLAGS) -c src/ghcuda.cpp

gppolyclip_mDIV512.o: src/gppolyclip_mDIV512.cu
	nvcc -Xptxas -O2 -o gppolyclip_mDIV512.o -c src/gppolyclip_mDIV512.cu

gppolyclip_mDIV512_withcount.o: src/gppolyclip_mDIV512_withcount.cu
	nvcc -x cu -w -m64  -o gppolyclip_mDIV512_withcount.o -c src/gppolyclip_mDIV512_withcount.cu

gppolyclip_mDIV512-CMF_CF.o: src/gppolyclip_mDIV512-CMF_CF.cu
	nvcc -Xptxas -O2 -o gppolyclip_mDIV512-CMF_CF.o -c src/gppolyclip_mDIV512-CMF_CF.cu

gppolyclip_mDIV512-LSMF.o: src/gppolyclip_mDIV512-LSMF.cu
	nvcc -Xptxas -O2 -o gppolyclip_mDIV512-LSMF.o -c src/gppolyclip_mDIV512-LSMF.cu

gppolyclip_mDIV512-noFilters.o: src/gppolyclip_mDIV512-noFilters.cu
	nvcc -Xptxas -O2 -o gppolyclip_mDIV512-noFilters.o -c src/gppolyclip_mDIV512-noFilters.cu

foster:
	$(GCC) -O2 -std=c++11 optimizedFostersAlgorithm/polyclip_time.cpp -o polyclip

clean:
	rm *.o program

cleanresults:
	rm results/*.poly

cleanfoster:
	rm polyclip

cleanfosterresults:
	rm optimizedFostersAlgorithm/results/*.poly