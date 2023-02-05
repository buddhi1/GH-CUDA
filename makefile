GCC = g++

LIBB = 
LIBCUDA = -L/usr/local/cuda/lib64
LIBRA = 
CFLAGS = -m64 -O2 -Wall -c -std=c++11 -I/usr/local/include
DEBUG =

all: gppolyclip_mDIV512.o ghcuda.o 
	$(GCC) -O2 $(LFLAGS) -o program ghcuda.o gppolyclip_mDIV512.o $(LIBB) $(LIBRA) $(LIBCUDA) -lcudart

count: gppolyclip_mDIV512_withcount.o ghcuda.o 
	$(GCC) -O2 $(LFLAGS) -o program ghcuda.o gppolyclip_mDIV512_withcount.o $(LIBB) $(LIBRA) $(LIBCUDA) -lcudart

ghcuda.o: src/ghcuda.cpp
	$(GCC) $(DEBUG) $(CFLAGS) -c src/ghcuda.cpp

gppolyclip_mDIV512.o: src/gppolyclip_mDIV512.cu
	nvcc -Xptxas -O2 -o gppolyclip_mDIV512.o -c src/gppolyclip_mDIV512.cu

gppolyclip_mDIV512_withcount.o: src/gppolyclip_mDIV512_withcount.cu
	nvcc -x cu -w -m64  -o gppolyclip_mDIV512_withcount.o -c src/gppolyclip_mDIV512_withcount.cu

clean:
	rm *.o program