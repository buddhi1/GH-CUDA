compile:

make

for a fresh make
make clean
make

usage:

./program data/Fig1-P.poly data/Fig1-Q.poly results/Fig1-R.poly
./program data/Fig8-P.poly data/Fig8-Q.poly results/Fig8-R.poly
./program data/Fig14-P.poly data/Fig14-Q.poly results/Fig14-R.poly
./program data/Fig14-P-clockWise.poly data/Fig14-Q.poly results/Fig14-R.poly
./program data/Fig15-P.poly data/Fig15-Q.poly results/Fig15-R.poly
./program data/Fig16-P.poly data/Fig16-Q.poly results/Fig16-R.poly

./program data/Fig17-P.poly data/Fig17-Q.poly results/Fig17-R.poly
./program data/Fig17-P-Clockwise.poly data/Fig17-Q.poly results/Fig17-R.poly

./program data/Fig18-P.poly data/Fig18-Q.poly results/Fig18-R.poly
./program data/Fig19-P.poly data/Fig19-Q.poly results/Fig19-R.poly


polyclip [-union] input1.poly input2.poly output.poly

./polyclip -union ../examples/Fig20-E1.poly ../examples/Fig20-E2.poly Fig20-E12.poly
./polyclip -union Fig20-E12.poly ../examples/Fig20-E3.poly Fig20-E123.poly
./polyclip -union Fig20-E123.poly ../examples/Fig20-E4.poly Fig20-E1234.poly
./polyclip -union Fig20-E1234.poly ../examples/Fig20-E5.poly Fig20-E.poly


./polyclip -union ../examples/Fig20-M1.poly ../examples/Fig20-M2.poly Fig20-M12.poly
./polyclip -union Fig20-M12.poly ../examples/Fig20-M3.poly Fig20-M.poly


./polyclip -union ../examples/Fig20-H1.poly ../examples/Fig20-H2.poly Fig20-H.poly


./polyclip Fig20-E.poly Fig20-M.poly Fig20-EM.poly
./polyclip Fig20-EM.poly Fig20-H.poly Fig20-R.poly


ERROR
The indicese are different in the neighbor arrays due to the complexity. 
The issue is, we assumed the order of intersection found using PxQ and QxP are same. 
But, it is diofferent and whole idea of neighbor is broken now

FIX
Need to establish PxQ and QxP intersections are found in the same order 
OR invent anothe index to keep track
Try to do the same intersection point is both PxQ and QxP and see

TO DO CODE
1. Shared memory usage
2. CMBR filter using GPU
3. Point-in-polygon test GPU