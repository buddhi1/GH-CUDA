compile:

make

for a fresh make
make clean
make

usage:

./program data/Fig8-P.poly data/Fig8-Q.poly results/Fig8-R.poly
./program data/Fig14-P.poly data/Fig14-Q.poly results/Fig14-R.poly
./program data/Fig15-P.poly data/Fig15-Q.poly results/Fig15-R.poly
./program data/Fig16-P.poly data/Fig16-Q.poly results/Fig16-R.poly


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