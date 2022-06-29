compile:

make

for a fresh make
make clean
make


shapefile to WKT
ogr2ogr -f CSV continents.csv continents/continent.shp -lco GEOMETRY=AS_WKT

Map Polygon in Google maps
https://arthur-e.github.io/Wicket/sandbox-gmaps3.html

Map a shape file
https://mygeodata.cloud/converter/shp-to-csv 

usage:

./program data/Fig1-P.poly data/Fig1-Q.poly results/Fig1-R.poly
./program data/Fig8-P.poly data/Fig8-Q.poly results/Fig8-R.poly
./program data/Fig14-P.poly data/Fig14-Q.poly results/Fig14-R.poly
./program data/Fig14-P-clockWise.poly data/Fig14-Q.poly results/Fig14-R.poly
./program data/Fig15-P.poly data/Fig15-Q.poly results/Fig15-R.poly

./program data/Fig16-P.poly data/Fig16-Q.poly results/Fig16-R.poly
./program data/Fig16-Q.poly data/Fig16-P.poly results/Fig16-R.poly

./program data/Fig17-P.poly data/Fig17-Q.poly results/Fig17-R.poly
./program data/Fig17-P-Clockwise.poly data/Fig17-Q.poly results/Fig17-R.poly

./program data/Fig18-P.poly data/Fig18-Q.poly results/Fig18-R.poly
./program data/Fig19-P.poly data/Fig19-Q.poly results/Fig19-R.poly

./program data/readPolygon/s.txt data/readPolygon/c.txt results/Fig-large-R.poly


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
1. The indicese are different in the neighbor arrays due to the complexity. 
    The issue is, we assumed the order of intersection found using PxQ and QxP are same. 
    But, it is different and whole idea of neighbor is broken now
    FIX:
        Need to establish PxQ and QxP intersections are found in the same order - correct assumption
        Invent another index to keep track - implemented this
        Try to do the same intersection point is both PxQ and QxP and see - did not try. this will not work due to complexity
    PROGRESS: this issue is fixed in new_mapping git branch

2. Found some errors in Init label. 
    REASONS: Still neighbor connection could be broken. Not sure yet
    FIX:   
        count2 updated was not made in the (id>=sizeP) section
        sorted alpha issue fixed
        sorted neighbors need a barrier. Currently, 2 kernel calls used as the fix. Not optimized
    PROGRESS: fixed

TO DO CODE
1. Shared memory usage
2. CMBR filter using GPU
3. Point-in-polygon test GPU
4. manage when size of PP or QQ >0 WORKING ON THIS
5. Since we have neighbor map and exact locations of neighbor location in the array, we should
    be able to run intersection calculation using max(m,n) processors in mim(m,n) time by 
    removing seperate section for Q info handling. 
6. Parallelize later parts using OpenMP or pthread
7. Optimization: If prefix sum of the given id is 0, no need to run the kernel for that pid



VERSIONS
1. gppolyclip_v1.cu: Single component menthods tetsed good. 
    Multi needs fixing in sorting and some other places



Special 6/29
PPVertexPointers[pi]=V;
// if(i==12381*2) cout<<PPVertexPointers[pi]->p.x<<", "<<PPVertexPointers[pi]->p.y<<" >> "<<V->p.x<<", "<<V->p.y<<endl;
pi++;

This was at the begining of the loop. But got error. Shifted to bottom

P->ne_10m_ocean
Q->ne_10m_land

6822 polygons found
6838 polygons found
P Polygons
100612 0
66475 36
15547 2742
10887 2741
9461 5978
3838 2854
3706 2737
3612 6017
3163 2723
3062 2740
2478 3346
2387 6805
2358 3293
2253 2726
2207 2385
1847 6785
1791 6020
1778 717
1738 6778
1568 2846
Q Polygons
81511 4
66475 1
15956 0
15547 33
10887 30
9461 3
3838 42
3706 25
3612 8
3163 19
3062 26
2478 17
2387 5
2358 12
2253 23
2207 28
1847 97
1791 10
1778 92
1738 96