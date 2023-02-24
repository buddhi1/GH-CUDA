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


----------------------------------------------------------------------------------------
Test cases command line
./program data/readPolygon/s.txt data/readPolygon/c.txt results/Fig-large-R.poly

./program data/syntheticO_n/lakes_851348.txt  data/syntheticO_n/lakes-synthetic_1_1.txt  results/Fig-large-R.poly
./program data/readPolygon/s.txt  data/syntheticO_n/s-synthetic_15_10.txt  results/Fig-large-R.poly
./program data/readPolygon/c.txt  data/syntheticO_n/c-synthetic_100_0.txt  results/Fig-large-R.poly

./program data/syntheticO_n/lakes_851348.txt  data/syntheticO_n/lakes_851348.txt  results/Fig-large-R.poly
./program data/readPolygon/s.txt  data/readPolygon/s.txt  results/Fig-large-R.poly
./program data/readPolygon/c.txt  data/readPolygon/c.txt results/Fig-large-R.poly


./program data/syntheticAlltoALL/worst-syntheticP_100.txt data/syntheticAlltoALL/worst-syntheticQ_100.txt results/Fig-syn.poly
./program data/syntheticAlltoALL/worst-syntheticP_500.txt data/syntheticAlltoALL/worst-syntheticQ_500.txt results/Fig-syn.poly
./program data/syntheticAlltoALL/worst-syntheticP_700.txt data/syntheticAlltoALL/worst-syntheticQ_700.txt results/Fig-syn.poly
./program data/syntheticAlltoALL/worst-syntheticP_1000.txt data/syntheticAlltoALL/worst-syntheticQ_1000.txt results/Fig-syn.poly
./program data/syntheticAlltoALL/worst-syntheticP_1500.txt data/syntheticAlltoALL/worst-syntheticQ_1500.txt results/Fig-syn.poly
./program data/syntheticAlltoALL/worst-syntheticP_2000.txt data/syntheticAlltoALL/worst-syntheticQ_2000.txt results/Fig-syn.poly
----------------------------------------------------------------------------------------


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
1. Shared memory usage - Done
2. CMBR filter using GPU - Have 2 version only bool array and reduced block size
3. Point-in-polygon test GPU - Not touched
4. manage when size of PP or QQ >0 not done
5. Since we have neighbor map and exact locations of neighbor location in the array, we should
    be able to run intersection calculation using max(m,n) processors in mim(m,n) time by 
    removing seperate section for Q info handling. - Done
6. Parallelize later parts using OpenMP or pthread - We have time to work on this just to reduce overall time. But cannot talk about it in the paper
7. Optimization: If prefix sum of the given id is 0, no need to run the kernel for that pid - Done



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



Results
-------------------------------------------------------
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 0
Shape file2: ../datasets/ne_10m_land.csv QQID: 1
PP Polygon size 100612 QQ Polygon size 66475
PP Count 100612 QQ Count 66475

P overlap count with CMBR 7107 Q overlap count with CMBR 66475 

Non-degen count P 0 *****--- Q 0
Intersection count P 0 *****--- Q 0

gpuCountIntersections kernel exe time(ms) 419.812653
prefixsum kernels exe time(ms) 2.825024
gpuNeighborMap kernel exe time(ms) 0.011616
gpuCalculateIntersections kernel exe time(ms) 0.029952
gpuSortPolyQ kernel exe time(ms) 0.008160
gpuCalculateInitLabel kernel exe time(ms) 0.009984

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 1 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 28 vertices removed

R has 1 component with 66447 vertices

All time in microseconds
Time: Total : 567794
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 0
Shape file2: ../datasets/ne_10m_land.csv QQID: 3
PP Polygon size 100612 QQ Polygon size 9461
PP Count 100612 QQ Count 9461

P overlap count with CMBR 2686 Q overlap count with CMBR 9461 

Non-degen count P 4 *****--- Q 4
Intersection count P 4 *****--- Q 4

gpuCountIntersections kernel exe time(ms) 90.043617
prefixsum kernels exe time(ms) 1.834816
gpuNeighborMap kernel exe time(ms) 6.232064
gpuCalculateIntersections kernel exe time(ms) 0.717536
gpuSortPolyQ kernel exe time(ms) 0.018336
gpuCalculateInitLabel kernel exe time(ms) 0.010240

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 0 vertices removed

R has 2 components with 4 + 9462 = 9466 vertices

All time in microseconds
Time: Total : 237538
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 36
Shape file2: ../datasets/ne_10m_land.csv QQID: 1
PP Polygon size 66475 QQ Polygon size 66475
PP Count 66475 QQ Count 66475
MBR_P [-168.137, -53.886, -34.7936, 72.0021
MBR_Q [-168.137, -53.886, 0, 72.0021
CMBR [-168.137, -53.886, -34.7936, 72.0021

P overlap count with CMBR 66475 Q overlap count with CMBR 66475 

Non-degen count P 0 *****--- Q 0
Intersection count P 66474 *****--- Q 66474

gpuCountIntersections kernel exe time(ms) 504.577057
prefixsum kernels exe time(ms) 2.227008
gpuNeighborMap kernel exe time(ms) 177.351868
gpuCalculateIntersections kernel exe time(ms) 178.300705
gpuSortPolyQ kernel exe time(ms) 0.011232
gpuCalculateInitLabel kernel exe time(ms) 0.016992

Copying completed
Labelling intersections...

... 0 delayed crossings and 1 delayed bouncings
... 1 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 28 vertices removed

R has 1 component with 66447 vertices

All time in microseconds
Time: Total : 991193
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2742
Shape file2: ../datasets/ne_10m_land.csv QQID: 33
PP Polygon size 15547 QQ Polygon size 15547
PP Count 15547 QQ Count 15547
MBR_P [-73.0572, 59.9879, -11.3768, 83.6341
MBR_Q [-73.0572, 59.9879, 0, 83.6341
CMBR [-73.0572, 59.9879, -11.3768, 83.6341

P overlap count with CMBR 15547 Q overlap count with CMBR 15547 

Non-degen count P 4 *****--- Q 4
Intersection count P 15550 *****--- Q 15550

gpuCountIntersections kernel exe time(ms) 36.898079
prefixsum kernels exe time(ms) 0.517280
gpuNeighborMap kernel exe time(ms) 12.572960
gpuCalculateIntersections kernel exe time(ms) 12.898240
gpuSortPolyQ kernel exe time(ms) 0.017408
gpuCalculateInitLabel kernel exe time(ms) 0.006144

Copying completed
Labelling intersections...

... 1 delayed crossings and 2 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 245 vertices removed

R has 2 components with 14751 + 555 = 15306 vertices

All time in microseconds
Time: Total : 193960
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2742 ******************
Shape file2: ../datasets/ne_10m_land.csv QQID: 30
PP Polygon size 15547 QQ Polygon size 10887
PP Count 15547 QQ Count 10887
MBR_P [-73.0572, 59.9879, -11.3768, 83.6341
MBR_Q [-90.1186, 61.8581, 0, 73.8553
CMBR [-73.0572, 61.8581, -11.3768, 73.8553

P overlap count with CMBR 9298 Q overlap count with CMBR 6261 

Non-degen count P 18 *****--- Q 18
Intersection count P 18 *****--- Q 18

gpuCountIntersections kernel exe time(ms) 18.415968
prefixsum kernels exe time(ms) 0.457184
gpuNeighborMap kernel exe time(ms) 0.995328
gpuCalculateIntersections kernel exe time(ms) 0.736000
gpuSortPolyQ kernel exe time(ms) 0.051232
gpuCalculateInitLabel kernel exe time(ms) 0.006784

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 0 vertices removed

R has 5 components with 7 + 4 + 5 + 3 + 18 = 37 vertices

All time in microseconds
Time: Total : 149978
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2742********************
Shape file2: ../datasets/ne_10m_land.csv QQID: 42
PP Polygon size 15547 QQ Polygon size 3838
PP Count 15547 QQ Count 3838
MBR_P [-73.0572, 59.9879, -11.3768, 83.6341
MBR_Q [-91.9566, 76.1291, 0, 83.1165
CMBR [-73.0572, 76.1291, -11.3768, 83.1165

P overlap count with CMBR 3227 Q overlap count with CMBR 594 

Non-degen count P 14 *****--- Q 14
Intersection count P 14 *****--- Q 14

gpuCountIntersections kernel exe time(ms) 5.686816
prefixsum kernels exe time(ms) 0.329408
gpuNeighborMap kernel exe time(ms) 0.992608
gpuCalculateIntersections kernel exe time(ms) 0.264128
gpuSortPolyQ kernel exe time(ms) 0.039776
gpuCalculateInitLabel kernel exe time(ms) 0.008128

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 3 vertices removed

R has 4 components with 6 + 3 + 5 + 7 = 21 vertices

All time in microseconds
Time: Total : 131063
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 5978
Shape file2: ../datasets/ne_10m_land.csv QQID: 3
PP Polygon size 9461 QQ Polygon size 9461
PP Count 9461 QQ Count 9461
MBR_P [113.157, -39.145, 153.631, -10.6877
MBR_Q [0, -39.145, 153.631, -10.6877
CMBR [113.157, -39.145, 153.631, -10.6877

P overlap count with CMBR 9461 Q overlap count with CMBR 9461 

Non-degen count P 10 *****--- Q 10
Intersection count P 9470 *****--- Q 9470

gpuCountIntersections kernel exe time(ms) 16.946207
prefixsum kernels exe time(ms) 0.322272
gpuNeighborMap kernel exe time(ms) 4.997280
gpuCalculateIntersections kernel exe time(ms) 5.239200
gpuSortPolyQ kernel exe time(ms) 0.036032
gpuCalculateInitLabel kernel exe time(ms) 0.006368

Copying completed
Labelling intersections...

... 6 delayed crossings and 3 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 4 vertices removed

R has 3 components with 8836 + 601 + 11 = 9448 vertices

All time in microseconds
Time: Total : 157752
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2854
Shape file2: ../datasets/ne_10m_land.csv QQID: 42
PP Polygon size 3838 QQ Polygon size 3838
PP Count 3838 QQ Count 3838
MBR_P [-91.9566, 76.1291, -61.0855, 83.1165
MBR_Q [-91.9566, 76.1291, 0, 83.1165
CMBR [-91.9566, 76.1291, -61.0855, 83.1165

P overlap count with CMBR 3838 Q overlap count with CMBR 3838 

Non-degen count P 0 *****--- Q 0
Intersection count P 3837 *****--- Q 3837

gpuCountIntersections kernel exe time(ms) 2.565120
prefixsum kernels exe time(ms) 0.136320
gpuNeighborMap kernel exe time(ms) 1.639936
gpuCalculateIntersections kernel exe time(ms) 1.685504
gpuSortPolyQ kernel exe time(ms) 0.004448
gpuCalculateInitLabel kernel exe time(ms) 0.006848

Copying completed
Labelling intersections...

... 0 delayed crossings and 1 delayed bouncings
... 1 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 145 vertices removed

R has 1 component with 3693 vertices

All time in microseconds
Time: Total : 132365
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2737
Shape file2: ../datasets/ne_10m_land.csv QQID: 25
PP Polygon size 3706 QQ Polygon size 3706
PP Count 3706 QQ Count 3706
MBR_P [-6.2364, 49.9591, 1.77117, 58.6771
MBR_Q [-6.2364, 49.9591, 1.77117, 58.6771
CMBR [-6.2364, 49.9591, 1.77117, 58.6771

P overlap count with CMBR 3706 Q overlap count with CMBR 3706 

Non-degen count P 4 *****--- Q 4
Intersection count P 3709 *****--- Q 3709

gpuCountIntersections kernel exe time(ms) 2.460352
prefixsum kernels exe time(ms) 0.125888
gpuNeighborMap kernel exe time(ms) 1.588192
gpuCalculateIntersections kernel exe time(ms) 1.635904
gpuSortPolyQ kernel exe time(ms) 0.017824
gpuCalculateInitLabel kernel exe time(ms) 0.008960

Copying completed
Labelling intersections...

... 3 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 3 vertices removed

R has 1 component with 12 vertices

All time in microseconds
Time: Total : 132952
____________________________________________________________________

======================================================================================
OCEAN vs Continents

======================================================================================


Shape file1: ../datasets/ne_10m_ocean.csv PPID: 0 **************************************
Shape file2: ../datasets/continents.csv QQID: 1048
PP Polygon size 100612 QQ Polygon size 15653
PP Count 100612 QQ Count 15653
MBR_P [-180, -85.2219, 180, 90
MBR_Q [-168.132, 7.20611, 0, 71.9969
CMBR [-168.132, 7.20611, 0, 71.9969

P overlap count with CMBR 6303 Q overlap count with CMBR 15653 

Non-degen count P 0 *****--- Q 0
Intersection count P 0 *****--- Q 0

gpuCountIntersections kernel exe time(ms) 124.958786
prefixsum kernels exe time(ms) 1.975264
gpuNeighborMap kernel exe time(ms) 0.009920
gpuCalculateIntersections kernel exe time(ms) 0.023776
gpuSortPolyQ kernel exe time(ms) 0.006016
gpuCalculateInitLabel kernel exe time(ms) 0.007776

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 1 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 1 vertices removed

R has 1 component with 15652 vertices

All time in microseconds
Time: Total : 260635
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 0---------------------contain
Shape file2: ../datasets/continents.csv QQID: 1202
PP Polygon size 100612 QQ Polygon size 12750
PP Count 100612 QQ Count 12750
MBR_P [-180, -85.2219, 180, 90
MBR_Q [-73.0536, 59.9814, 0, 83.6236
CMBR [-73.0536, 59.9814, 0, 83.6236

P overlap count with CMBR 0 Q overlap count with CMBR 12750 

Non-degen count P 0 *****--- Q 0
Intersection count P 0 *****--- Q 0

gpuCountIntersections kernel exe time(ms) 88.916801
prefixsum kernels exe time(ms) 1.923712
gpuNeighborMap kernel exe time(ms) 0.010592
gpuCalculateIntersections kernel exe time(ms) 0.023808
gpuSortPolyQ kernel exe time(ms) 0.006432
gpuCalculateInitLabel kernel exe time(ms) 0.008128

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 1 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 2 vertices removed

R has 1 component with 12748 vertices

All time in microseco
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 0
Shape file2: ../datasets/continents.csv QQID: 1661
PP Polygon size 100612 QQ Polygon size 12613
PP Count 100612 QQ Count 12613
MBR_P [-180, -85.2219, 180, 90
MBR_Q [-180, -90, 180, -63.2054
CMBR [-180, -85.2219, 180, -63.2054

P overlap count with CMBR 15294 Q overlap count with CMBR 12164 

Non-degen count P 1426 *****--- Q 1426
Intersection count P 1427 *****--- Q 1427

gpuCountIntersections kernel exe time(ms) 101.261826
prefixsum kernels exe time(ms) 1.900000
gpuNeighborMap kernel exe time(ms) 6.549792
gpuCalculateIntersections kernel exe time(ms) 1.474560
gpuSortPolyQ kernel exe time(ms) 0.033600
gpuCalculateInitLabel kernel exe time(ms) 0.012288

Copying completed
Labelling intersections...

... 1 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 222 vertices removed

R has 706 components with 83 + 4 + 4 + 7 + 29 + 8 + 136 + 16 + 4 + 12 + 3 + 13 + 21 + 43 + 40 + 16 + 27 + 12 + 34 + 8 + 59 + 49 + 18 + 20 + 24 + 30 + 8 + 7 + 30 + 3 + 20 + 11 + 11 + 14 + 4 + 21 + 5 + 30 + 25 + 24 + 3 + 15 + 8 + 24 + 9 + 5 + 7 + 93 + 7 + 12 + 80 + 157 + 321 + 69 + 190 + 88 + 13 + 26 + 3 + 36 + 19 + 15 + 39 + 58 + 4 + 16 + 36 + 10 + 97 + 4 + 28 + 3 + 28 + 18 + 16 + 23 + 55 + 59 + 14 + 9 + 21 + 7 + 8 + 6 + 5 + 72 + 24 + 11 + 12 + 59 + 51 + 35 + 30 + 6 + 44 + 47 + 7 + 21 + 90 + 68 + 9 + 41 + 36 + 190 + 50 + 20 + 140 + 82 + 19 + 132 + 7 + 10 + 9 + 18 + 42 + 5 + 13 + 12 + 58 + 11 + 20 + 22 + 14 + 11 + 4 + 5 + 8 + 6 + 15 + 5 + 20 + 13 + 17 + 8 + 23 + 13 + 3 + 6 + 8 + 32 + 18 + 9 + 6 + 13 + 3 + 18 + 35 + 14 + 9 + 8 + 8 + 3 + 5 + 25 + 38 + 26 + 18 + 31 + 9 + 15 + 39 + 21 + 20 + 29 + 6 + 8 + 5 + 6 + 11 + 19 + 18 + 4 + 6 + 11 + 13 + 3 + 34 + 3 + 5 + 20 + 6 + 3 + 12 + 22 + 30 + 6 + 7 + 18 + 3 + 11 + 9 + 12 + 21 + 29 + 61 + 12 + 25 + 10 + 3 + 16 + 7 + 12 + 10 + 6 + 4 + 9 + 10 + 6 + 5 + 11 + 10 + 8 + 7 + 22 + 3 + 3 + 16 + 26 + 13 + 11 + 19 + 10 + 8 + 4 + 5 + 5 + 6 + 7 + 11 + 8 + 10 + 5 + 6 + 10 + 9 + 4 + 14 + 3 + 10 + 10 + 5 + 10 + 24 + 7 + 14 + 24 + 7 + 4 + 8 + 4 + 4 + 17 + 9 + 9 + 18 + 10 + 9 + 5 + 5 + 21 + 8 + 10 + 8 + 8 + 9 + 22 + 9 + 7 + 6 + 21 + 15 + 7 + 13 + 16 + 7 + 12 + 20 + 9 + 6 + 33 + 23 + 12 + 15 + 5 + 13 + 7 + 5 + 6 + 9 + 17 + 19 + 4 + 23 + 14 + 9 + 24 + 5 + 7 + 8 + 19 + 9 + 7 + 8 + 11 + 13 + 5 + 6 + 17 + 15 + 7 + 8 + 6 + 65 + 4 + 8 + 8 + 3 + 18 + 3 + 4 + 4 + 15 + 10 + 3 + 5 + 3 + 7 + 3 + 4 + 7 + 10 + 5 + 10 + 73 + 5 + 90 + 8 + 8 + 9 + 25 + 24 + 7 + 9 + 9 + 8 + 4 + 4 + 8 + 6 + 6 + 12 + 5 + 12 + 3 + 8 + 10 + 3 + 5 + 6 + 10 + 28 + 6 + 14 + 20 + 4 + 8 + 11 + 9 + 4 + 10 + 4 + 30 + 10 + 14 + 29 + 9 + 3 + 10 + 3 + 19 + 6 + 14 + 8 + 9 + 29 + 14 + 41 + 14 + 31 + 12 + 54 + 40 + 5 + 43 + 17 + 21 + 10 + 14 + 7 + 6 + 5 + 9 + 3 + 23 + 7 + 9 + 24 + 6 + 4 + 13 + 132 + 19 + 30 + 23 + 3 + 8 + 13 + 6 + 5 + 42 + 6 + 83 + 103 + 113 + 18 + 42 + 9 + 121 + 37 + 12 + 8 + 5 + 94 + 38 + 20 + 38 + 8 + 34 + 3 + 22 + 24 + 74 + 103 + 77 + 60 + 20 + 14 + 29 + 24 + 8 + 38 + 17 + 19 + 15 + 8 + 3 + 14 + 7 + 4 + 8 + 12 + 8 + 7 + 5 + 7 + 10 + 26 + 36 + 69 + 9 + 18 + 19 + 81 + 5 + 43 + 12 + 27 + 49 + 27 + 6 + 7 + 3 + 16 + 6 + 22 + 9 + 28 + 25 + 7 + 15 + 20 + 12 + 7 + 16 + 6 + 12 + 20 + 29 + 12 + 21 + 8 + 21 + 15 + 14 + 5 + 17 + 13 + 5 + 13 + 9 + 4 + 30 + 31 + 8 + 21 + 11 + 12 + 55 + 6 + 3 + 5 + 15 + 7 + 17 + 3 + 9 + 29 + 12 + 32 + 9 + 10 + 12 + 15 + 41 + 13 + 23 + 63 + 19 + 17 + 3 + 11 + 32 + 4 + 13 + 32 + 7 + 4 + 42 + 72 + 55 + 18 + 4 + 15 + 19 + 65 + 20 + 159 + 3 + 178 + 157 + 30 + 5 + 108 + 167 + 14 + 84 + 33 + 12 + 32 + 3 + 39 + 12 + 10 + 38 + 5 + 6 + 3 + 9 + 36 + 23 + 53 + 48 + 20 + 4 + 20 + 3 + 5 + 16 + 5 + 60 + 4 + 3 + 7 + 32 + 7 + 143 + 35 + 46 + 10 + 9 + 4 + 15 + 17 + 10 + 5 + 8 + 18 + 5 + 170 + 88 + 16 + 142 + 21 + 11 + 50 + 24 + 11 + 103 + 8 + 48 + 210 + 28 + 15 + 17 + 5 + 46 + 17 + 27 + 5 + 54 + 49 + 50 + 5 + 80 + 378 + 6 + 101 + 11 + 5 + 226 + 49 + 14 + 63 + 15 + 6 + 29 + 11 + 10 + 11 + 8 + 7 + 39 + 12 + 141 + 7 + 4 + 7 + 240 + 3 + 17 + 20 + 30 + 23 + 116 + 3 + 30 + 10 + 28 + 9 + 7 + 5 + 16 + 3 + 17 + 16 + 12 + 11 + 10 + 4 + 195 + 22 + 5 + 50 + 3 + 3 + 62 + 18 + 9 + 6 + 8 + 17 + 24 + 19 + 55 + 4 + 26 + 19 + 16 + 5 + 21 + 18 + 6 + 9 + 12 + 56 + 4 = 16895 vertices

All time in microseconds
Time: Total : 261048
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 0
Shape file2: ../datasets/continents.csv QQID: 1886
PP Polygon size 100612 QQ Polygon size 11398
PP Count 100612 QQ Count 11398
MBR_P [-180, -85.2219, 180, 90
MBR_Q [-9.49083, 36.0061, 66.2105, 71.1131
CMBR [-9.49083, 36.0061, 66.2105, 71.1131

P overlap count with CMBR 31072 Q overlap count with CMBR 11398 

Non-degen count P 7162 *****--- Q 7162
Intersection count P 7162 *****--- Q 7162

gpuCountIntersections kernel exe time(ms) 110.923073
prefixsum kernels exe time(ms) 1.878304
gpuNeighborMap kernel exe time(ms) 9.418592
gpuCalculateIntersections kernel exe time(ms) 7.989856
gpuSortPolyQ kernel exe time(ms) 0.049152
gpuCalculateInitLabel kernel exe time(ms) 0.023392

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 26 vertices removed

R has 3563 components with 3 + 3 + 6 + 7 + 3 + 7 + 3 + 7 + 5 + 7 + 5 + 3 + 3 + 9 + 6 + 5 + 4 + 8 + 11 + 9 + 3 + 5 + 5 + 12 + 4 + 6 + 6 + 7 + 3 + 3 + 3 + 6 + 7 + 4 + 3 + 3 + 6 + 11 + 14 + 13 + 5 + 5 + 21 + 3 + 8 + 4 + 7 + 8 + 3 + 6 + 5 + 3 + 6 + 4 + 8 + 5 + 3 + 5 + 3 + 3 + 6 + 14 + 3 + 19 + 3 + 6 + 10 + 15 + 3 + 6 + 5 + 14 + 3 + 3 + 6 + 4 + 5 + 3 + 6 + 17 + 18 + 3 + 18 + 14 + 3 + 4 + 4 + 3 + 5 + 8 + 3 + 4 + 3 + 4 + 6 + 4 + 4 + 8 + 3 + 3 + 3 + 3 + 8 + 4 + 3 + 4 + 3 + 3 + 4 + 3 + 6 + 9 + 3 + 3 + 12 + 7 + 8 + 4 + 4 + 4 + 20 + 4 + 3 + 11 + 6 + 3 + 6 + 7 + 5 + 7 + 7 + 5 + 6 + 4 + 3 + 19 + 4 + 3 + 31 + 9 + 4 + 18 + 8 + 3 + 3 + 6 + 4 + 3 + 11 + 13 + 24 + 5 + 16 + 9 + 3 + 4 + 13 + 4 + 4 + 10 + 3 + 11 + 6 + 7 + 3 + 3 + 3 + 3 + 8 + 10 + 3 + 6 + 12 + 3 + 10 + 5 + 8 + 3 + 4 + 5 + 5 + 3 + 20 + 14 + 6 + 3 + 10 + 15 + 3 + 3 + 4 + 3 + 12 + 9 + 20 + 4 + 4 + 7 + 3 + 15 + 50 + 3 + 3 + 5 + 12 + 4 + 9 + 8 + 7 + 9 + 7 + 3 + 9 + 5 + 7 + 22 + 4 + 4 + 8 + 15 + 3 + 3 + 4 + 5 + 16 + 5 + 6 + 9 + 13 + 9 + 11 + 5 + 3 + 15 + 6 + 3 + 8 + 4 + 3 + 4 + 16 + 18 + 18 + 3 + 4 + 3 + 11 + 15 + 6 + 8 + 6 + 3 + 8 + 4 + 3 + 3 + 3 + 3 + 16 + 13 + 3 + 6 + 8 + 6 + 26 + 7 + 7 + 18 + 11 + 3 + 3 + 16 + 3 + 15 + 14 + 3 + 30 + 3 + 3 + 4 + 4 + 6 + 34 + 20 + 5 + 13 + 23 + 6 + 5 + 10 + 4 + 6 + 28 + 9 + 17 + 11 + 3 + 14 + 8 + 8 + 7 + 4 + 8 + 3 + 20 + 3 + 3 + 4 + 8 + 8 + 6 + 3 + 5 + 10 + 6 + 4 + 6 + 11 + 6 + 3 + 6 + 8 + 6 + 6 + 4 + 5 + 6 + 6 + 3 + 6 + 3 + 4 + 7 + 4 + 4 + 7 + 3 + 3 + 16 + 9 + 5 + 3 + 6 + 3 + 4 + 19 + 14 + 4 + 7 + 6 + 5 + 3 + 5 + 4 + 3 + 3 + 7 + 3 + 5 + 5 + 7 + 3 + 14 + 7 + 5 + 9 + 25 + 23 + 8 + 6 + 4 + 4 + 6 + 3 + 10 + 11 + 3 + 3 + 14 + 19 + 5 + 3 + 11 + 5 + 3 + 7 + 5 + 14 + 7 + 5 + 11 + 9 + 10 + 6 + 3 + 8 + 5 + 5 + 6 + 9 + 3 + 12 + 5 + 3 + 4 + 7 + 3 + 7 + 14 + 10 + 3 + 5 + 5 + 6 + 8 + 3 + 8 + 5 + 4 + 6 + 10 + 7 + 3 + 3 + 11 + 3 + 4 + 12 + 5 + 4 + 7 + 4 + 6 + 3 + 3 + 12 + 5 + 27 + 7 + 13 + 3 + 3 + 12 + 3 + 3 + 15 + 3 + 11 + 3 + 27 + 3 + 6 + 10 + 4 + 8 + 9 + 40 + 8 + 19 + 10 + 5 + 4 + 19 + 5 + 12 + 18 + 5 + 3 + 4 + 3 + 4 + 12 + 5 + 6 + 6 + 11 + 4 + 4 + 3 + 6 + 4 + 7 + 3 + 3 + 3 + 10 + 4 + 20 + 3 + 5 + 4 + 3 + 8 + 3 + 4 + 3 + 4 + 17 + 9 + 7 + 7 + 26 + 7 + 3 + 12 + 18 + 3 + 4 + 5 + 3 + 3 + 8 + 7 + 19 + 5 + 3 + 3 + 5 + 7 + 9 + 16 + 6 + 3 + 3 + 5 + 4 + 10 + 5 + 6 + 3 + 3 + 14 + 10 + 3 + 4 + 4 + 10 + 8 + 3 + 6 + 7 + 4 + 5 + 10 + 9 + 6 + 4 + 5 + 16 + 5 + 4 + 4 + 9 + 3 + 8 + 7 + 8 + 6 + 23 + 4 + 21 + 20 + 23 + 28 + 5 + 4 + 3 + 3 + 4 + 4 + 3 + 18 + 7 + 4 + 3 + 3 + 7 + 4 + 10 + 5 + 4 + 3 + 12 + 5 + 6 + 5 + 25 + 4 + 7 + 4 + 4 + 3 + 4 + 15 + 4 + 3 + 5 + 4 + 19 + 5 + 4 + 10 + 3 + 3 + 4 + 7 + 3 + 4 + 3 + 4 + 14 + 6 + 5 + 4 + 3 + 18 + 7 + 7 + 7 + 3 + 13 + 5 + 3 + 4 + 8 + 3 + 7 + 3 + 31 + 6 + 3 + 7 + 16 + 12 + 6 + 4 + 6 + 10 + 3 + 3 + 8 + 11 + 3 + 3 + 8 + 5 + 3 + 4 + 9 + 3 + 3 + 3 + 14 + 3 + 3 + 3 + 4 + 27 + 3 + 3 + 3 + 10 + 3 + 15 + 12 + 10 + 3 + 13 + 6 + 7 + 8 + 6 + 17 + 5 + 6 + 7 + 8 + 5 + 3 + 5 + 3 + 3 + 3 + 14 + 7 + 7 + 3 + 5 + 3 + 5 + 3 + 3 + 13 + 3 + 9 + 8 + 9 + 6 + 5 + 8 + 7 + 5 + 6 + 4 + 3 + 22 + 14 + 5 + 4 + 3 + 4 + 7 + 14 + 12 + 5 + 4 + 9 + 4 + 6 + 3 + 8 + 3 + 5 + 37 + 7 + 7 + 8 + 23 + 3 + 7 + 3 + 5 + 9 + 3 + 4 + 3 + 13 + 6 + 6 + 4 + 6 + 53 + 4 + 4 + 9 + 3 + 4 + 5 + 31 + 15 + 6 + 4 + 8 + 9 + 4 + 3 + 4 + 10 + 21 + 12 + 12 + 3 + 6 + 14 + 7 + 3 + 4 + 7 + 4 + 3 + 19 + 3 + 5 + 6 + 12 + 3 + 3 + 6 + 8 + 4 + 6 + 4 + 3 + 5 + 4 + 3 + 9 + 3 + 6 + 3 + 3 + 7 + 10 + 5 + 6 + 15 + 8 + 13 + 6 + 4 + 3 + 3 + 3 + 3 + 3 + 24 + 5 + 3 + 14 + 3 + 12 + 3 + 5 + 12 + 6 + 4 + 15 + 19 + 3 + 3 + 13 + 16 + 11 + 7 + 9 + 3 + 43 + 3 + 3 + 3 + 9 + 26 + 4 + 6 + 5 + 4 + 8 + 21 + 3 + 6 + 3 + 5 + 4 + 3 + 8 + 11 + 5 + 13 + 4 + 4 + 3 + 4 + 4 + 3 + 36 + 6 + 6 + 3 + 4 + 5 + 10 + 13 + 30 + 10 + 7 + 17 + 8 + 4 + 6 + 4 + 6 + 11 + 3 + 11 + 9 + 12 + 5 + 4 + 18 + 5 + 3 + 6 + 3 + 7 + 7 + 7 + 14 + 3 + 4 + 6 + 4 + 3 + 4 + 9 + 9 + 9 + 10 + 5 + 8 + 4 + 6 + 3 + 3 + 10 + 7 + 5 + 4 + 13 + 7 + 3 + 5 + 27 + 6 + 6 + 9 + 4 + 10 + 5 + 3 + 5 + 8 + 4 + 25 + 16 + 9 + 3 + 3 + 11 + 13 + 3 + 4 + 8 + 3 + 16 + 4 + 17 + 25 + 7 + 3 + 5 + 5 + 6 + 3 + 3 + 9 + 6 + 3 + 11 + 16 + 12 + 4 + 29 + 12 + 23 + 19 + 3 + 4 + 3 + 4 + 8 + 19 + 4 + 8 + 4 + 5 + 5 + 3 + 8 + 6 + 8 + 4 + 5 + 5 + 10 + 4 + 4 + 4 + 7 + 3 + 3 + 7 + 4 + 6 + 3 + 4 + 8 + 9 + 4 + 3 + 4 + 9 + 3 + 6 + 4 + 3 + 12 + 7 + 6 + 3 + 35 + 3 + 11 + 3 + 14 + 6 + 4 + 6 + 4 + 3 + 8 + 9 + 9 + 19 + 27 + 11 + 5 + 3 + 17 + 5 + 3 + 6 + 4 + 3 + 6 + 22 + 15 + 5 + 4 + 3 + 4 + 6 + 3 + 3 + 4 + 3 + 5 + 4 + 3 + 4 + 4 + 13 + 7 + 12 + 4 + 3 + 6 + 11 + 12 + 7 + 6 + 9 + 20 + 4 + 3 + 22 + 19 + 7 + 4 + 5 + 8 + 3 + 9 + 5 + 17 + 7 + 3 + 10 + 7 + 5 + 5 + 4 + 4 + 3 + 16 + 3 + 7 + 8 + 3 + 5 + 3 + 12 + 3 + 7 + 11 + 8 + 3 + 3 + 15 + 3 + 11 + 4 + 6 + 8 + 4 + 14 + 5 + 14 + 3 + 5 + 7 + 10 + 6 + 17 + 5 + 6 + 3 + 18 + 13 + 4 + 14 + 7 + 3 + 3 + 29 + 4 + 14 + 12 + 9 + 4 + 3 + 3 + 13 + 23 + 7 + 5 + 3 + 29 + 3 + 4 + 10 + 4 + 3 + 6 + 23 + 5 + 4 + 3 + 3 + 32 + 6 + 7 + 6 + 17 + 5 + 4 + 3 + 16 + 8 + 3 + 6 + 13 + 18 + 21 + 42 + 8 + 7 + 8 + 24 + 5 + 7 + 11 + 5 + 4 + 6 + 10 + 3 + 14 + 4 + 10 + 19 + 3 + 12 + 3 + 7 + 7 + 14 + 10 + 6 + 12 + 3 + 7 + 11 + 3 + 29 + 18 + 11 + 8 + 6 + 3 + 3 + 5 + 4 + 8 + 4 + 3 + 3 + 5 + 5 + 3 + 5 + 7 + 4 + 3 + 16 + 6 + 13 + 12 + 3 + 3 + 3 + 6 + 7 + 11 + 3 + 3 + 3 + 25 + 8 + 10 + 3 + 6 + 12 + 5 + 6 + 37 + 4 + 13 + 16 + 12 + 16 + 17 + 19 + 11 + 13 + 4 + 4 + 7 + 6 + 5 + 7 + 5 + 14 + 9 + 5 + 17 + 5 + 5 + 6 + 44 + 4 + 3 + 4 + 4 + 21 + 4 + 3 + 3 + 9 + 16 + 11 + 4 + 14 + 4 + 6 + 7 + 44 + 16 + 15 + 4 + 4 + 11 + 18 + 8 + 8 + 4 + 32 + 12 + 6 + 13 + 43 + 3 + 3 + 10 + 4 + 16 + 3 + 5 + 5 + 5 + 5 + 4 + 7 + 5 + 4 + 8 + 16 + 8 + 67 + 14 + 9 + 6 + 4 + 5 + 4 + 30 + 10 + 32 + 10 + 61 + 3 + 7 + 5 + 5 + 5 + 3 + 16 + 17 + 25 + 4 + 3 + 39 + 7 + 5 + 5 + 6 + 5 + 9 + 5 + 3 + 6 + 10 + 9 + 3 + 3 + 10 + 4 + 10 + 3 + 9 + 3 + 4 + 10 + 10 + 4 + 4 + 3 + 6 + 9 + 3 + 25 + 8 + 3 + 4 + 4 + 13 + 3 + 4 + 5 + 11 + 7 + 4 + 13 + 51 + 6 + 6 + 4 + 3 + 8 + 3 + 5 + 3 + 4 + 9 + 8 + 4 + 7 + 23 + 42 + 36 + 9 + 9 + 3 + 5 + 4 + 5 + 10 + 4 + 4 + 3 + 5 + 15 + 9 + 6 + 9 + 7 + 3 + 12 + 5 + 4 + 3 + 3 + 11 + 5 + 8 + 33 + 11 + 5 + 5 + 17 + 4 + 3 + 6 + 3 + 3 + 3 + 6 + 4 + 4 + 8 + 10 + 3 + 5 + 3 + 5 + 4 + 5 + 10 + 5 + 15 + 5 + 27 + 8 + 9 + 5 + 7 + 4 + 3 + 18 + 3 + 8 + 8 + 6 + 8 + 16 + 10 + 14 + 8 + 11 + 8 + 12 + 3 + 7 + 11 + 3 + 7 + 4 + 5 + 32 + 6 + 4 + 7 + 7 + 10 + 48 + 5 + 6 + 3 + 4 + 5 + 4 + 4 + 11 + 18 + 13 + 3 + 6 + 5 + 5 + 7 + 8 + 55 + 27 + 3 + 14 + 3 + 13 + 5 + 5 + 4 + 5 + 4 + 9 + 3 + 6 + 6 + 16 + 7 + 8 + 3 + 3 + 4 + 5 + 27 + 3 + 7 + 3 + 8 + 3 + 6 + 4 + 14 + 4 + 3 + 3 + 6 + 4 + 6 + 6 + 10 + 15 + 3 + 9 + 4 + 3 + 6 + 3 + 5 + 3 + 3 + 4 + 7 + 10 + 9 + 3 + 3 + 3 + 4 + 19 + 9 + 4 + 4 + 6 + 3 + 9 + 9 + 14 + 13 + 4 + 13 + 12 + 4 + 12 + 3 + 4 + 3 + 8 + 7 + 5 + 3 + 4 + 5 + 3 + 4 + 12 + 4 + 3 + 3 + 4 + 3 + 8 + 4 + 7 + 14 + 9 + 5 + 3 + 5 + 4 + 4 + 11 + 5 + 9 + 3 + 17 + 7 + 5 + 4 + 14 + 9 + 7 + 7 + 6 + 5 + 8 + 8 + 15 + 5 + 7 + 3 + 3 + 3 + 5 + 4 + 3 + 7 + 10 + 3 + 34 + 3 + 4 + 7 + 6 + 6 + 13 + 3 + 5 + 5 + 6 + 13 + 15 + 3 + 71 + 4 + 6 + 7 + 5 + 5 + 9 + 13 + 3 + 11 + 10 + 4 + 14 + 12 + 5 + 8 + 7 + 11 + 3 + 8 + 7 + 4 + 14 + 6 + 4 + 3 + 10 + 29 + 21 + 5 + 8 + 6 + 9 + 34 + 36 + 6 + 4 + 14 + 5 + 3 + 3 + 7 + 7 + 5 + 7 + 3 + 22 + 3 + 7 + 13 + 8 + 10 + 3 + 4 + 3 + 3 + 7 + 10 + 7 + 11 + 4 + 3 + 36 + 4 + 12 + 10 + 3 + 3 + 19 + 14 + 9 + 4 + 10 + 3 + 7 + 3 + 7 + 4 + 4 + 7 + 12 + 13 + 9 + 45 + 9 + 11 + 3 + 7 + 4 + 8 + 5 + 28 + 13 + 4 + 10 + 12 + 4 + 27 + 3 + 5 + 6 + 12 + 5 + 10 + 4 + 3 + 4 + 5 + 17 + 3 + 3 + 5 + 4 + 14 + 7 + 8 + 8 + 5 + 3 + 4 + 6 + 4 + 6 + 7 + 5 + 3 + 8 + 7 + 5 + 12 + 4 + 9 + 14 + 7 + 5 + 5 + 7 + 5 + 11 + 8 + 10 + 12 + 4 + 4 + 11 + 7 + 3 + 30 + 3 + 3 + 5 + 5 + 4 + 3 + 3 + 9 + 14 + 3 + 17 + 8 + 8 + 6 + 4 + 3 + 3 + 4 + 3 + 4 + 16 + 3 + 3 + 13 + 6 + 7 + 4 + 10 + 7 + 12 + 10 + 10 + 15 + 10 + 6 + 8 + 14 + 13 + 6 + 7 + 5 + 11 + 6 + 6 + 4 + 9 + 4 + 6 + 6 + 9 + 3 + 6 + 3 + 3 + 13 + 12 + 7 + 13 + 6 + 11 + 5 + 22 + 5 + 21 + 4 + 8 + 4 + 3 + 9 + 5 + 21 + 7 + 8 + 8 + 9 + 7 + 17 + 29 + 15 + 3 + 4 + 6 + 7 + 7 + 6 + 4 + 11 + 3 + 13 + 7 + 6 + 4 + 4 + 9 + 4 + 3 + 42 + 11 + 8 + 3 + 11 + 6 + 3 + 3 + 7 + 3 + 15 + 8 + 6 + 6 + 8 + 4 + 3 + 5 + 5 + 3 + 4 + 7 + 4 + 4 + 9 + 8 + 14 + 12 + 10 + 5 + 5 + 4 + 3 + 3 + 10 + 4 + 4 + 3 + 11 + 13 + 20 + 16 + 3 + 8 + 9 + 3 + 6 + 6 + 3 + 5 + 5 + 3 + 5 + 3 + 4 + 3 + 5 + 3 + 3 + 6 + 10 + 11 + 4 + 9 + 3 + 4 + 9 + 4 + 4 + 8 + 3 + 3 + 3 + 5 + 12 + 6 + 6 + 3 + 0 + 7 + 5 + 3 + 7 + 3 + 7 + 3 + 15 + 9 + 5 + 7 + 10 + 5 + 4 + 9 + 11 + 11 + 7 + 28 + 5 + 3 + 3 + 7 + 40 + 3 + 4 + 9 + 5 + 4 + 4 + 4 + 9 + 4 + 3 + 3 + 14 + 3 + 4 + 4 + 6 + 3 + 4 + 6 + 3 + 3 + 5 + 5 + 3 + 7 + 6 + 6 + 14 + 3 + 11 + 8 + 3 + 6 + 3 + 14 + 5 + 3 + 5 + 7 + 4 + 3 + 8 + 3 + 4 + 11 + 15 + 5 + 9 + 5 + 4 + 3 + 3 + 3 + 4 + 3 + 6 + 3 + 16 + 3 + 3 + 5 + 3 + 5 + 25 + 11 + 3 + 4 + 4 + 5 + 3 + 5 + 11 + 7 + 3 + 10 + 4 + 8 + 17 + 11 + 10 + 18 + 5 + 3 + 3 + 3 + 4 + 12 + 8 + 15 + 25 + 9 + 3 + 18 + 5 + 7 + 4 + 6 + 5 + 7 + 3 + 3 + 3 + 11 + 3 + 3 + 4 + 8 + 6 + 3 + 3 + 4 + 5 + 6 + 4 + 20 + 3 + 3 + 4 + 17 + 5 + 5 + 4 + 4 + 6 + 8 + 3 + 4 + 3 + 5 + 3 + 4 + 9 + 14 + 3 + 4 + 3 + 8 + 20 + 15 + 22 + 19 + 6 + 9 + 8 + 9 + 13 + 3 + 15 + 10 + 4 + 6 + 6 + 7 + 6 + 3 + 5 + 3 + 4 + 4 + 28 + 7 + 9 + 14 + 4 + 7 + 3 + 4 + 10 + 3 + 29 + 6 + 4 + 8 + 4 + 8 + 3 + 23 + 3 + 16 + 15 + 3 + 3 + 3 + 4 + 3 + 8 + 5 + 7 + 7 + 4 + 7 + 4 + 6 + 10 + 3 + 15 + 4 + 4 + 3 + 28 + 3 + 6 + 3 + 4 + 4 + 4 + 3 + 3 + 4 + 4 + 3 + 3 + 4 + 3 + 3 + 4 + 3 + 11 + 10 + 3 + 4 + 7 + 5 + 3 + 10 + 9 + 3 + 4 + 9 + 5 + 5 + 8 + 4 + 6 + 4 + 3 + 5 + 3 + 3 + 3 + 5 + 30 + 3 + 4 + 3 + 7 + 5 + 3 + 3 + 14 + 5 + 36 + 8 + 12 + 3 + 16 + 3 + 16 + 11 + 6 + 7 + 8 + 7 + 11 + 6 + 8 + 4 + 25 + 4 + 4 + 4 + 7 + 4 + 3 + 11 + 5 + 3 + 3 + 6 + 6 + 6 + 6 + 3 + 3 + 12 + 4 + 3 + 5 + 3 + 5 + 4 + 3 + 15 + 6 + 3 + 3 + 9 + 5 + 3 + 7 + 5 + 7 + 3 + 3 + 3 + 4 + 3 + 3 + 3 + 5 + 5 + 4 + 8 + 3 + 9 + 4 + 4 + 5 + 4 + 6 + 9 + 6 + 3 + 3 + 4 + 13 + 7 + 3 + 5 + 3 + 5 + 3 + 33 + 18 + 3 + 3 + 7 + 4 + 7 + 3 + 3 + 3 + 4 + 4 + 6 + 3 + 10 + 5 + 26 + 3 + 4 + 5 + 5 + 3 + 3 + 3 + 14 + 4 + 3 + 3 + 3 + 3 + 5 + 19 + 4 + 8 + 5 + 4 + 3 + 5 + 9 + 4 + 3 + 5 + 4 + 4 + 7 + 10 + 10 + 3 + 16 + 3 + 10 + 20 + 15 + 8 + 3 + 15 + 6 + 3 + 3 + 4 + 3 + 3 + 6 + 7 + 3 + 3 + 37 + 6 + 3 + 5 + 7 + 3 + 6 + 14 + 15 + 21 + 4 + 3 + 4 + 3 + 4 + 6 + 19 + 3 + 3 + 5 + 12 + 14 + 11 + 9 + 4 + 3 + 4 + 17 + 28 + 28 + 5 + 9 + 10 + 12 + 3 + 10 + 5 + 5 + 3 + 3 + 3 + 4 + 6 + 8 + 13 + 6 + 10 + 24 + 16 + 24 + 29 + 43 + 15 + 7 + 3 + 9 + 3 + 6 + 4 + 40 + 6 + 3 + 10 + 16 + 4 + 3 + 4 + 4 + 11 + 18 + 4 + 27 + 4 + 3 + 3 + 4 + 3 + 6 + 11 + 7 + 6 + 6 + 8 + 3 + 5 + 5 + 3 + 17 + 6 + 9 + 3 + 3 + 3 + 6 + 3 + 5 + 7 + 14 + 8 + 5 + 5 + 3 + 4 + 3 + 10 + 3 + 4 + 3 + 3 + 4 + 5 + 4 + 4 + 8 + 8 + 3 + 18 + 3 + 3 + 4 + 8 + 8 + 4 + 5 + 11 + 3 + 6 + 4 + 11 + 6 + 7 + 3 + 3 + 27 + 8 + 4 + 33 + 3 + 11 + 6 + 7 + 5 + 6 + 4 + 4 + 13 + 6 + 3 + 5 + 5 + 18 + 4 + 8 + 5 + 17 + 12 + 27 + 6 + 5 + 13 + 40 + 10 + 11 + 4 + 10 + 7 + 23 + 10 + 4 + 7 + 3 + 9 + 9 + 25 + 3 + 15 + 3 + 6 + 4 + 7 + 3 + 3 + 3 + 12 + 11 + 5 + 4 + 3 + 20 + 4 + 4 + 7 + 3 + 4 + 15 + 6 + 7 + 16 + 7 + 11 + 10 + 15 + 14 + 9 + 6 + 3 + 3 + 11 + 4 + 4 + 9 + 8 + 0 + 4 + 6 + 6 + 28 + 3 + 7 + 4 + 6 + 6 + 16 + 7 + 3 + 4 + 3 + 9 + 3 + 7 + 8 + 4 + 3 + 3 + 4 + 8 + 9 + 4 + 5 + 16 + 4 + 22 + 14 + 5 + 10 + 9 + 6 + 16 + 13 + 8 + 5 + 12 + 4 + 7 + 3 + 21 + 17 + 4 + 3 + 5 + 8 + 11 + 3 + 19 + 14 + 17 + 5 + 15 + 4 + 15 + 12 + 7 + 5 + 26 + 3 + 4 + 5 + 3 + 8 + 4 + 5 + 4 + 7 + 16 + 8 + 3 + 4 + 5 + 4 + 3 + 7 + 4 + 4 + 5 + 3 + 3 + 5 + 14 + 6 + 3 + 4 + 6 + 11 + 3 + 4 + 12 + 3 + 4 + 7 + 8 + 3 + 3 + 5 + 3 + 5 + 12 + 7 + 4 + 3 + 7 + 8 + 4 + 6 + 4 + 3 + 3 + 3 + 8 + 3 + 3 + 6 + 3 + 12 + 6 + 3 + 4 + 6 + 4 + 13 + 6 + 6 + 3 + 10 + 4 + 7 + 7 + 14 + 3 + 3 + 37 + 5 + 15 + 6 + 4 + 4 + 3 + 8 + 5 + 3 + 5 + 3 + 5 + 4 + 8 + 11 + 15 + 6 + 5 + 20 + 4 + 3 + 13 + 11 + 14 + 4 + 3 + 3 + 4 + 13 + 8 + 4 + 5 + 4 + 6 + 3 + 10 + 7 + 5 + 3 + 8 + 11 + 5 + 14 + 4 + 4 + 6 + 4 + 4 + 3 + 14 + 9 + 9 + 6 + 6 + 9 + 3 + 8 + 5 + 4 + 8 + 9 + 8 + 6 + 6 + 4 + 4 + 12 + 6 + 4 + 11 + 6 + 7 + 14 + 4 + 6 + 5 + 13 + 9 + 3 + 3 + 18 + 8 + 9 + 4 + 13 + 3 + 3 + 4 + 3 + 3 + 8 + 5 + 7 + 5 + 7 + 3 + 7 + 4 + 3 + 3 + 6 + 14 + 3 + 4 + 4 + 3 + 7 + 3 + 4 + 4 + 7 + 5 + 3 + 16 + 8 + 6 + 4 + 17 + 12 + 6 + 4 + 6 + 4 + 3 + 8 + 3 + 3 + 5 + 4 + 5 + 3 + 5 + 3 + 3 + 4 + 5 + 3 + 3 + 3 + 3 + 3 + 16 + 3 + 5 + 7 + 4 + 3 + 3 + 4 + 4 + 4 + 4 + 7 + 3 + 5 + 6 + 4 + 3 + 5 + 3 + 3 + 12 + 3 + 11 + 4 + 3 + 3 + 13 + 8 + 13 + 14 + 6 + 7 + 4 + 3 + 4 + 3 + 4 + 3 + 7 + 13 + 6 + 6 + 6 + 4 + 6 + 3 + 8 + 17 + 5 + 5 + 5 + 3 + 3 + 4 + 15 + 4 + 7 + 3 + 3 + 11 + 7 + 21 + 3 + 6 + 4 + 8 + 31 + 14 + 6 + 6 + 12 + 4 + 4 + 7 + 23 + 9 + 10 + 34 + 8 + 5 + 5 + 5 + 4 + 39 + 5 + 5 + 3 + 6 + 4 + 3 + 3 + 6 + 6 + 12 + 14 + 9 + 9 + 25 + 4 + 7 + 16 + 4 + 5 + 3 + 17 + 4 + 17 + 3 + 17 + 3 + 4 + 26 + 7 + 20 + 11 + 12 + 4 + 8 + 4 + 6 + 6 + 3 + 3 + 6 + 9 + 10 + 5 + 8 + 4 + 4 + 5 + 3 + 4 + 28 + 5 + 3 + 0 + 4 + 12 + 8 + 13 + 6 + 3 + 9 + 3 + 4 + 3 + 4 + 5 + 3 + 4 + 3 + 3 + 8 + 5 + 4 + 7 + 4 + 7 + 3 + 9 + 5 + 3 + 6 + 8 + 6 + 9 + 6 + 7 + 23 + 9 + 4 + 4 + 11 + 10 + 7 + 4 + 3 + 10 + 6 + 6 + 5 + 3 + 17 + 7 + 8 + 5 + 15 + 19 + 3 + 7 + 6 + 5 + 3 + 7 + 3 + 4 + 5 + 3 + 50 + 5 + 3 + 5 + 7 + 6 + 5 + 8 + 5 + 19 + 4 + 5 + 6 + 3 + 28 + 11 + 5 + 8 + 3 + 4 + 3 + 4 + 5 + 4 + 6 + 3 + 3 + 3 + 3 + 4 + 4 + 4 + 4 + 4 + 5 + 7 + 5 + 11 + 3 + 3 + 4 + 3 + 6 + 5 + 7 + 6 + 4 + 4 + 8 + 3 + 4 + 5 + 6 + 13 + 11 + 8 + 3 + 5 + 20 + 5 + 12 + 8 + 13 + 3 + 5 + 4 + 3 + 3 + 3 + 3 + 8 + 4 + 3 + 7 + 19 + 3 + 5 + 3 + 6 + 7 + 8 + 3 + 5 + 3 + 4 + 11 + 8 + 6 + 3 + 7 + 5 + 3 + 5 + 4 + 9 + 7 + 3 + 3 + 3 + 7 + 6 + 5 + 24 + 6 + 4 + 5 + 4 + 9 + 4 + 11 + 8 + 5 + 4 + 4 + 8 + 3 + 10 + 7 + 4 + 4 + 10 + 13 + 12 + 9 + 8 + 6 + 5 + 11 + 3 + 20 + 7 + 5 + 7 + 8 + 4 + 5 + 3 + 3 + 3 + 4 + 5 + 7 + 13 + 9 + 4 + 4 + 3 + 17 + 5 + 5 + 5 + 23 + 8 + 7 + 4 + 16 + 3 + 3 + 4 + 10 + 0 + 6 + 5 + 3 + 11 + 18 + 3 + 5 + 3 + 3 + 11 + 3 + 3 + 6 + 7 + 12 + 11 + 8 + 3 + 6 + 4 + 6 + 3 + 5 + 7 + 3 + 3 + 3 + 4 + 8 + 6 + 8 + 4 + 6 + 4 + 7 + 5 + 5 + 20 + 6 + 5 + 9 + 7 + 4 + 3 + 10 + 5 + 18 + 4 + 10 + 3 + 4 + 3 + 4 + 4 + 5 + 10 + 12 + 6 + 3 + 7 + 10 + 4 + 5 + 8 + 6 + 6 + 5 + 4 + 3 + 4 + 5 + 8 + 13 + 5 + 4 + 7 + 8 + 4 + 3 + 0 + 4 + 4 + 5 + 3 + 4 + 10 + 3 + 3 + 3 + 13 + 4 + 4 + 4 + 7 + 15 + 8 + 13 + 9 + 6 + 8 + 4 + 5 + 5 + 3 + 4 + 8 + 5 + 7 + 6 + 13 + 8 + 3 + 8 + 4 + 8 + 6 + 5 + 3 + 4 + 3 + 10 + 4 + 6 + 3 + 4 + 10 + 14 + 12 + 3 + 11 + 6 + 3 + 8 + 30 + 14 + 11 + 3 + 8 + 5 + 9 + 13 + 7 + 6 + 28 + 22 + 5 + 13 + 5 + 9 + 14 + 40 + 17 + 11 + 4 + 6 + 3 + 9 + 4 + 6 + 8 + 11 + 3 + 7 + 5 + 12 + 6 + 18 + 4 + 3 + 4 + 4 + 3 + 5 + 3 + 7 + 3 + 3 + 8 + 8 + 9 + 3 + 8 + 9 + 5 + 14 + 5 + 5 + 10 + 7 + 21 + 7 + 5 + 5 + 5 + 4 + 5 + 9 + 9 + 4 + 8 + 10 + 7 + 4 + 20 + 12 + 6 + 7 + 17 + 12 + 5 + 4 + 3 + 5 + 3 + 5 + 8 + 6 + 4 + 3 + 6 + 12 + 19 + 19 + 4 + 6 + 7 + 5 + 5 + 10 + 10 + 5 + 8 + 3 + 15 + 6 + 16 + 3 + 3 + 13 + 4 + 3 + 6 + 11 + 3 + 6 + 18 + 12 + 3 + 15 + 3 + 3 + 7 + 3 + 6 + 3 + 9 + 18 + 4 + 3 + 9 + 4 + 9 + 4 + 5 + 9 + 5 + 3 + 7 + 8 + 4 + 3 + 13 + 7 + 7 + 12 + 22 + 4 + 22 + 5 + 7 + 3 + 3 + 7 + 9 + 6 + 6 + 6 + 5 + 5 + 4 + 6 + 15 + 5 + 4 + 3 + 8 + 9 + 13 + 7 + 15 + 3 + 3 + 11 + 16 + 17 + 4 + 9 + 5 + 3 + 4 + 6 + 5 + 3 + 10 + 8 + 3 + 13 + 6 + 9 + 4 + 3 + 10 + 9 + 3 + 3 + 8 + 9 + 3 + 13 + 8 + 3 + 27 + 28 + 3 + 3 + 11 + 4 + 3 + 21 + 20 + 5 + 3 + 9 + 6 + 5 + 9 = 26901 vertices

All time in microseconds
Time: Total : 261518
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 0
Shape file2: ../datasets/continents.csv QQID: 1524
PP Polygon size 100612 QQ Polygon size 5574
PP Count 100612 QQ Count 5574
MBR_P [-180, -85.2219, 180, 90
MBR_Q [-81.3551, -53.8864, 0, 12.4639
CMBR [-81.3551, -53.8864, 0, 12.4639

P overlap count with CMBR 2337 Q overlap count with CMBR 5574 

Non-degen count P 2 *****--- Q 2
Intersection count P 2 *****--- Q 2

gpuCountIntersections kernel exe time(ms) 59.813122
prefixsum kernels exe time(ms) 1.814848
gpuNeighborMap kernel exe time(ms) 6.232320
gpuCalculateIntersections kernel exe time(ms) 0.431808
gpuSortPolyQ kernel exe time(ms) 0.005920
gpuCalculateInitLabel kernel exe time(ms) 0.010240

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 0 vertices removed

R has 1 component with 5577 vertices

All time in microseconds
Time: Total : 226842
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 0
Shape file2: ../datasets/continents.csv QQID: 54
PP Polygon size 100612 QQ Polygon size 4987
PP Count 100612 QQ Count 4987
MBR_P [-180, -85.2219, 180, 90
MBR_Q [-17.5328, -34.822, 51.4113, 37.3404
CMBR [-17.5328, -34.822, 51.4113, 37.3404

P overlap count with CMBR 18525 Q overlap count with CMBR 4987 

Non-degen count P 3832 *****--- Q 3832
Intersection count P 3832 *****--- Q 3832

gpuCountIntersections kernel exe time(ms) 65.672287
prefixsum kernels exe time(ms) 1.804320
gpuNeighborMap kernel exe time(ms) 7.976864
gpuCalculateIntersections kernel exe time(ms) 2.818048
gpuSortPolyQ kernel exe time(ms) 0.059936
gpuCalculateInitLabel kernel exe time(ms) 0.022528

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 25 vertices removed

R has 1904 components with 0 + 6 + 7 + 16 + 16 + 12 + 6 + 3 + 36 + 11 + 3 + 3 + 5 + 9 + 3 + 3 + 3 + 4 + 6 + 11 + 13 + 10 + 3 + 8 + 6 + 6 + 5 + 3 + 13 + 4 + 7 + 9 + 8 + 4 + 4 + 5 + 6 + 5 + 17 + 6 + 4 + 3 + 3 + 3 + 3 + 4 + 7 + 5 + 3 + 4 + 7 + 3 + 6 + 8 + 11 + 5 + 3 + 4 + 5 + 5 + 3 + 3 + 4 + 4 + 4 + 6 + 5 + 5 + 7 + 3 + 5 + 6 + 3 + 3 + 9 + 4 + 4 + 3 + 5 + 3 + 4 + 4 + 7 + 5 + 7 + 3 + 3 + 3 + 3 + 3 + 3 + 3 + 8 + 3 + 3 + 6 + 5 + 5 + 3 + 3 + 4 + 10 + 3 + 5 + 6 + 8 + 14 + 4 + 3 + 3 + 6 + 6 + 3 + 8 + 6 + 5 + 3 + 5 + 4 + 5 + 5 + 3 + 4 + 3 + 3 + 3 + 3 + 3 + 4 + 3 + 5 + 3 + 9 + 5 + 7 + 3 + 3 + 5 + 9 + 5 + 6 + 19 + 3 + 17 + 4 + 6 + 4 + 4 + 5 + 4 + 5 + 16 + 3 + 3 + 7 + 3 + 4 + 13 + 3 + 8 + 4 + 3 + 4 + 7 + 4 + 3 + 4 + 4 + 4 + 4 + 10 + 4 + 5 + 5 + 4 + 3 + 3 + 3 + 3 + 8 + 7 + 4 + 3 + 4 + 3 + 5 + 3 + 4 + 3 + 3 + 3 + 3 + 4 + 3 + 3 + 5 + 4 + 3 + 6 + 5 + 3 + 4 + 5 + 3 + 3 + 4 + 3 + 3 + 3 + 4 + 4 + 7 + 5 + 4 + 6 + 3 + 5 + 3 + 4 + 8 + 3 + 3 + 4 + 3 + 4 + 3 + 4 + 3 + 3 + 5 + 4 + 5 + 4 + 3 + 4 + 5 + 4 + 3 + 3 + 3 + 8 + 11 + 3 + 0 + 3 + 6 + 6 + 3 + 6 + 9 + 6 + 3 + 3 + 5 + 5 + 3 + 6 + 4 + 4 + 3 + 4 + 3 + 5 + 3 + 4 + 7 + 4 + 5 + 6 + 6 + 5 + 3 + 4 + 6 + 3 + 3 + 8 + 6 + 5 + 4 + 3 + 4 + 6 + 5 + 3 + 4 + 4 + 3 + 5 + 5 + 3 + 4 + 3 + 3 + 5 + 3 + 3 + 3 + 3 + 3 + 3 + 3 + 3 + 3 + 5 + 10 + 3 + 3 + 3 + 4 + 3 + 4 + 7 + 7 + 10 + 5 + 7 + 6 + 3 + 5 + 5 + 5 + 11 + 8 + 3 + 3 + 5 + 3 + 3 + 3 + 5 + 4 + 3 + 9 + 3 + 3 + 3 + 4 + 7 + 3 + 6 + 5 + 11 + 9 + 5 + 5 + 9 + 5 + 5 + 13 + 4 + 7 + 3 + 3 + 5 + 5 + 10 + 4 + 9 + 12 + 3 + 7 + 6 + 4 + 6 + 4 + 5 + 5 + 3 + 5 + 5 + 5 + 3 + 3 + 6 + 10 + 10 + 3 + 6 + 3 + 10 + 5 + 7 + 3 + 3 + 5 + 5 + 4 + 6 + 3 + 6 + 9 + 4 + 3 + 5 + 6 + 4 + 3 + 4 + 3 + 5 + 6 + 5 + 6 + 11 + 7 + 4 + 11 + 32 + 7 + 11 + 4 + 8 + 4 + 3 + 4 + 3 + 4 + 6 + 18 + 7 + 8 + 5 + 11 + 30 + 4 + 27 + 3 + 11 + 6 + 4 + 4 + 20 + 7 + 3 + 3 + 3 + 7 + 4 + 5 + 10 + 8 + 8 + 14 + 3 + 3 + 4 + 5 + 4 + 5 + 7 + 5 + 3 + 3 + 4 + 7 + 4 + 5 + 5 + 8 + 5 + 14 + 11 + 6 + 3 + 3 + 5 + 21 + 41 + 23 + 6 + 3 + 3 + 7 + 4 + 6 + 4 + 3 + 3 + 4 + 3 + 4 + 3 + 3 + 3 + 3 + 6 + 4 + 3 + 7 + 12 + 3 + 11 + 9 + 11 + 8 + 23 + 6 + 12 + 9 + 8 + 45 + 7 + 6 + 3 + 9 + 3 + 5 + 4 + 3 + 3 + 29 + 12 + 6 + 5 + 8 + 22 + 6 + 4 + 5 + 15 + 3 + 3 + 6 + 3 + 7 + 9 + 10 + 12 + 3 + 22 + 8 + 3 + 9 + 79 + 4 + 11 + 10 + 17 + 0 + 9 + 17 + 17 + 20 + 7 + 4 + 7 + 3 + 25 + 5 + 3 + 3 + 10 + 7 + 12 + 4 + 6 + 7 + 4 + 3 + 11 + 12 + 18 + 5 + 9 + 18 + 18 + 3 + 4 + 8 + 61 + 5 + 25 + 15 + 10 + 3 + 5 + 19 + 9 + 33 + 3 + 9 + 34 + 31 + 20 + 19 + 5 + 3 + 3 + 6 + 6 + 10 + 10 + 6 + 8 + 3 + 12 + 48 + 28 + 3 + 5 + 3 + 13 + 15 + 23 + 5 + 11 + 3 + 36 + 13 + 9 + 3 + 11 + 7 + 5 + 4 + 4 + 8 + 3 + 7 + 4 + 4 + 10 + 16 + 4 + 8 + 8 + 9 + 80 + 5 + 7 + 3 + 3 + 3 + 3 + 21 + 8 + 10 + 17 + 6 + 3 + 4 + 12 + 26 + 8 + 10 + 7 + 3 + 19 + 5 + 9 + 5 + 6 + 8 + 6 + 9 + 3 + 3 + 4 + 14 + 5 + 7 + 24 + 7 + 3 + 3 + 8 + 22 + 10 + 3 + 16 + 5 + 3 + 3 + 3 + 3 + 13 + 3 + 3 + 3 + 3 + 6 + 7 + 4 + 3 + 5 + 14 + 7 + 7 + 4 + 5 + 14 + 17 + 6 + 4 + 19 + 7 + 5 + 47 + 6 + 3 + 6 + 6 + 3 + 3 + 4 + 4 + 6 + 3 + 16 + 5 + 3 + 4 + 3 + 5 + 8 + 7 + 3 + 5 + 7 + 4 + 3 + 6 + 5 + 19 + 9 + 16 + 12 + 9 + 12 + 8 + 7 + 7 + 6 + 4 + 4 + 8 + 6 + 9 + 13 + 6 + 5 + 8 + 8 + 7 + 5 + 9 + 5 + 6 + 7 + 10 + 26 + 4 + 3 + 27 + 5 + 4 + 5 + 4 + 8 + 14 + 26 + 6 + 16 + 3 + 10 + 3 + 5 + 7 + 4 + 13 + 11 + 4 + 6 + 3 + 11 + 10 + 6 + 3 + 6 + 8 + 48 + 9 + 6 + 5 + 4 + 4 + 3 + 25 + 3 + 7 + 3 + 6 + 17 + 4 + 3 + 7 + 0 + 4 + 6 + 3 + 4 + 16 + 14 + 78 + 5 + 11 + 10 + 9 + 7 + 3 + 3 + 15 + 21 + 4 + 3 + 8 + 3 + 3 + 7 + 8 + 3 + 3 + 5 + 44 + 5 + 3 + 14 + 4 + 4 + 4 + 5 + 3 + 3 + 4 + 5 + 5 + 4 + 7 + 3 + 11 + 7 + 7 + 3 + 8 + 8 + 28 + 4 + 11 + 4 + 20 + 7 + 10 + 13 + 5 + 3 + 8 + 4 + 3 + 4 + 56 + 19 + 13 + 15 + 12 + 3 + 3 + 3 + 19 + 10 + 3 + 3 + 8 + 5 + 3 + 3 + 10 + 4 + 7 + 10 + 13 + 25 + 8 + 3 + 44 + 12 + 4 + 17 + 9 + 3 + 5 + 7 + 7 + 5 + 5 + 5 + 11 + 3 + 3 + 6 + 3 + 3 + 6 + 3 + 14 + 5 + 7 + 3 + 3 + 7 + 4 + 11 + 5 + 4 + 3 + 11 + 7 + 3 + 8 + 5 + 11 + 18 + 7 + 3 + 3 + 10 + 4 + 22 + 18 + 6 + 8 + 3 + 7 + 12 + 14 + 4 + 3 + 11 + 6 + 4 + 4 + 5 + 15 + 3 + 8 + 16 + 7 + 4 + 12 + 16 + 3 + 10 + 7 + 24 + 20 + 8 + 5 + 4 + 12 + 27 + 6 + 5 + 3 + 6 + 8 + 7 + 3 + 4 + 3 + 3 + 5 + 9 + 3 + 3 + 7 + 6 + 3 + 3 + 3 + 3 + 3 + 7 + 9 + 3 + 3 + 6 + 3 + 9 + 5 + 3 + 11 + 21 + 8 + 21 + 3 + 14 + 5 + 3 + 4 + 6 + 12 + 5 + 3 + 6 + 3 + 5 + 3 + 4 + 3 + 6 + 7 + 9 + 20 + 4 + 37 + 4 + 6 + 4 + 7 + 5 + 7 + 6 + 10 + 3 + 6 + 4 + 4 + 11 + 41 + 11 + 4 + 7 + 3 + 11 + 4 + 5 + 3 + 3 + 21 + 3 + 4 + 4 + 3 + 7 + 3 + 7 + 6 + 3 + 3 + 6 + 3 + 3 + 5 + 4 + 5 + 5 + 10 + 4 + 10 + 3 + 5 + 7 + 8 + 6 + 3 + 4 + 16 + 3 + 16 + 3 + 6 + 3 + 3 + 4 + 3 + 3 + 3 + 7 + 4 + 13 + 4 + 6 + 6 + 3 + 5 + 3 + 13 + 4 + 5 + 6 + 3 + 6 + 4 + 23 + 4 + 3 + 3 + 3 + 3 + 5 + 3 + 4 + 6 + 4 + 5 + 4 + 4 + 8 + 24 + 4 + 12 + 9 + 5 + 6 + 8 + 11 + 3 + 6 + 3 + 6 + 10 + 13 + 3 + 7 + 5 + 5 + 4 + 5 + 0 + 4 + 5 + 3 + 3 + 5 + 3 + 5 + 3 + 4 + 9 + 4 + 6 + 3 + 9 + 12 + 5 + 4 + 5 + 3 + 5 + 4 + 6 + 9 + 7 + 8 + 5 + 4 + 15 + 3 + 4 + 5 + 3 + 7 + 6 + 3 + 5 + 6 + 8 + 8 + 7 + 3 + 8 + 3 + 5 + 4 + 6 + 3 + 11 + 5 + 10 + 7 + 7 + 4 + 7 + 4 + 5 + 10 + 3 + 5 + 6 + 3 + 3 + 5 + 3 + 15 + 14 + 5 + 7 + 6 + 4 + 12 + 16 + 3 + 3 + 3 + 14 + 5 + 3 + 9 + 19 + 5 + 6 + 3 + 6 + 8 + 3 + 5 + 3 + 3 + 17 + 12 + 8 + 3 + 7 + 4 + 3 + 3 + 7 + 7 + 7 + 3 + 6 + 13 + 6 + 4 + 14 + 3 + 6 + 3 + 8 + 3 + 10 + 5 + 3 + 5 + 4 + 8 + 6 + 5 + 7 + 4 + 6 + 9 + 14 + 12 + 16 + 5 + 4 + 8 + 16 + 13 + 5 + 25 + 6 + 12 + 3 + 5 + 3 + 3 + 5 + 8 + 5 + 8 + 7 + 3 + 9 + 5 + 13 + 4 + 3 + 3 + 3 + 3 + 3 + 4 + 4 + 9 + 13 + 6 + 3 + 28 + 12 + 7 + 14 + 20 + 7 + 3 + 4 + 11 + 5 + 8 + 4 + 30 + 16 + 5 + 14 + 12 + 7 + 3 + 3 + 3 + 3 + 5 + 6 + 23 + 17 + 19 + 5 + 9 + 6 + 4 + 12 + 7 + 3 + 5 + 7 + 5 + 11 + 11 + 33 + 4 + 9 + 4 + 5 + 5 + 28 + 14 + 6 + 9 + 4 + 5 + 4 + 5 + 3 + 3 + 4 + 3 + 6 + 3 + 4 + 5 + 4 + 9 + 5 + 4 + 3 + 6 + 4 + 7 + 4 + 12 + 8 + 6 + 4 + 9 + 3 + 3 + 3 + 3 + 6 + 3 + 8 + 14 + 19 + 3 + 6 + 13 + 5 + 3 + 3 + 11 + 6 + 6 + 13 + 9 + 15 + 3 + 5 + 13 + 3 + 5 + 23 + 6 + 4 + 6 + 6 + 4 + 16 + 5 + 3 + 14 + 10 + 4 + 20 + 13 + 16 + 3 + 8 + 3 + 4 + 15 + 8 + 5 + 3 + 7 + 5 + 7 + 3 + 4 + 16 + 7 + 5 + 28 + 6 + 9 + 14 + 7 + 4 + 10 + 3 + 35 + 5 + 8 + 7 + 12 + 16 + 22 + 24 + 3 + 4 + 3 + 4 + 10 + 6 + 9 + 21 + 6 + 25 + 5 + 17 + 5 + 5 + 4 + 6 + 7 + 3 + 5 + 9 + 5 + 6 + 3 + 8 + 3 + 3 + 3 + 3 + 4 + 3 + 6 + 4 + 6 + 8 + 3 + 9 + 5 + 15 + 4 + 4 + 9 + 5 + 3 + 3 + 5 + 11 + 5 + 4 + 11 + 14 + 5 + 8 + 6 + 12 + 9 + 13 + 5 + 6 + 8 + 6 + 5 + 3 + 12 + 4 + 5 + 7 + 3 + 4 + 4 + 4 + 23 + 4 + 7 + 6 + 34 + 6 + 3 + 26 + 3 + 6 + 3 + 14 + 3 + 5 + 15 + 17 + 8 + 3 + 5 + 3 + 7 + 10 + 9 + 4 + 7 + 4 + 5 + 8 + 4 + 4 + 4 + 6 + 5 + 3 + 6 + 6 + 7 + 3 + 24 + 9 + 8 + 8 + 16 + 3 + 7 + 3 + 22 + 4 + 6 + 4 + 32 + 5 + 16 + 31 + 6 + 15 + 6 + 7 + 7 + 4 + 5 + 11 + 8 + 8 + 3 + 3 + 8 + 9 + 4 + 6 + 4 + 5 + 9 + 8 + 14 + 13 + 6 + 23 + 10 + 15 + 11 + 21 + 29 + 7 + 6 + 14 + 8 + 6 + 4 + 7 + 7 + 3 + 15 + 6 + 8 + 3 + 9 + 3 + 4 + 5 + 4 + 8 + 15 + 4 + 5 + 4 + 10 + 3 + 5 + 15 + 4 + 4 + 3 + 7 + 6 + 7 + 10 + 6 + 3 + 9 + 13 + 4 + 5 + 3 + 6 + 3 + 4 + 9 + 4 + 3 + 5 + 4 + 15 + 26 + 5 + 3 + 8 + 3 + 5 + 3 + 4 + 5 + 4 + 4 + 5 + 4 + 7 + 3 + 3 + 3 + 4 + 3 + 4 + 4 + 3 + 4 + 11 + 4 + 4 + 4 + 3 + 3 + 4 + 4 + 3 + 9 + 4 + 5 + 4 + 3 + 30 + 3 + 3 + 3 + 6 + 3 + 5 + 3 + 8 + 3 + 3 + 3 + 9 + 12 + 3 + 7 + 4 + 3 + 6 + 18 + 3 + 3 + 4 + 7 + 3 + 4 + 6 + 6 + 3 + 8 + 8 + 4 + 8 + 3 + 8 + 5 + 3 + 3 + 3 + 5 + 3 + 5 + 10 + 4 + 3 + 8 + 5 + 10 + 6 + 9 + 14 + 9 + 8 + 5 + 17 + 7 + 6 + 3 + 3 + 5 + 10 + 3 + 4 + 4 + 5 + 5 + 6 + 7 + 7 + 5 + 6 + 3 + 4 + 4 + 4 + 4 + 5 + 5 + 4 + 3 + 3 + 5 + 5 + 10 + 12 + 3 + 4 + 9 + 5 + 7 + 3 + 3 + 7 + 6 + 8 + 6 + 4 + 7 + 6 + 4 + 3 + 5 + 7 + 7 + 3 + 4 + 3 + 9 + 4 + 4 + 10 + 8 + 3 + 7 + 5 + 7 + 4 + 5 + 5 + 3 + 3 + 6 + 4 + 5 + 7 + 4 + 11 + 3 + 3 + 16 + 6 + 8 + 4 + 3 + 5 + 4 + 18 + 3 + 12 + 6 + 3 + 5 + 19 + 4 + 0 + 7 + 6 + 3 + 9 + 3 + 3 + 4 + 3 + 7 + 3 + 4 + 6 + 7 + 4 + 11 + 3 + 3 + 4 + 11 + 5 + 4 + 6 + 5 + 5 + 10 + 3 + 3 + 3 + 4 + 3 + 9 + 6 + 3 + 5 + 3 + 4 + 3 + 7 + 15 + 4 + 11 + 8 + 6 + 5 + 4 + 5 + 4 + 23 + 5 + 4 + 3 + 5 + 7 + 22 + 5 + 39 + 3 + 14 + 16 + 3 + 9 + 9 + 7 + 5 + 10 + 9 + 20 + 11 + 3 + 5 + 3 + 3 + 11 + 3 + 5 + 6 + 5 + 5 + 4 + 10 + 3 + 8 = 13595 vertices

All time in microseconds
Time: Total : 221138
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 36
Shape file2: ../datasets/continents.csv QQID: 1048
PP Polygon size 66475 QQ Polygon size 15653
PP Count 66475 QQ Count 15653
MBR_P [-168.137, -53.886, -34.7936, 72.0021
MBR_Q [-168.132, 7.20611, 0, 71.9969
CMBR [-168.132, 7.20611, -34.7936, 71.9969

P overlap count with CMBR 51442 Q overlap count with CMBR 15653 

Non-degen count P 11412 *****--- Q 11412
Intersection count P 11412 *****--- Q 11412

gpuCountIntersections kernel exe time(ms) 130.346466
prefixsum kernels exe time(ms) 1.416928
gpuNeighborMap kernel exe time(ms) 9.154080
gpuCalculateIntersections kernel exe time(ms) 14.794432
gpuSortPolyQ kernel exe time(ms) 0.059072
gpuCalculateInitLabel kernel exe time(ms) 0.032192

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 23 vertices removed

R has 43 components with 44836 + 3 + 4 + 4 + 3 + 4 + 3 + 3 + 7 + 7 + 4 + 4 + 4 + 4 + 3 + 4 + 3 + 8 + 3 + 3 + 4 + 5 + 5 + 4 + 3 + 4 + 3 + 3 + 5 + 6 + 3 + 4 + 3 + 11 + 6 + 6 + 5 + 6 + 5 + 5 + 3 + 4 + 3 = 45020 vertices

All time in microseconds
Time: Total : 300693
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 36
Shape file2: ../datasets/continents.csv QQID: 1524
PP Polygon size 66475 QQ Polygon size 5574
PP Count 66475 QQ Count 5574
MBR_P [-168.137, -53.886, -34.7936, 72.0021
MBR_Q [-81.3551, -53.8864, 0, 12.4639
CMBR [-81.3551, -53.886, -34.7936, 12.4639

P overlap count with CMBR 18345 Q overlap count with CMBR 5574 

Non-degen count P 4440 *****--- Q 4440
Intersection count P 4440 *****--- Q 4440

gpuCountIntersections kernel exe time(ms) 45.530785
prefixsum kernels exe time(ms) 1.214368
gpuNeighborMap kernel exe time(ms) 6.002688
gpuCalculateIntersections kernel exe time(ms) 3.359040
gpuSortPolyQ kernel exe time(ms) 0.054720
gpuCalculateInitLabel kernel exe time(ms) 0.020512

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 9 vertices removed

R has 6 components with 16346 + 3 + 3 + 3 + 4 + 5 = 16364 vertices

All time in microseconds
Time: Total : 191608
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2742**************************
Shape file2: ../datasets/continents.csv QQID: 1048
PP and QQ swapped since QQ> PP
New PP Polygon size 15653 New QQ Polygon size 15547
PP Count 15653 QQ Count 15547
MBR_P [-168.132, 7.20611, 0, 71.9969
MBR_Q [-73.0572, 59.9879, -11.3768, 83.6341
CMBR [-73.0572, 59.9879, -11.3768, 71.9969

P overlap count with CMBR 179 Q overlap count with CMBR 10170 

Non-degen count P 18 *****--- Q 18
Intersection count P 18 *****--- Q 18

gpuCountIntersections kernel exe time(ms) 18.158209
prefixsum kernels exe time(ms) 0.538528
gpuNeighborMap kernel exe time(ms) 0.984928
gpuCalculateIntersections kernel exe time(ms) 1.094592
gpuSortPolyQ kernel exe time(ms) 0.016384
gpuCalculateInitLabel kernel exe time(ms) 0.014240

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 1 vertices removed

R has 5 components with 7 + 10 + 10 + 10 + 6 = 43 vertices

All time in microseconds
Time: Total : 151605
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2742
Shape file2: ../datasets/continents.csv QQID: 1202
PP Polygon size 15547 QQ Polygon size 12750
PP Count 15547 QQ Count 12750
MBR_P [-73.0572, 59.9879, -11.3768, 83.6341
MBR_Q [-73.0536, 59.9814, 0, 83.6236
CMBR [-73.0536, 59.9879, -11.3768, 83.6236

P overlap count with CMBR 15547 Q overlap count with CMBR 12749 

Non-degen count P 4448 *****--- Q 4448
Intersection count P 4448 *****--- Q 4448

gpuCountIntersections kernel exe time(ms) 28.858624
prefixsum kernels exe time(ms) 0.480064
gpuNeighborMap kernel exe time(ms) 3.147936
gpuCalculateIntersections kernel exe time(ms) 2.874464
gpuSortPolyQ kernel exe time(ms) 0.045120
gpuCalculateInitLabel kernel exe time(ms) 0.022208

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 142 vertices removed

R has 14 components with 19250 + 4 + 4 + 4 + 34 + 20 + 7 + 17 + 5 + 3 + 4 + 6 + 9 + 7 = 19374 vertices

All time in microseconds
Time: Total : 167320
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2742*******************
Shape file2: ../datasets/continents.csv QQID: 1081
PP Polygon size 15547 QQ Polygon size 4562
PP Count 15547 QQ Count 4562
MBR_P [-73.0572, 59.9879, -11.3768, 83.6341
MBR_Q [-90.0486, 61.8583, 0, 73.85
CMBR [-73.0572, 61.8583, -11.3768, 73.85

P overlap count with CMBR 9294 Q overlap count with CMBR 2398 

Non-degen count P 14 *****--- Q 14
Intersection count P 14 *****--- Q 14

gpuCountIntersections kernel exe time(ms) 9.125248
prefixsum kernels exe time(ms) 0.335872
gpuNeighborMap kernel exe time(ms) 0.989696
gpuCalculateIntersections kernel exe time(ms) 0.322336
gpuSortPolyQ kernel exe time(ms) 0.040960
gpuCalculateInitLabel kernel exe time(ms) 0.008416

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 0 vertices removed

R has 3 components with 10 + 5 + 5 = 20 vertices

All time in microseconds
Time: Total : 139149
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2742********************************+++++++++++++********
Shape file2: ../datasets/continents.csv QQID: 1193
PP Polygon size 15547 QQ Polygon size 4028
PP Count 15547 QQ Count 4028
MBR_P [-73.0572, 59.9879, -11.3768, 83.6341
MBR_Q [-91.953, 76.128, 0, 83.1139
CMBR [-73.0572, 76.128, -11.3768, 83.1139

P overlap count with CMBR 3224 Q overlap count with CMBR 780 

Non-degen count P 32 *****--- Q 32
Intersection count P 32 *****--- Q 32

gpuCountIntersections kernel exe time(ms) 5.746112
prefixsum kernels exe time(ms) 0.419520
gpuNeighborMap kernel exe time(ms) 1.003360
gpuCalculateIntersections kernel exe time(ms) 0.290656
gpuSortPolyQ kernel exe time(ms) 0.073792
gpuCalculateInitLabel kernel exe time(ms) 0.009952

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 0 vertices removed

R has 8 components with 4 + 4 + 4 + 4 + 4 + 5 + 4 + 4 = 33 vertices

All time in microseconds
Time: Total : 157997
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2741**************************
Shape file2: ../datasets/continents.csv QQID: 1048
PP and QQ swapped since QQ> PP
New PP Polygon size 15653 New QQ Polygon size 10887
PP Count 15653 QQ Count 10887
MBR_P [-168.132, 7.20611, 0, 71.9969
MBR_Q [-90.1186, 61.8581, -61.2675, 73.8553
CMBR [-90.1186, 61.8581, -61.2675, 71.9969

P overlap count with CMBR 928 Q overlap count with CMBR 9694 

Non-degen count P 32 *****--- Q 32
Intersection count P 32 *****--- Q 32

gpuCountIntersections kernel exe time(ms) 16.931423
prefixsum kernels exe time(ms) 0.441600
gpuNeighborMap kernel exe time(ms) 0.992256
gpuCalculateIntersections kernel exe time(ms) 0.868352
gpuSortPolyQ kernel exe time(ms) 0.005184
gpuCalculateInitLabel kernel exe time(ms) 0.024192

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 0 vertices removed

R has 8 components with 20 + 18 + 13 + 18 + 18 + 38 + 27 + 4 = 156 vertices

All time in microseconds
Time: Total : 151373
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2741
Shape file2: ../datasets/continents.csv QQID: 1081
PP Polygon size 10887 QQ Polygon size 4562
PP Count 10887 QQ Count 4562
MBR_P [-90.1186, 61.8581, -61.2675, 73.8553
MBR_Q [-90.0486, 61.8583, 0, 73.85
CMBR [-90.0486, 61.8583, -61.2675, 73.85

P overlap count with CMBR 10879 Q overlap count with CMBR 4562 

Non-degen count P 3138 *****--- Q 3138
Intersection count P 3138 *****--- Q 3138

gpuCountIntersections kernel exe time(ms) 7.584576
prefixsum kernels exe time(ms) 0.258816
gpuNeighborMap kernel exe time(ms) 1.771616
gpuCalculateIntersections kernel exe time(ms) 1.361920
gpuSortPolyQ kernel exe time(ms) 0.058496
gpuCalculateInitLabel kernel exe time(ms) 0.012640

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 9 vertices removed

R has 5 components with 11519 + 3 + 3 + 4 + 5 = 11534 vertices

All time in microseconds
Time: Total : 148847
____________________________________________________________________
Shape file1: ../datasets/ne_10m_ocean.csv PPID: 2854*******************
Shape file2: ../datasets/continents.csv QQID: 1193
PP and QQ swapped since QQ> PP
New PP Polygon size 4028 New QQ Polygon size 3838
PP Count 4028 QQ Count 3838
MBR_P [-91.953, 76.128, 0, 83.1139
MBR_Q [-91.9566, 76.1291, -61.0855, 83.1165
CMBR [-91.953, 76.1291, -61.0855, 83.1139

P overlap count with CMBR 4027 Q overlap count with CMBR 3837 

Non-degen count P 1524 *****--- Q 1524
Intersection count P 1524 *****--- Q 1524

gpuCountIntersections kernel exe time(ms) 2.585248
prefixsum kernels exe time(ms) 0.141280
gpuNeighborMap kernel exe time(ms) 0.693376
gpuCalculateIntersections kernel exe time(ms) 0.803040
gpuSortPolyQ kernel exe time(ms) 0.022336
gpuCalculateInitLabel kernel exe time(ms) 0.013088

Copying completed
Labelling intersections...

... 0 delayed crossings and 0 delayed bouncings
... 0 interior and 0 identical components added to result
... 0 bouncing vertex pairs split

Creating result...

Post-processing...

... 77 vertices removed

R has 1 component with 5369 vertices

All time in microseconds
Time: Total : 127744
____________________________________________________________________

==========================================================================================


==========================================================================================




____________________________________________________________________
____________________________________________________________________
____________________________________________________________________
____________________________________________________________________
____________________________________________________________________
____________________________________________________________________




-----------------------******************* Parks Dataset **********************-------------------------
11548832 polygons found
P Polygons
127953 684356
123848 10613514
101164 1106503
100453 876440
91736 1065174
91046 1243058
87763 943452
84530 844165
81897 169840
79686 321571
79602 140315
77367 173408
75543 1789086
70179 1312358
69283 1502030
63802 729699
60397 34579
56751 679995
54992 34622
52175 1512242


------------------------********************* lakes Dataset ***********************---------------------
8330134 polygons found
P Polygons
206429 851348
196726 132372
174386 7890016
102721 174690
89713 7997546
78362 332224
76946 8043478
73907 492879
71739 7907419
70430 148346
68511 2418
66533 41116
66114 360489
62213 593561
62138 58746
61070 111678
59708 78019
55111 2422
54413 695318
51664 803987