# Efficient PRAM and Practical GPU Algorithms for Large Polygon Clipping with Degenerate Cases

This repository contains source code of our latest work, which is accepted for CCGRID 2023. 

## Installation

```
$ git clone https://github.com/buddhi1/GHCUDA.git
$ cd GH-CUDA
$ make
```

Use the following commands to unzip data.
```
$ cd data/
$ unzip poly data.zip
$ cd ..
```


Usage: Use save flag only when the resulting polygon needs to be saved to a file.
```
$ ./program ⟨base polygon file⟩ ⟨overlay polygon file⟩ ⟨result polygon path⟩ save
```
Demonstrated experiments have 3 cases. For convenience, we have included bash scripts to test each case. We recommend testing each case at a time.

## Evaluate Performance
This case requires both sequential and parallel algorithms to execute both real-world and synthetic datasets.

* Compile GPU code:
```
$ make clean
$ make
```
* Real-world dataset experiments:
```
/bin/bash real-world dataset run.sh
```
* Synthetic dataset experiments:
```
$ /bin/bash synthetic dataset run.sh
```


* Compile sequential code:
```
$ make cleanfoster
$ make foster
```
7. Real-world dataset experiments:
```
$ /bin/bash real-world dataset foster.sh
```
8. Synthetic dataset experiments:
```
$ /bin/bash synthetic dataset foster.sh
```

These executions generate *real-world_tests.csv*, *synthetic_tests.csv*, *real-world_sequential.csv*, and *syn-
thetic_sequential.csv* spreadsheets in **results/** folder with runtime breakdowns.

Copy real-world parallel execution data in *real-world_tests.csv* and real-world sequential data in *real-world_sequential.csv* to Table 1 in **Real-world_dataset** workbook of the template Excel file. Copy synthetic parallel execution data in *synthetic_tests.csv* and synthetic sequential data in *synthetic-sequential.csv* to the Table in **Synthetic_dataset** workbook of the template Excel file. The template Excel file is **graphs.xlsx** in the **results/** folder.

## Evaluate Filter Performance

Three experiments need to be completed to compare different filter configurations.
1. CMF and CF filters only:
```
$ make clean
$ make cmfcf
```
* Real-world dataset experiments:
```
$ /bin/bash real-world dataset run.sh
```
These executions generate *real-world_tests.csv*, and *synthetic_tests.csv* spreadsheets in **results/** folder. Save data in table 3 case ii of the template Excel.

2. LSMF filter only:
```
$ make clean
$ make lsmf
```
* Real-world dataset experiments:
```
$ /bin/bash real-world dataset run.sh
```
Save data in table 3 case iii of the template Excel.

3. No filters:
```
$ make clean
$ make nofilters
```
* Rreal-world dataset experiments:
```
$ /bin/bash real-world dataset run.sh
```
Save data in table 3 case iv of the template Excel. Case i uses the data from section A where all filters are employed. Use that data in Table 3 case i. Table 4 and the reduction percentage graphs are auto-generated using Table 3 data.

## Evaluate Edge Pair Reduction by Filters
We record the number of edge pairs after each filter which is employed in Intersection Count and Create M ap Q List kernels.
```
$ make clean
$ make count
```
* real-world dataset experiments:
```
$ /bin/bash real-world dataset run.sh
```
This execution generates *real-world_tests_counts.csv* spreadsheet in **results/** folder with the edge pair counts after each filter for the above-mentioned kernels. Table 2 in the template excel represents this data

## Reference
Reference to cite when you use our work:
```
@article{
  comming soon
}
```
