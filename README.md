# Efficient PRAM and Practical GPU Algorithms for Large Polygon Clipping with Degenerate Cases

This repository contains source code of our latest work, which is accepted for CCGRID 2023. 

[![DOI](https://zenodo.org/badge/480433755.svg)](https://zenodo.org/badge/latestdoi/480433755)

## Installation

```
$ git clone https://github.com/buddhi1/GHCUDA.git
$ cd GH-CUDA
$ make
```

*make* file reads the cuda path from *$LD_LIBRARY_PATH* variable. If the path is not correctly found, replace *LIBCUDA* variable in the *make* file with the correct path. 
Use the following commands to unzip data.
```
$ cd data/
$ unzip poly_data.zip
$ cd ..
```
Now data folder should have realworld and synthetic folders with data files in them.

**Debug mode:** Following steps can be used to quickly check if everything is set up correctly. case is the test case id of real-world or synthetic dataset. Results are saved in the *results/ folder*.

* Compile GPU code:
```
$ make clean
$ make
```

* Real-world dataset experiments: (Fastest case=8)
```
$ /bin/bash real-world_dataset_run.sh ⟨case⟩
```
* Synthetic dataset experiments: (Fastest case=4)
```
$ /bin/bash synthetic_dataset_run.sh ⟨case⟩
```
* Compile sequential code:
```
$ make cleanfoster
$ make foster
```
* Real-world dataset experiments: (Fastest case=8)
```
$/bin/bash real-world_dataset_foster.sh ⟨case⟩
```
* Synthetic dataset experiments: (Fastest case=4)
```
$ /bin/bash synthetic_dataset_foster.sh ⟨case⟩
```

**Usage:** Use save flag only when the resulting polygon needs to be saved to a file.
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
/bin/bash real-world_dataset_run.sh
```
* Synthetic dataset experiments:
```
$ /bin/bash synthetic_dataset_run.sh
```


* Compile sequential code:
```
$ make cleanfoster
$ make foster
```
* Real-world dataset experiments:
```
$ /bin/bash real-world_dataset_foster.sh
```
* Synthetic dataset experiments:
```
$ /bin/bash synthetic_dataset_foster.sh
```

These executions generate *real-world_tests.csv*, *synthetic_tests.csv*, *real-world_sequential.csv*, and *syn-
thetic_sequential.csv* spreadsheets in **results/** folder with runtime breakdowns.

Copy real-world parallel execution data in *real-world_tests.csv* and real-world sequential data in *real-world_sequential.csv* to Table 1 in **Real-world_dataset** workbook of the template Excel file. Copy synthetic parallel execution data in *synthetic_tests.csv* and synthetic sequential data in *synthetic-sequential.csv* to the Table in **Synthetic_dataset** workbook of the template Excel file. The template Excel file is **graphs.xlsx** in the **results/** folder.

**Save results:** Use *save* flag at the end of the script calls above to save clipped polygon coordinates in files. Foster’s results are saved in *optimizedFostersAlgorithm/results/* folder and GPU results are saved in *results/* folder. Saved results can be used to compare outputs and we do not recommend using *save* flag for performance evaluation.

## Evaluate Filter Performance

Three experiments need to be completed to compare different filter configurations.
1. CMF and CF filters only:
```
$ make clean
$ make cmfcf
```
* Real-world dataset experiments:
```
$ /bin/bash real-world_dataset_run.sh
```
These executions generate *real-world_tests.csv*, and *synthetic_tests.csv* spreadsheets in **results/** folder. Save data in table 3 case ii of the template Excel.

2. LSMF filter only:
```
$ make clean
$ make lsmf
```
* Real-world dataset experiments:
```
$ /bin/bash real-world_dataset_run.sh
```
Save data in table 3 case iii of the template Excel.

3. No filters:
```
$ make clean
$ make nofilters
```
* Rreal-world dataset experiments:
```
$ /bin/bash real-world_dataset_run.sh
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
$ /bin/bash real-world_dataset_run.sh
```
This execution generates *real-world_tests_counts.csv* spreadsheet in **results/** folder with the edge pair counts after each filter for the above-mentioned kernels. Table 2 in the template excel represents this data

## Reference
Reference to cite when you use our work:
```
@article{
  comming soon
}
```
