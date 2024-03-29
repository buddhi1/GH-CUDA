# script run all test cases in synthetic dataset

bDataset=("worst-syntheticP_")
oDataset=("worst-syntheticQ_")
base=("lakes_851348" "s" "c" ${bDataset[0]}"500" ${bDataset[0]}"700" ${bDataset[0]}"1000" ${bDataset[0]}"1500" ${bDataset[0]}"2000")
overlay=("lakes-synthetic_1_1" "s-synthetic_15_10" "c-synthetic_100_0" ${oDataset[0]}"500" ${oDataset[0]}"700" ${oDataset[0]}"1000" ${oDataset[0]}"1500" ${oDataset[0]}"2000")

# instructions
echo "Usage: /bin/bash real-world_dataset_run.sh [option] [save]"
echo -e "Run all: option=<>\nDebug mode: option= [1-9]\n"

re='^[0-9]+$'

if ! [[ $1 =~ $re ]] ; then
    # change following to select the test cases
    start=0
    end=7
    skip=10
    echo "Synthetic dataset execution Starting.."

    # echo -e '#\tDataset\tCount Intersections\tPrefixsum\tCretae Map Q List\tSave Intersection\tSort Q\tTotal time' 
    str="#, Dataset, Count Intersections, Prefixsum, Cretae Map Q List, Save Intersection, Sort Q, Initial Lableing, Total parallel time\n" 
    for i in $( eval echo {$start..$end} ); 
    do
        echo "Test case $(($i+1)) running.."
        if [ $i -eq $skip ];
        then
            echo "Test case: $(($i+1)) skipped!"
            echo "----------------------------------------------------"
            echo
            continue
        fi
        str="${str} $(($i+1)), ${base[$i]} ${overlay[$i]}, "
        if [ "$1" == "save" ]; then
            str="${str} $(./program data/synthetic/${base[$i]}.txt  data/synthetic/${overlay[$i]}.txt results/results_${base[$i]}_${overlay[$i]}.poly save)\n"
        else
            str="${str} $(./program data/synthetic/${base[$i]}.txt  data/synthetic/${overlay[$i]}.txt results/results_${base[$i]}_${overlay[$i]}.poly)\n"
        fi

        # str="${str} $(./program data/synthetic/${base[$i]}.txt  data/synthetic/${overlay[$i]}.txt results/results_${base[$i]}_${overlay[$i]}.poly)\n"
        # ./program data/synthetic/${base[$i]}.txt  data/synthetic/${overlay[$i]}.txt results/results_${base[$i]}_${overlay[$i]}.poly
    done
else
    echo "Real-world test case "$1" execution starting!"
    ./program data/synthetic/${base[$(($1-1))]}.txt  data/synthetic/${overlay[$(($1-1))]}.txt results/results_${base[$(($1-1))]}_${overlay[$(($1-1))]}.poly save debug
fi
echo -e "$str" > results/synthetic_tests.csv
echo "Synthetic dataset execution completed!"