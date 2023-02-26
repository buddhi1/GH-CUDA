# script run all test cases in real-world dataset

bDataset=("worst-syntheticP_")
oDataset=("worst-syntheticQ_")
base=("lakes_851348" "s" "c" ${bDataset[0]}"500" ${bDataset[0]}"700" ${bDataset[0]}"1000" ${bDataset[0]}"1500" ${bDataset[0]}"2000")
overlay=("lakes-synthetic_1_1" "s-synthetic_15_10" "c-synthetic_100_0" ${oDataset[0]}"500" ${oDataset[0]}"700" ${oDataset[0]}"1000" ${oDataset[0]}"1500" ${oDataset[0]}"2000")

# change following to select the test cases
start=0
end=7
skip=10
echo "Synthetic dataset execution Starting.."

# echo -e '#\tDataset\tCount Intersections\tPrefixsum\tCretae Map Q List\tSave Intersection\tSort Q\tTotal time' 
str="#, Dataset, Calculate Intersections, Lablel Intersections, Trace Results, Total sequential time\n" 
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
        str="${str} $(./polyclip data/synthetic/${base[$i]}.txt  data/synthetic/${overlay[$i]}.txt optimizedFostersAlgorithm/results/results_${base[$i]}_${overlay[$i]}.poly save)\n"
    else
        str="${str} $(./polyclip data/synthetic/${base[$i]}.txt  data/synthetic/${overlay[$i]}.txt optimizedFostersAlgorithm/results/results_${base[$i]}_${overlay[$i]}.poly)\n"
    fi
    # str="${str} $(./polyclip data/synthetic/${base[$i]}.txt  data/synthetic/${overlay[$i]}.txt results/results_${base[$i]}_${overlay[$i]}.poly)\n"
    # ./polyclip data/synthetic/${base[$i]}.txt  data/synthetic/${overlay[$i]}.txt results/results_${base[$i]}_${overlay[$i]}.poly
done
echo -e "$str" > results/synthetic_sequential.csv
echo "Synthetic dataset execution completed!"