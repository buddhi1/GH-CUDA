# script run all test cases in real-world dataset

bDataset=("ne_10m_ocean_")
oDataset=("ne_10m_land_" "continents_")
base=("s" ${bDataset[0]}"36" ${bDataset[0]}"0" ${bDataset[0]}"0" ${bDataset[0]}"0" ${oDataset[1]}"1048" ${bDataset[0]}"2742" ${bDataset[0]}"2742" ${oDataset[1]}"1048")
overlay=("c" ${oDataset[0]}"1" ${oDataset[0]}"4" ${oDataset[1]}"521" ${oDataset[1]}"1661" ${bDataset[0]}"2742" ${oDataset[1]}"1081" ${oDataset[1]}"1193" ${bDataset[0]}"2741")


# change following to select the test cases
start=0
end=8
skip=10
echo "Real-world dataset execution starting!"

# echo -e '#\tDataset\tCount Intersections\tPrefixsum\tCretae Map Q List\tSave Intersection\tSort Q\tTotal time' 
str="#, Dataset, totoal candidate edge pairs, IC after CMF, IC after LSMF, IS after CF, IS after LSMF, Total parallel time\n" 

for i in $( eval echo {$start..$end}); 
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

    str="${str} $(./program data/realworld/${base[$i]}.txt  data/realworld/${overlay[$i]}.txt results/results_${base[$i]}_${overlay[$i]}.poly)\n"
    # ./program data/realworld/${base[$i]}.txt  data/realworld/${overlay[$i]}.txt results/results_${base[$i]}_${overlay[$i]}.poly
done
echo -e "$str" > results/real-world_tests_counts.csv
echo "Real-world dataset execution completed!"