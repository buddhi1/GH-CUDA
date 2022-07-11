lakes=("851348" "132372" "7890016" "174690" "7997546" "332224" "8043478" "492879" "7907419" "148346" "2418" "41116" "360489" "593561" "58746" "111678" "78019" "2422" "695318" "803987")
parks=("684356" "10613514" "1106503" "876440" "1065174" "1243058" "943452" "844165" "169840" "321571" "140315" "173408" "1789086" "1312358" "1502030" "729699" "34579" "679995" "34622" "1512242")
# lakes=("695318" "803987")
# parks=("34622" "1512242")


for lake in ${lakes[@]}; do
    for park in ${parks[@]}; do
        echo $lake
        echo $park
        ./program ../datasets/lakes_OSM_new_tiger/lakes_$lake.txt  ../datasets/parks_OSM_new_tiger/parks_$park.txt results/Fig-large-R.poly
    done
done

# ./program ../datasets/lakes_OSM_new_tiger/lakes_803987.txt  ../datasets/parks_OSM_new_tiger/parks_10613514.txt results/Fig-large-R.poly