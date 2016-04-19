#!/bin/bash
# Reads value from a specified cell and writes it with its variable name

echo "//Initial values" > 0.oxidation/initC
time="0.23.flow"
dir=$(pwd)/0.oxidation

init()
{
    for filename in * 
    do 
        if [ -f "$filename" ] 
        then 
            head -122 $filename | echo -e "scalar $"$filename"("$(tail -1)");" >> $dir/initC
        fi 
    done 
}

(cd $time && init)
