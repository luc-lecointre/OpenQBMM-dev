#!/bin/bash
# Reads value from a specified cell and writes it with its variable name

echo "//Initial values" > 0/initC
time="../combustion/0.23.flow"
dir=$(pwd)/0

init()
{
    for filename in * 
    do 
        if [ $filename != "phi" ]
        then
            if [ $filename != "T" ]
            then
                if [ -f "$filename" ]
                then
                    head -122 $filename | echo -e "scalar "$filename"("$(tail -1)");" >> $dir/initC
                fi                    
            fi
        fi 
    done 
}

(cd $time && init)
