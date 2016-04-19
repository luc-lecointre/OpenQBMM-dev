#!/bin/bash
# Reads value from a specified cell and writes it with its variable name

#echo "//Initial values" > 0/initC
time="../combustion/0.23.flow"
dir=$(pwd)/0

init()
{
    for filename in *
    do 
        if [ $filename != "phi" ] && [ $filename != "T" ] && [ $filename != "U" ] && [ $filename != "O2" ] && [ $filename != "N2" ] && [ $filename != "alphat" ] && [ $filename != "p" ] && [ $filename != "p_rgh" ]; then
            if [ -f "$filename" ]
            then
                message=`head -122 $filename | echo "        value           uniform "$(tail -1)";"`
                cd $dir
                cp O2 $filename
                sed -i "27s/.*/$message/" $filename
                cd ../../combustion/0.23.flow
            fi                    
        fi 
    done 
}

(cd $time && init)
