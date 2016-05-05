#!/bin/bash
# Reads value from a specified cell and writes it with its variable name

#echo "//Initial values" > 0/initC
time=$(pwd)/../combustion/0.1.pbe
dir=$(pwd)/0

cp -f 0.org/* 0/

init()
{
    for filename in *
    do 
        if [ $filename != "phi" ] && [ $filename != "T" ] && [ $filename != "U" ] && [ $filename != "O2" ] && [ $filename != "N2" ] && [ $filename != "alphat" ] && [ $filename != "p" ] && [ $filename != "p_rgh" ]; then
            if [ -f "$filename" ]
            then
                message=`head -122 $filename | echo "        value           uniform "$(tail -1)";"`
                cd $dir
                cp Ydefault $filename
                sed -i "27s/.*/$message/" $filename
                cd ../../combustion/0.1.pbe
            fi                    
        fi 
    done 
}


(cd $time && init)
