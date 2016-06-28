#!/bin/bash
# Reads value from a specified cell and writes it with its variable name

#echo "//Initial values" > 0/initC
time=$(pwd)/../combustion/0.2.pbe
dir=$(pwd)/0

init()
{
    for filename in *
    do 
        if [ -f "$filename" ] && [ "$filename" != "phi" ] && [ "$filename" != "T" ] && [ "$filename" != "U" ] && [ "$filename" != "O2" ] && [ "$filename" != "N2" ] && [ "$filename" != "alphat" ] && [ "$filename" != "p" ] && [ "$filename" != "p_rgh" ]; then
            message=`head -122 $filename | echo "        value           uniform "$(tail -1)";"`
            cd $dir
            cp Ydefault $filename
            sed -i "27s/.*/$message/" $filename
            cd ../../combustion/0.2.pbe
        fi
        cd $dir
        if [ "$filename" = "moment.0.populationBalance" ]
        then
            message='dimensions      [0 -3 0 0 0 0 0];'
            sed -i "18s/.*/$message/" $filename
        elif [ "$filename" = "moment.1.populationBalance" ]
        then
            message='dimensions      [0 0 0 0 0 0 0];'
            sed -i "18s/.*/$message/" $filename
        elif [ "$filename" = "moment.2.populationBalance" ]
        then
            message='dimensions      [0 3 0 0 0 0 0];'
            sed -i "18s/.*/$message/" $filename
        elif [ "$filename" = "moment.3.populationBalance" ]
        then
            message='dimensions      [0 6 0 0 0 0 0];'
            sed -i "18s/.*/$message/" $filename
        elif [ "$filename" = "moment.4.populationBalance" ]
        then
            message='dimensions      [0 9 0 0 0 0 0];'
            sed -i "18s/.*/$message/" $filename
        elif [ "$filename" = "moment.5.populationBalance" ]
        then
            message='dimensions      [0 12 0 0 0 0 0];'
            sed -i "18s/.*/$message/" $filename
        elif [ "$filename" = "moment.6.populationBalance" ]
        then
            message='dimensions      [0 15 0 0 0 0 0];'
            sed -i "18s/.*/$message/" $filename
        fi
        cd ../../combustion/0.2.pbe
    done 
}

cp -f 0.org/* 0/
(cd $time && init)
