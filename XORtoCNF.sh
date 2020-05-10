#!/bin/bash
offset=$(($1+1))
offsetNeg=$((0-$offset))
cat $2 | grep -v 'x ' > F_3_norm.dimacs
while read -r c
do
	clause=$(echo $c)
    i=0
	for l in $clause
	do
        arr[$i]=$l
        arrNeg[$i]=$((0-$l))
        i=$(($i+1))
	done
    len=$i
    if [ $len -eq 4 ]
    then
        echo "${arr[1]} ${arr[2]} 0" >> F_3_norm.dimacs
        echo "${arrNeg[1]} ${arrNeg[2]} 0" >> F_3_norm.dimacs
    fi
    
    if [ $len -eq 5 ]
    then
        echo "${arrNeg[1]} ${arrNeg[2]} ${arr[3]} 0" >> F_3_norm.dimacs
        echo "${arr[1]} ${arrNeg[2]} ${arrNeg[3]} 0" >> F_3_norm.dimacs
        echo "${arrNeg[1]} ${arr[2]} ${arrNeg[3]} 0" >> F_3_norm.dimacs
        echo "${arr[1]} ${arr[2]} ${arr[3]} 0" >> F_3_norm.dimacs
    fi
    
    if [ $len -gt 5 ]
    then
        echo "${arrNeg[1]} ${arrNeg[2]} $offset 0" >> F_3_norm.dimacs
        echo "${arr[1]} ${arrNeg[2]} $offsetNeg 0" >> F_3_norm.dimacs
        echo "${arrNeg[1]} ${arr[2]} $offsetNeg 0" >> F_3_norm.dimacs
        echo "${arr[1]} ${arr[2]} $offset 0" >> F_3_norm.dimacs
        prev=$offsetNeg
        prevNeg=$offset
        offset=$(($offset+1))
        offsetNeg=$((0-$offset))
        i=3
        len=$(($len-2))
        while [ $len -gt 4 ]
        do
            echo "$prevNeg ${arrNeg[$i]} $offset 0" >> F_3_norm.dimacs
            echo "$prev ${arrNeg[$i]} $offsetNeg 0" >> F_3_norm.dimacs
            echo "$prevNeg ${arr[$i]} $offsetNeg 0" >> F_3_norm.dimacs
            echo "$prev ${arr[$i]} $offset 0" >> F_3_norm.dimacs
            prev=$offsetNeg
            prevNeg=$offset
            offset=$(($offset+1))
            offsetNeg=$((0-$offset))
            i=$(($i+1))
            len=$(($len-1))
        done
        j=$(($i+1))
        echo "$prevNeg ${arrNeg[$i]} ${arr[$j]} 0" >> F_3_norm.dimacs
        echo "$prev ${arrNeg[$i]} ${arrNeg[$j]} 0" >> F_3_norm.dimacs
        echo "$prevNeg ${arr[$i]} ${arrNeg[$j]} 0" >> F_3_norm.dimacs
        echo "$prev ${arr[$i]} ${arr[$j]} 0" >> F_3_norm.dimacs
    fi
done <<< "$(grep 'x ' $2)"
echo $offset
