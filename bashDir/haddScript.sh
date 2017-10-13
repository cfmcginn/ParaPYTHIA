#!/bin/bash

if [ "$#" -ne 1 ]
then
    echo "Need input"
    exit
fi

count=0
looper=`ls $1`
for i in $looper
do
    count=$(($count + 1))
done

echo $count

divval=$(($count/12))
echo $divval
count=0

stringout=''

mergeval=0

for i in $looper
do
    if [ $(($count%$divval)) -eq 0 ]
    then
	if [ $count -ne 0 ]
	then	
	    echo hadd $1/outFile_MERGED_$mergeval.root $stringout >& mergeLog_$mergeval.log
	    hadd $1/outFile_MERGED_$mergeval.root $stringout >> mergeLog_$mergeval.log &
	    stringout=''
	    
	    mergeval=$(($mergeval + 1))
	fi
    fi

    stringout="$stringout $1/$i"

    count=$(($count + 1))
done

echo "Left: "
echo $stringout

hadd $1/outFile_MERGED_$mergeval.root $stringout >& mergeLog_$mergeval.log &
