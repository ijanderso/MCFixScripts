#!/bin/bash

inputDIR=$1
lifetime=$2

ls $inputDIR*lhe > "filelist.txt"

while read line; do
	echo "lifetimefix.C+(\"${line}\",${lifetime})"
	root -l -n -q "lifetimefix.C+(\"${line}\",${lifetime})"
done < "filelist.txt"

echo "Done converting"