#!/bin/bash

inputDIR=$1

ls $inputDIR*lhe > "filelist.txt"

while read line; do
	root -l -n -q "scalefix.C+(\"${line}\")"
done < "filelist.txt"

echo "Done converting"