#!/bin/bash

inputDIR=$1

ls $inputDIR*lhe > "filelist.txt"

while read line; do
	root -l -n -q "leptonmassfix.C+(\"${line}\")"
done < "filelist.txt"

echo "Done converting"