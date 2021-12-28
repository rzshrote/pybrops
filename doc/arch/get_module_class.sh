#!/usr/bin/env bash

# output file name
outfile="file_list.txt"

# write header
printf "%s\t%s\n" "Module" "Class" > $outfile

# 
find ../pybropt/ -name '*.py' -print0 | while read -d $'\0' f
do
fpath="${f%/*}"
tmp="${fpath##..\/pybropt\/}"
mname="${tmp//\//.}"
fname="${f##*/}"
bname="${fname%.*}"
printf "%s\t%s\n" $bname $mname
done | grep --invert-match '__init__' >> $outfile

