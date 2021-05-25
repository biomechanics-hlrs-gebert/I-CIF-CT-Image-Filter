#!/bin/bash
#
declare -a rry=("aux" /
"fdb_latexmk" /
"toc" /
"log" /
"fls" /
"snm" /
"synctex.gz" /
"nav" /
"out" /
"bbl" /
"blg" /
"pdflatex_compile__.log" /
"texput.log")
#
for i in "${rry[@]}"
do
    rm $PWD/*"."$i  >/dev/null 2>/dev/null &
done

