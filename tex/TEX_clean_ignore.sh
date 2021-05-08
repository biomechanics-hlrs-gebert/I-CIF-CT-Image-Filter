#!/bin/bash
#
declare -a rry=("aux" "fdb_latexmk" "toc" "log" "fls" "snm" "synctex.gz" "nav" "out" "bbl" "blg")
#
for i in "${rry[@]}"
do
    rm *"."$i  >/dev/null 2>/dev/null &
done

