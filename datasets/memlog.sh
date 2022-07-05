#!/bin/bash
# ----------------------------------------------------------------------------------------
# Memory logging
# Author:      Johannes Gebert gebert@hlrs.de
# Created: on: 30.01.2022
# Last edit:   13.02.2022
# ----------------------------------------------------------------------------------------
# You may need to switch between languages.
# ----------------------------------------------------------------------------------------
if [[ -z $1 ]]; then
    echo "Please specify a file to log to."
else
    # if [[ -f $1 ]]; then
    #     echo "The file $1 already exists. Please delete manually."
    # else
        while sleep 2
        do
            {
            # Grep for memory or swap (!)
#             free -h | grep "Mem" | tr -d '\n'
              free -h | grep "Speicher" | tr -d '\n'
              echo "     " | tr -d '\n'
              date
            } >> "$1"
        done
    # fi
fi

