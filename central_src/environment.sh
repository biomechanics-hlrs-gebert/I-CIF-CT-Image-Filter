#!/bin/bash
# ----------------------------------------------------------------------------------------
# Johannes Gebert - Doctoral project - set environment
#
# Author:          Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:         10.05.2021
# Last edit:       28.12.2021
# ----------------------------------------------------------------------------------------
# Update the subtree
#
which git > /dev/null 2> /dev/null
if  [ $? -eq 0 ] ; then 
    git subtree pull --prefix \
    central_src git@github.com:biomechanics-hlrs-gebert/A-CESO-Central_Sources.git \
    main --squash
fi
# ----------------------------------------------------------------------------------------
# Set colors for text output
#
if [ -z $SITE_NAME ]; then
    red='\033[0;31m'
    green='\033[0;32m'
    yellow='\033[0;33m'
    nc='\033[0m'
else
# Expected "HLRS" to be given for "SITE_NAME". 
# If SITE_NAME is not given, it is assumend, that the terminal understands colorizing :-)
    red=''
    green=''
    yellow=''
    nc=''
fi
# ----------------------------------------------------------------------------------------
# System environment files
#
sys_env_path='/central_src/auxiliaries/system_environments/'
# ----------------------------------------------------------------------------------------
#
usage ()
{
    echo "-- "
    echo "-- Usage: "
    echo "--      source environment.sh <system>"
    echo "-- "
    echo "-- Environments available:"
    echo "--      ${green}hawk${nc}   - HLRS HPE Apollo"
    echo "--      ${green}vulcan${nc} - HLRS NEC Cluster"
    echo "--      ${green}julius${nc} - A Whiskey Lake Notebook, 4 cores, 16Gb memory, APU"
    echo "-- "
    echo "--      Appending --no-output suppresses all output."
    echo "-- "
    echo "--------------------------------------------------------------------------------"
}
#
#------------------------------------------------------------------------------
prefix=$PWD
#
if [ -z $1 ]; then
    usage
else
    #
    if [ "$2" != "--no-output" ]; then
        echo "--------------------------------------------------------------------------------"
        echo "-- ${green}Setting environment${nc} for system: "$1
        echo "--"
        export NO_OUTPUT=NO
    else
        export NO_OUTPUT=YES
    fi
    #
    sys_set=0
    for sys_file in $(ls --color=never ${prefix}${sys_env_path})
    do
        system=$(basename -s .sh $sys_file)
        #
        test $system = $1 && source ${prefix}${sys_env_path}${sys_file} && sys_set=1
        #
        # System
        export SYS_ENV=$1
    done

    if [ $sys_set -eq 0 ]; then
       echo "--"
       echo "-- ${yellow}System ${red}$1 ${yellow}currently is not supported.${nc}"
       usage
    else
	#
	# ----------------------------------------
	# PATH extensions
	export PATH=${prefix}/bin:$PATH
	#
    if [ "NO_OUTPUT" != "YES" ]; then
	    echo "-- ${green}Done${nc}"
	    echo "--------------------------------------------------------------------------------"
	fi
    fi
fi