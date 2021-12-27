#!/bin/bash
# ----------------------------------------------------------------------------------------
# Johannes Gebert - Doctoral project - initialize the subtree of the centralized sources
#
# Author:          Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:         27.12.2021
# Last edit:       27.12.2021
# ----------------------------------------------------------------------------------------
git subtree add --prefix central_src git@github.com:biomechanics-hlrs-gebert/A-CESO-Central_Sources.git  main --squash
