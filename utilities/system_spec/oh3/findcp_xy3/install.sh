#!/bin/bash

# set environmental variables
source ../../../../bin/setsgenvars.sh

echo Compiling findcp
echo $SGENFC -g findcp.f90 -o findcp.x $SGENFLAG $SGENLIB
$SGENFC -g findcp.f90 -o findcp.x $SGENFLAG $SGENLIB  
echo Cleaning up
rm opttools.mod
echo Done 
