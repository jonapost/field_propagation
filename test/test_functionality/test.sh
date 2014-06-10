#!/bin/sh
#
# A simple script to run all the tests in this directory and check
# their results against the expected (previous) results
#
# $Id: test.sh 66241 2012-12-13 18:34:42Z gunter $
# $Name: geant4-09-02-ref-05 $
#

echo "Running on `hostname`, which is a `uname -a` machine" 
host=`hostname`

#
MAKE=make
target=testPropagateMagField
echo  "Compiling $target ... "
$MAKE -s  $target   
echo  "Executing $target ..."
for n in 8 
do
  ./$target $n > $target.newout$n \
		     2> $target.newerr$n
  echo  ".. difference from expected output: "
  diff -wb $target.out$n $target.newout$n
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $target.err$n $target.newerr$n
  sleep 1;
  echo  " "
done

exit
