#!/bin/sh
#
# A simple script to run all the tests in this directory and check
# their results against the expected (previous) results
#
# $Id: test.sh 105866 2017-08-24 13:29:35Z japost $
# $Name: geant4-09-02-ref-05 $
#

# Find the Geant4 environment
if [[ "x" != "x$G4BUILD_DIR" ]]; then
   . $G4BUILD_DIR/geant4make.sh
fi

echo "Running on `hostname`, which is a `uname -a` machine" 
host=`hostname`

#  Choice of build engine
 MAKE=make
#MAKE=ninja
#####

#  Choice of build target directory BUILD_DIR
BUILD_DIR=$G4BUILD_DIR
echo $BUILD_DIR
#G4BIN=$G4WORKDIR/bin/
#G4BIN=$G4BIN
#G4BIN=$G4BIN_GMAKE/bin/
G4BIN=$G4BIN_CMAKE/

#   Location of executable binaries  BINDIR
#BINDIR=$G4BIN/$G4SYSTEM/
BINDIR=$G4BIN/

target=testPropagateMagField
echo  "Compiling $target ... "
( cd $BUILD_DIR ; $MAKE $target )  || exit
if [[ ! -x $BINDIR/$target ]] ; then
  echo "Could not find executable $target in directory $BINDIR" 
  exit 1
fi
echo  "Executing $target ..."
for n in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 23 45 56 78 145 745
do
  echo "Executing with stepper choice $n ..."
  $BINDIR/$target $n  > $target.newout$n \
		     2> $target.newerr$n
  if [[ -f $target.out$n ]]; then
     if [[ `cmp --silent $target.out$n $target.newout$n` ]]; then
	 echo  ".. difference(s) from expected output: "
	 diff -wb $target.out$n $target.newout$n
	 sleep 2;
     fi
  else
     echo  " Expected output *not* found. Time to create " $target.out$n
  fi
  # if [[ -f $target.err$n ]]; then
  if [[ -f $target.newerr$n ]]; then
      if [[ `wc -l $target.newerr$n | awk ' { print $1; } '` != 0 ]]; then
	  echo  ".. Unexpected error output: "
	  if [[ -f $target.err$n ]]; then
	      diff -wb $target.err$n $target.newerr$n
	  else
	      head -50 $target.newerr$n
	  fi
	  sleep 5;
      else
	  ## echo " Deleting empty error file ... "
	  rm $target.newerr$n
      fi
  else
      echo "WARNING> 'Error' output file $target.newerr$n not found."
  fi      
  echo  " "
done

exit

target=testProElectroMagField
echo  "Compiling $target ... "
( cd $BUILD_DIR ; $MAKE $target )   || exit
echo  "Executing $target ..."
for n in 1 2 3 4 8   23 45 56 78 145 745
do
  echo "Executing with stepper choice $n ..."
  $BINDIR/$target $n  > $target.newout$n \
		     2> $target.newerr$n
  echo  ".. difference from expected output: "
  diff -wb $target.out$n $target.newout$n
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $target.err$n $target.newerr$n
  sleep 1;
  echo  " "
done

for i in *Spin.cc
do
  target=`basename $i .cc`
  echo  "Compiling $target ... "
  ( cd $BUILD_DIR ; $MAKE $target )  || exit
  echo  "Executing $target ..."
  for n in  4 3 2 1 0  # 23 45 56 78 145 745
  do
    echo "Executing with stepper choice $n ..."
    $BINDIR/$target $n  > $target.newout$n \
		       2> $target.newerr$n
    echo  ".. difference from expected output: "
    diff -wb $target.out$n $target.newout$n
    sleep 1;
    echo  ".. difference from expected error: "
    diff -wb $target.err$n $target.newerr$n
    sleep 1;
    echo  " "
  done
done

exit

for i in *Spin.cc
do
  target=`basename $i .cc`
  echo  "Compiling $target ... "
  gmake G4TARGET=$target || exit
  echo  "Executing $target ..."
  $BINDIR/$target $n  > $target.newout$n \
		     2> $target.newerr$n
  echo  ".. difference from expected output: "
  diff -wb $target.out $target.newout
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $target.err $target.newerr
  sleep 1;
  echo  " "
done

exit
