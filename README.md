# VecEM

## Added files:

* add files class diagram to `documents/`
* add file `TCashKarpRKF45.hh` to `em/include/`
* add file `TMagFieldEquation.hh` to `em/include/`

## Test:

* compile `make testPropagateMagField` at `test/`, 
	and run `testPropagateMagField 14`
* note that revise `G4BASE` and `G4GITBASE` before compiling

## Profiling:

for the testing example, testPropagateMagField:

* TCashKarpRKF45::Stepper(): 66%;
* G4CachedMagneticField::GetFieldValue():11%;

time (run 10 times main() function several times):

* TCashKarpRKF45: 0.448s +- 0.003s
* G4CashKarpRKF45: 0.470s +- 0.004s
(note: I turned off major STDOUT, i.e verbose and report.)


## Issues:

