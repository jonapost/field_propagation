# VecEM

## Added files:

* add files class diagram to `documents/`

* add file `TCashKarpRKF45.hh` to `em/include/`

* add file `TMagFieldEquation.hh` to `em/include/`

## Test:

* compile `make testPropagateMagField` at `test/`, 
	and run `testPropagateMagField 14`

* note that revise `G4BASE` and `G4GITBASE` before compiling


## Bugs:

* the test program crashed while using callgrind, but run normally without it

* I don't know the method to create `.cc` files for template classes. I simply put the source code in the header file. `.cc` files are in `em/src/backup`.
