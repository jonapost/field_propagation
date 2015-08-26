# Field Propagation for Geant4 (and soon GeantV)

Recent Contributors:
* John Apostolakis (CERN)
* Somnath Banerjee (GSoC 2015 contributor)
* Jason Suagee (GSoC 2015 contributor)
* Sandro Wenzel (CERN)
* Quichen Xie (GSoC 2014 contributor)

Summary of Goals (GSoC 2015):
* Add integration Runge-Kutta methods (integrators), with additional capabilities including FSAL (first same as last) and interpolants (method to evalue the integrand at an intermediate point);
* Implement modern integrators with higher efficiency (DormandPrince 5/4, Bogacki Shampine 5/4, Verner);
* Introduce new general Nystrom RK methods - ie ones for ODEs with the first derivative in the right hand side;
* Enable vectorisation for one or multiple tracks.

An extended description of the 2015 GSoC project goals can be found at the [Wiki page](https://github.com/jonapost/field_propagation/wiki)

Started from two sources:
* the Geant4 magneticfield category / directory
* the refinements / improvements created by Quichen Xie (GSoC 2014), supervised by Sandro Wenzel, aided by John Apostolakis

## Nystrom-dev branch

* See documents/Status.md for an update on the status of the work on this branch, and 
* documents/caveats.md for a lists of important things to be aware of when using or modifying this work.

## Test:

* /test directory
* create a [functionality report](https://bitbucket.org/jonapost/field_propagation/wiki/functionality) 

## BenchMark:

* /benchmark directory
* read the [benchmark report](https://bitbucket.org/jonapost/field_propagation/wiki/benchmark)
