# Field Propagation for Geant4 (and soon GeantV)

Started from two sources:
* the Geant4 magneticfield category / directory
* the refinements / improvements created by Quichen Xie (GSoC 2014), supervised by Sandro Wenzel, aided by John Apostolakis

Recent Contributors:
* John Apostolakis (CERN)
* Somnath Banerjee (GSoC 2015 contributor)
* Jason Suagee (GSoC 2015 contributor)
* Sandro Wenzel (CERN)
* Quichen Xie (GSoC 2014 contributor)

Goals (GSoC 2015):
* Add new integration steppers (integrators), which implement FSAL (first same as last) Runge-Kutta methods
* Add modern integrators with higher efficiency (DormandPrince 5/4, Bogacki Shampine 5/4, Verner)
* Add integrators with interpolants - i.e. which are capable to interpolate values of the integrand(s) at intermediate points in an interval
* Improve the facility / design to enable the use of general Nystrom RK methods - ie ones ODEs which have involve the first derivative in the right hand side;
* Add embedded Nystrom steppers with error estimation
* Add Nystrom methods with interpolants.

## Test:

* /test directory
* create a [functionality report](https://bitbucket.org/jonapost/field_propagation/wiki/functionality) 

## BenchMark:

* /benchmark directory
* read the [benchmark report](https://bitbucket.org/jonapost/field_propagation/wiki/benchmark)
