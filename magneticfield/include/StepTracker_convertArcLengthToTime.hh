// Nystrom stepper implemenation by Jason Suagee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 27 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
#ifndef MAGNETICFIELD_INCLUDE_STEPTRACKER_CONVERTARCLENGTHTOTIME_HH_
#define MAGNETICFIELD_INCLUDE_STEPTRACKER_CONVERTARCLENGTHTOTIME_HH_

#include "StepTracker.hh"

class StepTracker_convertArcLengthToTime: public StepTracker {
public:
   StepTracker_convertArcLengthToTime();
   virtual ~StepTracker_convertArcLengthToTime();

   void add_to_current_time( G4double h_to_add );


private:
   G4double velocity;

};

#endif /* MAGNETICFIELD_INCLUDE_STEPTRACKER_CONVERTARCLENGTHTOTIME_HH_ */
