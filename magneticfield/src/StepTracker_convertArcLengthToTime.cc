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
#include "StepTracker_convertArcLengthToTime.hh"
#include <assert.h>

#include "G4ThreeVector.hh"

#define BUFFER_COLUMN_LEN 10


StepTracker_convertArcLengthToTime::StepTracker_convertArcLengthToTime() : StepTracker() {
}

StepTracker_convertArcLengthToTime::~StepTracker_convertArcLengthToTime() {
}

void StepTracker_convertArcLengthToTime::add_to_current_time( G4double h_to_add ) {

   // Right now assuming pure magnetic field, so velocity doesn't change.
   std::vector<double> last_vector;
   last_vector = get_buffer_ptr() -> back();

   velocity = G4ThreeVector(last_vector[4],last_vector[5],last_vector[6]).mag();

   assert( velocity != 0. );

   // h_to_add is passed as arc length to add, but we want to add time.
   // so we divide by velocity.
   StepTracker::add_to_current_time( h_to_add / velocity );


}
