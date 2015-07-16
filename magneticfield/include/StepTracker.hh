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

#ifndef MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_
#define MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_

#include "G4Types.hh"

#include <vector>
using namespace std;

#define BUFFER_COLUMN_LEN 10


class StepTracker {
public:
   StepTracker();
   virtual ~StepTracker();

   void RecordResultOfStepper( G4double yIn[],
                           G4double dydx[]);

   void ReportCurveLength(G4double current_curve_length, G4double htry );

   void StepsAccepted( G4double newCurveLength );

   inline G4int getBufferLength() { return buffer_ptr -> size(); }

   void outputBuffer(char *outfile_name);

private:
   vector< vector<G4double> > *buffer_ptr;
   //G4double last_y[BUFFER_COLUMN_LEN];
};

#endif /* MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_ */
