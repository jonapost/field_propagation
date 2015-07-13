/*
 * StepTracker.hh
 *
 *  Created on: Jul 12, 2015
 *      Author: jason
 */

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
                           G4double dydx[],
                           G4double hstep,
                           G4double yOutput[],
                           G4double yError[] );

   void ReportCurveLength(G4double current_curve_length, G4double htry );

   void StepsAccepted( G4double newCurveLength );


private:
   vector< vector<G4double> > *buffer_ptr;
   G4double last_y[BUFFER_COLUMN_LEN];
};

#endif /* MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_ */
