/*
 * StepTracker.cc
 *
 *  Created on: Jul 12, 2015
 *      Author: jason
 */

#include <StepTracker.hh>

#include <iostream>

#include <vector>
using namespace std;

#define BUFFER_COLUMN_LEN 10
#define TIME_SLOT 9

StepTracker::StepTracker() {

   buffer_ptr = new vector< vector<G4double> > ();
   cout << "allocated" << endl;
}

StepTracker::~StepTracker() {

   delete buffer_ptr;
}


void StepTracker::ReportCurveLength(G4double current_curve_length, G4double htry) {
   cout << "ReportCurveLength: " << current_curve_length
         << " htry: " << htry << endl;

}


void StepTracker::RecordResultOfStepper( G4double yIn[],
      G4double dydx[],
      G4double hstep,
      G4double yOutput[],
      G4double yError[] ) {


   cout << "yIn: ";
   for (int i = 0; i < BUFFER_COLUMN_LEN; i ++) {
      if (i < 3)
         cout << yIn[i] << "," ;
      last_y[i] = yIn[i];
   }
   cout << endl;

}

void StepTracker::StepsAccepted( G4double newCurveLength ) {

   vector< vector<G4double> > &buffer = *buffer_ptr;

   buffer_ptr -> push_back( vector<G4double> (BUFFER_COLUMN_LEN) );
   G4int last_index = buffer_ptr -> size() - 1;
   for (int i = 0; i < 3; i ++) {
      buffer[last_index][i] = last_y[i];
   }
   buffer[last_index][TIME_SLOT] = newCurveLength;




   cout << "StepAccepted, New Curve Length: " << newCurveLength << endl;

}
