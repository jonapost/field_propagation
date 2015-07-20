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

#include "StepTracker.hh"

#include <iostream>
#include <fstream>

#include <vector>
using namespace std;

#define BUFFER_COLUMN_LEN 10
#define TIME_SLOT 7

StepTracker::StepTracker() {

   buffer_ptr = new vector< vector<G4double> > ();
   //cout << "allocated" << endl;
}

StepTracker::~StepTracker() {

   delete buffer_ptr;
}


void StepTracker::outputBuffer(char *outfile_name) {

   vector< vector<G4double> > &buffer = *buffer_ptr;

   ofstream outfile;
   outfile.open(outfile_name, ios::binary | ios::out);
   for (int i = 0; i < getBufferLength(); i ++) {
      for (int j = 0; j < BUFFER_COLUMN_LEN; j ++){
         outfile.write( reinterpret_cast<char*>( &(buffer[i][j]) ), sizeof( buffer[i][j] ) );
      }
   }
   outfile.close();
   /*
   for (int i = 0; i < getBufferLength(); i ++) {
      for (int j = 0; j < BUFFER_COLUMN_LEN; j ++)
         cout << buffer[i][j] << ",";
      cout << endl;
   }
   */
}

void StepTracker::ReportCurveLength(G4double current_curve_length, G4double htry) {
   cout << "ReportCurveLength: " << current_curve_length
         << " htry: " << htry << endl;

}


void StepTracker::RecordResultOfStepper( G4double yIn[], G4double dydx[]) {

   vector< vector<G4double> > &buffer = *buffer_ptr;

   // Only record data from the beginning of the interval.
   // This is because there is no pre-set way to pass the last RHS
   // evaluation to StepTracker from the BaseStepper of MagIntegratorStepper_byTime. (Yet.)


   //Time is stored in first component.

   buffer_ptr -> push_back( vector<G4double> (BUFFER_COLUMN_LEN) );
   G4int last_index = buffer_ptr -> size() - 1;
   buffer[last_index][0] = yIn[TIME_SLOT];
   for (int i = 0; i < 6; i ++) {
      buffer[last_index][i + 1] = yIn[i];
   }
   for (int i = 3; i < 6; i ++) {
      buffer[last_index][i + 4] = dydx[i];
   }
   /*
   cout << "yIn: ";
   for (int i = 0; i < BUFFER_COLUMN_LEN; i ++) {
      if (i < 3)
         cout << yIn[i] << "," ;
      last_y[i] = yIn[i];
   }
   cout << endl;
   */
}

void StepTracker::StepsAccepted( G4double newCurveLength ) {

   /*
   vector< vector<G4double> > &buffer = *buffer_ptr;

   buffer_ptr -> push_back( vector<G4double> (BUFFER_COLUMN_LEN) );
   G4int last_index = buffer_ptr -> size() - 1;
   for (int i = 0; i < 3; i ++) {
      buffer[last_index][i] = last_y[i];
   }
   buffer[last_index][TIME_SLOT] = newCurveLength;




   cout << "StepAccepted, New Curve Length: " << newCurveLength << endl;
   */
}
