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

#include <assert.h>

#include "G4ThreeVector.hh"

#define BUFFER_COLUMN_LEN 11
#define POSITION_SLOT 2
#define MOMENTUM_SLOT 5
#define RHS_SLOT 8

StepTracker::StepTracker(G4double beginning[BUFFER_COLUMN_LEN]) { // beginning must be of the form
                                                 // {time, position[3], momentum[3], RHS[3]}

   buffer_ptr = new vector< vector<G4double> > ();

   last_time_val_was_accepted = false;
   armed = false;

   thrown_away_steps = 0;

   buffer_ptr -> push_back( vector<G4double> (BUFFER_COLUMN_LEN) );
   vector< vector<G4double> > &buffer = *buffer_ptr;

   for (int i = 0; i < BUFFER_COLUMN_LEN; i ++)
      buffer[0][i] = beginning[i];

}

StepTracker::~StepTracker() {

   delete buffer_ptr;
}

vector< vector<G4double> > * StepTracker::get_buffer_ptr() { return buffer_ptr;}

void StepTracker::outputBuffer(char *outfile_name, char *meta_outfile_name) {
   vector< vector<G4double> > &buffer = *buffer_ptr;

   ofstream meta_outfile(meta_outfile_name, ios::out);
   meta_outfile << getBufferLength();
   meta_outfile.close();

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

   // Gets called from

   // cout << "ReportCurveLength: " << current_curve_length
   //      << " htry: " << htry << endl;

}


//void StepTracker::start_next_at_beginning() {
//}

void StepTracker::RecordResultOfStepper( G4double yIn[], G4double dydx[]) {

   //cout << "RecordResultOfStepper" << endl;

   vector< vector<G4double> > &buffer = *buffer_ptr;

   // Only record data from the beginning of the interval.
   // This is because there is no pre-set way to pass the last RHS
   // evaluation to StepTracker from the BaseStepper of MagIntegratorStepper_byTime. (Yet.)

   //Time is stored in first component.

   if ( last_time_val_was_accepted ) {       // Create new space for next time/position/momentum values

      buffer_ptr -> push_back( vector<G4double> (BUFFER_COLUMN_LEN) );
      last_time_val_was_accepted = false;
   }
   else {
      thrown_away_steps ++;
   }
                                             // Otherwise, just write over the last position/momentum values

   if ( buffer_ptr -> size() == 1 )          // Take care of beginning. (Don't overwrite the first row.)
      buffer_ptr -> push_back( vector<G4double> (BUFFER_COLUMN_LEN) );

   G4int last_index = buffer_ptr -> size() - 1;

   // buffer[last_index][0] = yIn[TIME_SLOT];

   for (int i = 0; i < 6; i ++) {
      buffer[last_index][i + POSITION_SLOT] = yIn[i];
   }
   for (int i = 0; i < 3; i ++) {
      buffer[last_index][i + RHS_SLOT] = dydx[i + 3];
   }

}

G4double StepTracker::last_velocity() {
   vector< vector<G4double> > &buffer = *buffer_ptr;
   G4int last_index = getBufferLength() - 1;
   return G4ThreeVector( buffer[last_index][MOMENTUM_SLOT],
                         buffer[last_index][MOMENTUM_SLOT  +1],
                         buffer[last_index][MOMENTUM_SLOT + 2] ).mag();

}


void StepTracker::add_to_current_time( G4double time_to_add, G4double arclength_to_add) {

   // Right now assuming pure magnetic field, so velocity doesn't change.
   //std::vector<double> last_vector;
   //last_vector = get_buffer_ptr() -> back();

   // add_to_current_time() is called from context of position/momentum coordinates,
   // not position/velocity coordinates. So we have to divide by the mass.
   // G4double velocity = G4ThreeVector(last_vector[4],last_vector[5],last_vector[6]).mag() / get_mass();

   // assert( velocity != 0. );

   vector< vector<G4double> > &buffer = *buffer_ptr;

   last_time_val_was_accepted = true;
   G4int last_index = buffer_ptr -> size() - 1;

   /* Don't need this anymore:
   if (last_index == 0) { // Takes care of case when we are at the very beginning.
      buffer[last_index][0] = 0 + h_to_add; // / velocity;
   }
   else {
      buffer[last_index][0] = buffer[last_index - 1][0] + h_to_add; // / velocity;
   }
   */

   buffer[last_index][0] = buffer[last_index - 1][0] + time_to_add;        // time
   buffer[last_index][1] = buffer[last_index - 1][1] + arclength_to_add;   // arc_length
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
