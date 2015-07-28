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

//#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "G4ThreeVector.hh"

#define BUFFER_COLUMN_LEN 11
#define POSITION_SLOT 2
#define MOMENTUM_SLOT 5
#define RHS_SLOT 8

StepTracker::StepTracker(G4double beginning[BUFFER_COLUMN_LEN]) {
   // beginning must be of the form
   // {time, arclength, position[1..3], momentum[1..3], RHS[1..3]}

   last_time_val_was_accepted = false;
   armed = false;

   thrown_away_steps = 0;

   buffer.push_back( vector<G4double> (BUFFER_COLUMN_LEN) );
   no_function_calls_buffer.push_back( 0 );

   for (int i = 0; i < BUFFER_COLUMN_LEN; i ++)
      buffer[0][i] = beginning[i];

}

StepTracker::~StepTracker() {
}

void StepTracker::outputBuffer(char *outfile_name,
                               char *no_function_calls_outfile_name,
                               char *meta_outfile_name) {

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

   if (no_function_calls_outfile_name != 0) {
      ofstream no_function_calls_outfile( no_function_calls_outfile_name, ios::out );
      for (int i = 0; i < getBufferLength(); i ++)
         no_function_calls_outfile << no_function_calls_buffer[i] << endl;

      no_function_calls_outfile.close();
   }
}

void StepTracker::RecordResultOfStepper( G4double yIn[], G4double dydx[], G4int no_function_calls) {
   // Time is stored in first component.

   // Create new space for next time/position/momentum values:
   if ( last_time_val_was_accepted ) {

      buffer.push_back( vector<G4double> (BUFFER_COLUMN_LEN) );
      if (no_function_calls != -1)
         no_function_calls_buffer.push_back( no_function_calls );

      last_time_val_was_accepted = false;
   }
   // Otherwise, just write over the last position/momentum values
   else {
      thrown_away_steps ++;
   }

   // Take care of beginning. (Don't overwrite the first row.):
   if ( buffer.size() == 1 ) {
      buffer.push_back( vector<G4double> (BUFFER_COLUMN_LEN) );
      if (no_function_calls != -1)
         no_function_calls_buffer.push_back( no_function_calls );
   }

   G4int last_index = buffer.size() - 1;

   for (int i = 0; i < 6; i ++) {
      buffer[last_index][i + POSITION_SLOT] = yIn[i];
   }
   for (int i = 0; i < 3; i ++) {
      buffer[last_index][i + RHS_SLOT] = dydx[i + 3];
   }
}

G4bool StepTracker::check_that_wasnt_disgarded_by_Propagator(
      const G4double yIn[] ) {

   G4int last_index = getBufferLength() - 1;

   for (int i = 0; i < 3; i ++) {
      if (yIn[i] != buffer[last_index][i + POSITION_SLOT]) {
         return false;
      }
   }
   return true;
}

G4double StepTracker::last_velocity() {
   G4int last_index = getBufferLength() - 1;
   return G4ThreeVector( buffer[last_index][MOMENTUM_SLOT],
                         buffer[last_index][MOMENTUM_SLOT + 1],
                         buffer[last_index][MOMENTUM_SLOT + 2] ).mag();
}


void StepTracker::update_time_arclength( G4double time_to_add, G4double arclength_to_add) {

   // Right now assuming pure magnetic field, so velocity doesn't change.

   // If update_time_arclength() is called from context of position/momentum coordinates,
   // not position/velocity coordinates, it is the responsibility of the caller to
   // convert to position/velocity coordinates ( divide momentum coordinates
   // through by velocity).

   last_time_val_was_accepted = true;
   G4int last_index = buffer.size() - 1;

   buffer[last_index][0] = buffer[last_index - 1][0] + time_to_add;        // time
   buffer[last_index][1] = buffer[last_index - 1][1] + arclength_to_add;   // arc_length
}


