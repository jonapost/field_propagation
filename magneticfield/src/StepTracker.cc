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

#ifndef G4CACHED_MAGNETIC_FIELD_DEF
#include "G4CachedMagneticField.hh"
#endif


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

   last_curve_length = 0.;

}

StepTracker::~StepTracker() {
}


void StepTracker::record_if_post_intersection_point( G4FieldTrack& possible_post_intersection_point,
                                                     G4double passed_curve_length ) {

   if ( passed_curve_length < last_curve_length ) {
      // passed_curve_length is the curve_length of the prev
      // Here is where we do something to record in the buffer that we had an intersection pt.

      // Note: this does not treat the case where the last point returned by AdvanceChordLimited()
      // was actually an intersection point (exactly). We might want to also consider the used
      // number of function evaluations to distinguish if there was no intersection in
      // PropagatorInField::ComputeStep(). If there was the code dealing with finding the intersection
      // points uses function evaluations (or at least function calls) (Is this true??)

      cout << " Possible Intersection pt" << endl;


      G4double y[G4FieldTrack::ncompSVEC];
      G4double dydx[G4FieldTrack::ncompSVEC];


      possible_post_intersection_point.DumpToArray( y );

      // Get second derivatives (RHS function evaluation):
      myStepper -> ComputeRightHandSide( y, dydx );

      const G4CachedMagneticField *myField = (G4CachedMagneticField*)
                                 ( myStepper -> GetEquationOfMotion() -> GetFieldObj() );

      if ( integrating_by_velocity ) { // If using MagIntegratorStepper_byTime we have to scale
                                       // by 1/mass. This usually gets done from within
                                       // MagIntegratorStepper_byTime, but
                                       // record_if_post_intersection_point() will always be
                                       // called from AdvanceChordLimited().
         for (int i = 3; i < 6; i ++) {
            y[i] /= mass;
            dydx[i] /= mass;
            // TODO: Eventually change mass to myField->FMass(), but not necessary now
            // since relatavistic mass doesn't change in a pure magnetic field.
         }
      }

      G4int no_function_calls = myField -> GetCountCalls();

      buffer.pop_back(); // Get rid of overshoot point (overshot the intersection point).
      no_function_calls_buffer.pop_back(); // Similar as line above.

      RecordResultOfStepper( y, dydx, no_function_calls );

      update_time_arclength(
            (1. / last_velocity() ) *
            ( passed_curve_length - buffer[buffer.size() - 1][1] ),  // Must convert arc length difference
                                                                     // to time units, so divide by velocity
                passed_curve_length - buffer[buffer.size() - 1][1]   // Arc length to add is passed arc length
                                                                     // minus previous arc length (which is in
                                                                     // the last entry of buffer).
                                                            ) ;

      indices_of_intersection_points.push_back( getBufferLength() );

   }
   // Either way we don't want to erase the last row of buffer data.
   last_time_val_was_accepted = true;
}

void StepTracker::outputBuffer(char *outfile_name,
                               char *meta_outfile_name,
                               char *no_function_calls_outfile_name,
                               char *indices_intersection_pts_filename) {

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

   if (indices_intersection_pts_filename != 0) {
      ofstream indices_intersection_pts_outfile( indices_intersection_pts_filename, ios::out );
      for (int i = 0; i < indices_of_intersection_points.size(); i ++)
         indices_intersection_pts_outfile << indices_of_intersection_points[i] << endl;

      indices_intersection_pts_outfile.close();
   }

}

void StepTracker::RecordResultOfStepper( G4double yIn[], G4double dydx[], G4int no_function_calls) {
   // Time is stored in first component.

   if ( ! within_AdvanceChordLimited ) {
      cout << "Not within AdvancedChordLimited" << endl;
      return;
   }

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

   if ( ! within_AdvanceChordLimited ) {
      cout << "Not within AdvancedChordLimited" << endl;
      return;
   }

   last_time_val_was_accepted = true;
   G4int last_index = buffer.size() - 1;

   buffer[last_index][0] = buffer[last_index - 1][0] + time_to_add;        // time
   buffer[last_index][1] = buffer[last_index - 1][1] + arclength_to_add;   // arc_length
}


