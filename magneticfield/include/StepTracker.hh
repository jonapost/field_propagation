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
   StepTracker(G4double beginning[BUFFER_COLUMN_LEN]);
   virtual ~StepTracker();

   G4double last_velocity();

   inline void set_mass(G4double mass_of_particle ) { mass = mass_of_particle; }
   inline G4double get_mass() { return mass; }

   inline void ArmTracker() { armed = true; }
   inline void UnArmTracker() {armed = false; }
   inline bool isArmed() {return armed; }

   void RecordResultOfStepper( G4double yIn[],
                           G4double dydx[]);

   void ReportCurveLength(G4double current_curve_length, G4double htry );

   virtual void add_to_current_time( G4double time_to_add, G4double arclength_to_add );

   void StepsAccepted( G4double newCurveLength );

   inline G4int getBufferLength() { return buffer_ptr -> size(); }

   void outputBuffer(char *outfile_name, char *meta_outfile_name);

   inline G4int get_thrown_away_steps() { return thrown_away_steps; }
   inline G4int get_used_steps() { return getBufferLength(); } // Used in a different context than getBufferLength()

   vector< vector<G4double> > * get_buffer_ptr();

private:
   vector< vector<G4double> > *buffer_ptr;

   G4double mass;

   G4bool last_time_val_was_accepted, armed;

   G4int thrown_away_steps;
   //G4double last_y[BUFFER_COLUMN_LEN];
};

#endif /* MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_ */
