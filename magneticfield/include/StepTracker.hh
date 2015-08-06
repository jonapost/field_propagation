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
#include "G4FieldTrack.hh"

#ifndef G4MAGIntegratorSTEPPER
#include "G4MagIntegratorStepper.hh"
#else

class G4MagIntegratorStepper;

#endif

#include "isTracking.hh"

#include <vector>
using namespace std;

#define BUFFER_COLUMN_LEN 10


class StepTracker {
public:
   StepTracker(G4double beginning[BUFFER_COLUMN_LEN]);

   virtual ~StepTracker();

   G4bool check_that_wasnt_disgarded_by_Propagator(const G4double yIn[] );


   // Functions to handle book keeping of intersection points (thrown by the G4PropagatorInField class).
   inline void set_within_AdvanceChordLimited(G4bool status) { within_AdvanceChordLimited = status; }
   inline G4bool get_within_AdvanceChordLimited() { return within_AdvanceChordLimited; }

   void record_if_post_intersection_point(G4FieldTrack& possible_post_intersection_point,
                                          G4double passed_curve_length);

   inline void set_last_curve_length(G4double curve_length) { last_curve_length = curve_length; }
   inline void set_stepper_pointer( G4MagIntegratorStepper * stepper ) { myStepper = stepper; }

   inline void set_integrating_by_velocity(G4bool mbool) { integrating_by_velocity = mbool; }

   // End of intersection pt. book keeping Functions.


   void RecordResultOfStepper( G4double yIn0[], G4double dydx0[],
                               G4double yIn1[], G4double dydx1[],
                               G4int no_function_calls = -1);

   void update_time_arclength( G4double time_to_add, G4double arclength_to_add );

   G4double last_velocity(); // This is used by ChordFinder and MagIntDriver to perform
                             // the conversions from momentum coordinates to velocity
                             // coordinates. It's needed currently because neither
                             // ChordFinder or MagIntDriver do any record keeping.

   void outputBuffer(char *outfile_name,
                     char *meta_outfile_name,
                     // recording function call history is optional:
                     char *no_function_calls_outfile_name = 0,
                     char *indices_intersection_pts_filename = 0);

   // Inline Functions:

   inline void set_mass(G4double mass_of_particle ) { mass = mass_of_particle; }
   inline G4double get_mass() { return mass; }

   inline void ArmTracker() { armed = true; }
   // alert StepTracker that the next step should be stored. This is to avoid storing
   // steps taking within DistChord() by an aux stepper, or extra steps involved with
   // using Richardson Extrapolation. (Currently Richardson Extrapolation is not an
   // issue because the two half steps are performed before the real step, and thus are
   // overwritten by the real step. However it's probably best to just be careful.)

   inline void UnArmTracker() {armed = false; }
   inline bool isArmed() {return armed; }

   inline vector< vector<G4double> > * get_buffer_ptr() { return &buffer; }
   inline G4int getBufferLength() { return buffer.size(); }
   inline vector<G4int> *get_no_function_calls_buffer() { return &no_function_calls_buffer; }

   inline G4int get_thrown_away_steps() { return thrown_away_steps; }
   inline G4int get_used_steps() { return getBufferLength(); } // Used in a different context than getBufferLength()

   inline void set_last_time_val_was_accepted(G4bool val) { last_time_val_was_accepted = val; }

private:

   // These two bools don't need to be initialized (The first call
   //of the stepper is always from within AdvanceChordLimited() )

   G4bool integrating_by_velocity; // True if we are using MagIntegratorStepper_byTime
                                   // False if we are using MagIntegratorStepper_byArcLength
   G4bool within_AdvanceChordLimited;
   G4double last_curve_length;

   G4double time_left_over_from_intersection_pt_overshoot,
            arclength_left_over_from_intersection_pt_overshoot;

   G4MagIntegratorStepper *myStepper;

   vector< vector<G4double> > buffer;
   vector<G4int> no_function_calls_buffer;
   vector<G4int> indices_of_intersection_points;

   G4double mass;
   G4bool last_time_val_was_accepted, armed;
   G4int thrown_away_steps;
};

#endif /* MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_ */
