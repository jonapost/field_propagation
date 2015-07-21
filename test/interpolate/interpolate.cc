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

#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


#include "G4QuadrupoleMagField.hh"
#include "G4CachedMagneticField.hh"

#include "G4CashKarpRKF45.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "BogackiShampine23.hh"
#include "G4LineSection.hh"
#include "G4MagIntegratorStepper.hh"
#include "DormandPrince745.hh"
#include "BogackiShampine45.hh"
#include "G4ClassicalRK4.hh"
#include "G4SimpleHeum.hh"
#include "G4ChargeState.hh"

#include "G4NystromRK4.hh"

#include "Mag_UsualEqRhs_IntegrateByTime.hh"
#include "ChawlaSharmaRKNstepper.hh"

#include "MagIntegratorStepperbyTime.hh"
#include "FineRKNG34.hh"
#include "FineRKNG45.hh"

#include <iostream>
#include <fstream>

#include "G4ThreeVector.hh"

#include "Interpolant.hh"

using namespace std;
using namespace CLHEP;

#define BUFFER_COLUMN_LEN 10


int main(int argc, char *args[]) {

   G4int stepper_to_interpolate_no, stepper_to_compare_no, partition_size, no_input_pairs;
   char *infile_to_interpolate_name, *outfile_name;

   // Not implemented yet
   // char *infile_to_compare_name;

   // Parameter input:
   if (argc < 5) {
      cout << "You must supply 4 arguments. See the source file." << endl;
      return 1;
   }
   // Isn't used:
   //stepper_to_interpolate_no = atoi(args[1]);
   //stepper_to_compare_no = atoi(args[2]);

   no_input_pairs = atoi(args[1]);
   partition_size = atoi(args[2]);

   infile_to_interpolate_name = args[3];
   outfile_name = args[4];

   // Not implemented yet.
   // infile_to_compare_name = args[5];

   G4double **y0 = new G4double* [no_input_pairs];
   G4double **y1 = new G4double* [no_input_pairs];
   G4double **F0 = new G4double* [no_input_pairs];
   G4double **F1 = new G4double* [no_input_pairs];
   for (int i = 0; i < no_input_pairs; i ++) {
      y0[i] = new G4double[6];
      y1[i] = new G4double[6];
      F0[i] = new G4double[3];
      F1[i] = new G4double[3];
   }
   G4double *step = new G4double[no_input_pairs];

   G4double t0, t1; //, step_discard;

   ifstream infile_to_interpolate(infile_to_interpolate_name, ios::binary | ios::in);


   for (int i = 0; i < no_input_pairs; i ++) {
      infile_to_interpolate.read( (char*) &t0, sizeof(G4double) );

      infile_to_interpolate.read( (char*) y0[i], 6 * sizeof(G4double) );
      infile_to_interpolate.read( (char*) F0[i], 3 * sizeof(G4double) );

      infile_to_interpolate.read( (char*) &(t1), sizeof(G4double) );

      //infile_to_interpolate.read( (char*) &step_discard, sizeof(G4double) );
      infile_to_interpolate.read( (char*) y1[i], 6 * sizeof(G4double) );
      infile_to_interpolate.read( (char*) F1[i], 3 * sizeof(G4double) );

      step[i] = t1 - t0;
   }

   infile_to_interpolate.close();

   G4double ***interpolated_values = new G4double** [no_input_pairs];
   for (int i = 0; i < no_input_pairs; i ++) {
      interpolated_values[i] = new G4double* [partition_size];
      for (int j = 0; j < partition_size; j ++)
         interpolated_values[i][j] = new G4double[BUFFER_COLUMN_LEN];
   }

   Interpolant *mInterpolant = new Interpolant();
   G4double accumulated_time = 0.;
   G4double xi, step_len;

   for (int i = 0; i < no_input_pairs; i ++) {

      mInterpolant -> Initialize(y0[i], y1[i], F0[i], F1[i], step[i]);

      for (int j = 0; j < partition_size; j ++) {
         step_len = step[i] / partition_size;
         xi = j * step_len / step[i];
         // Time stored in first index ( [0] ).
         interpolated_values[i][j][0] = xi * step_len + accumulated_time;
         mInterpolant -> InterpolatePosition(xi, &( interpolated_values[i][j][1] ) );
         mInterpolant -> InterpolateMomentum(xi, &( interpolated_values[i][j][4] ) );
      }
      accumulated_time += step[i];
   }

   ofstream outfile(outfile_name, ios::binary | ios::out);

   for (int i = 0; i < no_input_pairs; i ++) {
      for (int j = 0; j < partition_size; j ++) {
         outfile.write( (char*) interpolated_values[i][j], BUFFER_COLUMN_LEN * sizeof(G4double) );
      }
   }

   outfile.close();


   for (int i = 0; i < no_input_pairs; i ++) {
         delete y0[i];
         delete y1[i];
         delete F0[i];
         delete F1[i];
   }
   delete y0; delete y1; delete F0; delete F1;
   delete step;

   for (int i = 0; i < no_input_pairs; i ++) {
      for (int j = 0; j < partition_size; j ++)
         delete interpolated_values[i][j];
      delete interpolated_values[i];
   }
   delete interpolated_values;

   return 0;
}




