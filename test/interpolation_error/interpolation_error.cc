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


#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"

#include "ErrorComputer.hh"

using namespace std;
using namespace CLHEP;

#define NUMBER_INTERPOLATION_VARIABLES 6
#define BUFFER_COLUMN_LEN 11


int main(int argc, char *args[]) {

   G4int len_bufferA, len_bufferB;
   char *bufferA_filename, *bufferB_filename, *output_filename;

   // Parameter input:
   if (argc < 6) {
      cout << "You must supply 5 arguments. See the source file." << endl;
      return 1;
   }

   len_bufferA = atoi(args[1]);
   len_bufferB = atoi(args[2]);

   bufferA_filename = args[3];
   bufferB_filename = args[4];

   output_filename = args[5];

   G4double **bufferA = new G4double* [len_bufferA];
   for (int i = 0; i < len_bufferA; i ++) {
      bufferA[i] = new G4double[BUFFER_COLUMN_LEN];
   }
   G4double **bufferB = new G4double* [len_bufferB];
   for (int i = 0; i < len_bufferB; i ++) {
      bufferB[i] = new G4double[BUFFER_COLUMN_LEN];
   }

   ifstream bufferA_file(bufferA_filename, ios::binary | ios::in);

   for (int i = 0; i < len_bufferA; i ++) {
      bufferA_file.read( reinterpret_cast<char*>( bufferA[i] ),
                        BUFFER_COLUMN_LEN * sizeof(G4double) );
      }
   bufferA_file.close();

   ifstream bufferB_file(bufferB_filename, ios::binary | ios::in);

   for (int i = 0; i < len_bufferB; i ++) {
      bufferB_file.read( reinterpret_cast<char*>( bufferB[i] ),
                        BUFFER_COLUMN_LEN * sizeof(G4double) );
      }
   bufferB_file.close();

   ErrorComputer *mErrorComputer = new ErrorComputer( bufferA, len_bufferA,
                                                      bufferB, len_bufferB );

   G4double **err = new G4double* [len_bufferB];
   for (int i = 0; i < len_bufferB; i ++) {
      err[i] = new G4double[2 + NUMBER_INTERPOLATION_VARIABLES];
      // Right now, just record error for position
      // (time goes in the first component).
   }

   G4int no_interpolated_values = mErrorComputer -> ErrorArray(err);


   ofstream meta_outfile("interpolation_error_metafile", ios::out);
   meta_outfile << no_interpolated_values;
   meta_outfile.close();

   ofstream output(output_filename, ios::binary | ios::out);
   for (int i = 0; i < no_interpolated_values; i ++) {
      output.write( reinterpret_cast<char*>(err[i]),
                    (2 + NUMBER_INTERPOLATION_VARIABLES) * sizeof(G4double) );
   }
   output.close();

   // Clean up:

   for (int i = 0; i < len_bufferA; i ++)
         delete bufferA[i];
   for (int i = 0; i < len_bufferB; i ++) {
         delete bufferB[i];
         delete err[i];
   }
   delete bufferA; delete bufferB; delete err;

   delete mErrorComputer;

   return 0;
}






