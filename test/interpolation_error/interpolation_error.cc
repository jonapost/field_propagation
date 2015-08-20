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

#include <assert.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "G4ThreeVector.hh"

#include "ErrorComputer2.hh"

using namespace std;
using namespace CLHEP;


#define TIME_SLOT 0
#define ARCLENGTH_SLOT 1

#ifdef INTENDED_FOR_ERROR_BY_STEPPER_PROGRAM

#define BUFFER_COLUMN_LEN 28 // room for start point and end point of each step
                             // plus time/arclength entries for each.

#define ENDPOINT_BASE_INDEX 14
#define POSITION_SLOT 2
#define MOMENTUM_SLOT 5
#define RHS_SLOT 8


#else
#define BUFFER_COLUMN_LEN 22 // room for start point and end point of each step
                             // plus time/arclength entries for each.

#define ENDPOINT_BASE_INDEX 11
#define POSITION_SLOT 2
#define MOMENTUM_SLOT 5
#define RHS_SLOT 8
#endif

#define NUMBER_INTERPOLATION_VARIABLES 6


int main(int argc, char *args[]) {

   G4int len_bufferA, len_alt_buffer, len_bufferB;
   char *bufferA_filename, *alt_buffer_filename, *bufferB_filename, *output_filename;

   // Parameter input:
   if (argc < 6) {
      cout << "You must supply 7 arguments. See the source file." << endl;
      return 1;
   }

   len_bufferA = atoi(args[1]);
   len_alt_buffer = atoi(args[2]);
   len_bufferB = atoi(args[3]);

   bufferA_filename = args[4];
   bufferB_filename = args[5];

   output_filename = args[6];
   if (argc > 7)
      alt_buffer_filename = args[7];


   // Creating buffer arrays:
   /*
   G4double **bufferA = new G4double* [len_bufferA];
   for (int i = 0; i < len_bufferA; i ++) {
      bufferA[i] = new G4double[BUFFER_COLUMN_LEN];
   }

   G4double **alt_buffer = new G4double* [len_alt_buffer];
   for (int i = 0; i < len_bufferA; i ++) {
      alt_buffer[i] = new G4double[BUFFER_COLUMN_LEN];
   }

   G4double **bufferB = new G4double* [len_bufferB];
   for (int i = 0; i < len_bufferB; i ++) {
      bufferB[i] = new G4double[BUFFER_COLUMN_LEN];
   }
   */
   ///end creation of buffer arrays

   // Begin input of buffer data from files:

   //vector< vector<G4double> > *bufferA, *bufferB, *alt_buffer;

   vector< vector<G4double> > bufferA = vector< vector<G4double> >();
   vector< vector<G4double> > bufferB = vector< vector<G4double> >();
   vector< vector<G4double> > alt_buffer = vector< vector<G4double> >();


   vector<G4double> *v;
   G4double d;

   ifstream bufferA_file(bufferA_filename, ios::binary | ios::in);

   for (int i = 0; i < len_bufferA; i ++) {
      v =  new vector<G4double>(BUFFER_COLUMN_LEN);
      bufferA.push_back( *v );
      delete v;
      for ( int j = 0; j < BUFFER_COLUMN_LEN; j ++) {
         bufferA_file.read( reinterpret_cast<char*>( &( d ) ), sizeof(G4double) );
         bufferA.at(i).at(j) = d;
      }
   }
   bufferA_file.close();

   if ( len_alt_buffer > 0 ) {
      ifstream alt_buffer_file(alt_buffer_filename, ios::binary | ios::in);

      for (int i = 0; i < len_alt_buffer; i ++) {

         v = new vector<G4double>(BUFFER_COLUMN_LEN);
         alt_buffer.push_back( *v );
         delete v;
         for ( int j = 0; j < BUFFER_COLUMN_LEN; j ++) {
            alt_buffer_file.read( reinterpret_cast<char*>( &( d ) ), sizeof(G4double) );
            alt_buffer.at(i).at(j) = d;
         }
      }
      alt_buffer_file.close();
   }


   ifstream bufferB_file(bufferB_filename, ios::binary | ios::in);

   for (int i = 0; i < len_bufferB; i ++) {
      v =  new vector<G4double>(BUFFER_COLUMN_LEN);
      bufferB . push_back( *v );
      delete v;
      for ( int j = 0; j < BUFFER_COLUMN_LEN; j ++) {
         bufferB_file.read( reinterpret_cast<char*>( &( d ) ), sizeof(G4double) );
         bufferB . at(i).at(j) = d;
      }
   }
   bufferB_file.close();

   //////////////////////////////////////////////////////////////
   // Input is done
   //////////////////////////////////////////////////////////////

   //ErrorComputer2 *mErrorComputer = new ErrorComputer2( bufferA,
   //                                                     alt_buffer,
   //                                                     bufferB );


   Interpolant * minterpolant = new Interpolant();

   G4double F0[3];
   G4double F1[3];

   G4double t, h, xi;

   G4double left_endpoint_data[6];
   G4double right_endpoint_data[6];

   G4double interpolant[6];

   G4bool was_located_in_a_regular_interval, was_located_in_an_overshoot_interval;

   G4int len_er;
   //vector< vector<G4double> > err = *( err_ptr );

   vector< vector<G4double> > *er_ptr = new vector< vector<G4double> >();
   //vector< vector<G4double> > er = *er_ptr;


   for (int i = 0; i < bufferB.size(); i ++) {

      was_located_in_a_regular_interval = false;
      was_located_in_an_overshoot_interval = false;

      t = bufferB[i][TIME_SLOT];

      if ( t > bufferA.back().at(ENDPOINT_BASE_INDEX + TIME_SLOT) ) {
         if ( alt_buffer.size() > 0 ) {
            if ( t > alt_buffer.back().at(ENDPOINT_BASE_INDEX + TIME_SLOT) ) {
               // Cannot interpolate since we are outside of the appropriate range.
               len_er = er_ptr -> size();
               //return er_ptr; // the size of err.
               break;
            }
         }
         else {
            len_er = er_ptr -> size();
            //return er_ptr; // the size of err.
            break;
         }
      }

      v =  new vector<G4double>(BUFFER_COLUMN_LEN);
      er_ptr -> push_back( *v );
      //er.push_back( vector<G4double>(2 + NUMBER_INTERPOLATION_VARIABLES) ); // 2 for time & arclength
      delete v;

      for (int j = 0; j < bufferA.size(); j ++) {

         if ( bufferA[j][TIME_SLOT] <= t and t <= bufferA[j][ENDPOINT_BASE_INDEX + TIME_SLOT] ) {

            h = bufferA[j][ENDPOINT_BASE_INDEX + TIME_SLOT] - bufferA[j][TIME_SLOT];

            for (int k = 0; k < 3; k ++) {
               F0[k] = bufferA[j][RHS_SLOT + k];
               F1[k] = bufferA[j][ENDPOINT_BASE_INDEX + RHS_SLOT + k];
            }
            for (int k = 0; k < 6; k ++) {
               left_endpoint_data[k] = bufferA[j][POSITION_SLOT + k];
               right_endpoint_data[k] = bufferA[j][ENDPOINT_BASE_INDEX + POSITION_SLOT + k];
            }
            minterpolant -> Initialize( left_endpoint_data,
                                        right_endpoint_data,
                                        F0, F1, h );

            xi = ( t - bufferA[j][TIME_SLOT] ) / h;

            assert( 0 <= xi and xi <= 1);
            minterpolant -> InterpolatePosition( xi, interpolant );
            minterpolant -> InterpolateVelocity( xi, &(interpolant[3]) );

            // Recording part:
            er_ptr -> at(i). at(0) = bufferB[i][0]; // Record time.
            er_ptr -> at(i). at(1) = bufferB[i][1]; // Record arclength.

            for (int k = 0; k < NUMBER_INTERPOLATION_VARIABLES; k ++) { // 3 because we are only computing error for position.

               er_ptr -> at(i). at(POSITION_SLOT + k) =
                     bufferB[i][POSITION_SLOT + k] - interpolant[k];
            }


            // To see what's going wrong with velocity interpolation:

            if (i < 10) {
               cout << bufferB[i][TIME_SLOT + 1] << ", j:" << j << ", bufferA:" << bufferA[j][TIME_SLOT + 1] << ", " << bufferA[j][ENDPOINT_BASE_INDEX + TIME_SLOT + 1] << endl;
               for (int k = 0; k < 6; k ++)
                  cout << setw(10) << bufferB[i][POSITION_SLOT + k] << ", ";
               cout << endl;
               for (int k = 0; k < 6; k ++)
                     cout << setw(10) << interpolant[k] << ", ";
                  cout << endl;
               for (int k = 0; k < 6; k ++)
                  cout << setw(10) << bufferA[j][POSITION_SLOT + k] << ", ";
               cout << endl;
               for (int k = 0; k < 6; k ++)
                  cout << setw(10) << bufferA[j][ENDPOINT_BASE_INDEX + POSITION_SLOT + k] << ", ";
               cout << endl;
            }



            was_located_in_a_regular_interval = true;
            break;
         }
      }

      if ( was_located_in_a_regular_interval == true )
         continue;

      if ( alt_buffer.size() == 0 )
         continue;

      for (int j = 0; j < alt_buffer.size(); j ++) {

         if ( alt_buffer[j][TIME_SLOT] <= t and t <= alt_buffer[j][ENDPOINT_BASE_INDEX + TIME_SLOT] ) {

            h = alt_buffer[j][ENDPOINT_BASE_INDEX + TIME_SLOT] - alt_buffer[j][TIME_SLOT];

            for (int k = 0; k < 3; k ++) {
               F0[k] = alt_buffer[j][RHS_SLOT + k];
               F1[k] = alt_buffer[j][ENDPOINT_BASE_INDEX + RHS_SLOT + k];
            }
            for (int k = 0; k < 6; k ++) {
               left_endpoint_data[k] = alt_buffer[j][POSITION_SLOT + k];
               right_endpoint_data[k] = alt_buffer[j][ENDPOINT_BASE_INDEX + POSITION_SLOT + k];
            }
            minterpolant -> Initialize( left_endpoint_data,
                                        right_endpoint_data,
                                        F0, F1, h );

            xi = ( t - alt_buffer[j][TIME_SLOT] ) / h;

            assert( 0 <= xi and xi <= 1);

            minterpolant -> InterpolatePosition( xi, interpolant );
            minterpolant -> InterpolateVelocity( xi, &(interpolant[3]) );


            // Recording part:
            er_ptr -> at(i). at(0) = bufferB[i][0]; // Record time.
            er_ptr -> at(i). at(1) = bufferB[i][1]; // Record arclength.

            for (int k = 0; k < NUMBER_INTERPOLATION_VARIABLES; k ++) { // 3 because we are only computing error for position.

               er_ptr -> at(i). at(POSITION_SLOT + k) =
                     bufferB[i][POSITION_SLOT + k] - interpolant[k];
            }
            was_located_in_an_overshoot_interval =  true;
            break;
         }
      }
      if ( was_located_in_an_overshoot_interval )
         continue;

      cout << "Shouldn't be here. Wasn't in any interval??" << endl;
      assert( false );

   }

   ofstream meta_outfile("interpolation_error_metafile", ios::out);
   meta_outfile << er_ptr -> size() << endl;
   meta_outfile.close();

   ofstream output(output_filename, ios::binary | ios::out);
   for (int i = 0; i < er_ptr -> size(); i ++) {
      for (int j = 0; j < 2 + NUMBER_INTERPOLATION_VARIABLES; j ++) {
         d = er_ptr -> at(i). at(j);
         output.write( reinterpret_cast<char*>( &d ), sizeof(G4double) );
      }
   }
   output.close();

   delete er_ptr;

   return 0;
}







