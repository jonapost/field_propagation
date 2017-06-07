//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4NewChordFinder.hh 97285 2016-05-31 19:52:18Z japost $
// 
// Class G4NewChordFinder
//
// class description:
//
//    Revised version of G4ChordFinder, able to use either the
//  existing G4MagIntegratorDriver or the new G4FSALIntegrationDriver.
//  Each object may have either of these drivers (not both.)
//
//  Purpose: first trial integration of FSAL drivers at ChordFinder
//             level.
//  To Do:   check performance, look into replacing with
//              template < class Driver > G4TChordFinder
//           which could use either type of driver without 'if's
//
// History:
// - 31.05.16 - John Apostolakis - Design and implementation 
// -------------------------------------------------------------------

#ifndef G4NEWCHORDFINDER_HH
#define G4NEWCHORDFINDER_HH

#include "G4MagIntegratorDriver.hh"
#include "G4FSALIntegrationDriver.hh"
class G4VFSALIntegrationStepper;
// #include "G4VFSALIntegrationStepper.hh"

#include "G4FieldTrack.hh"
class G4MagneticField;  


class G4NewChordFinder
{ 
   public:  // with description

      G4NewChordFinder( G4MagInt_Driver* pIntegrationDriver );

      G4NewChordFinder( G4MagneticField* itsMagField,
                        G4double         stepMinimum = 1.0e-2, // * mm 
                        G4MagIntegratorStepper* pItsStepper = 0 );  
        // A constructor that creates defaults for all "children" classes.

      G4NewChordFinder( G4FSALIntegrationDriver* pFSALIntegrationDriver);
        // Experimental for FSAL drivers and steppers

      G4NewChordFinder( G4MagneticField*           itsMagField,
                        G4VFSALIntegrationStepper* pFSALStepper,
                        G4double                   stepMinimum = 1.0e-2 // * mm 
         );  
        // A constructor that creates defaults for all "children" classes.
   
      virtual ~G4NewChordFinder();

      G4double    AdvanceChordLimited( G4FieldTrack& yCurrent,
                                       G4double stepInitial,
                                       G4double epsStep_Relative,
                                       const G4ThreeVector latestSafetyOrigin,
                                       G4double lasestSafetyRadius);
        // Uses ODE solver's driver to find the endpoint that satisfies 
        // the chord criterion: that d_chord < delta_chord
        // -> Returns Length of Step taken.
     
      G4FieldTrack ApproxCurvePointS( const G4FieldTrack&  curveAPointVelocity,
                                      const G4FieldTrack&  curveBPointVelocity,
                                      const G4FieldTrack&  ApproxCurveV,
                                      const G4ThreeVector& currentEPoint,
                                      const G4ThreeVector& currentFPoint,
                                      const G4ThreeVector& PointG,
                                            G4bool first,  G4double epsStep);
 
      G4FieldTrack ApproxCurvePointV( const G4FieldTrack&  curveAPointVelocity,
                                      const G4FieldTrack&  curveBPointVelocity,
                                      const G4ThreeVector& currentEPoint,
                                            G4double       epsStep);

      inline G4double InvParabolic( const G4double xa, const G4double ya,
                                    const G4double xb, const G4double yb,
                                    const G4double xc, const G4double yc );

      inline G4double  GetDeltaChord() const;
      inline void      SetDeltaChord(G4double newval);

      inline void SetIntegrationDriver(G4MagInt_Driver* IntegrationDriver);
      inline G4MagInt_Driver* GetIntegrationDriver();
        // Access and set Driver.
  
      inline void SetFSALIntegrationDriver(G4FSALIntegrationDriver* fsalDriver);
      inline G4FSALIntegrationDriver* GetFSALIntegrationDriver();
         // Access and set Driver.
  
      inline G4EquationOfMotion* GetEquationOfMotion();
        //  New - to cope with different type for FSAL steppers
  
      inline void ResetStepEstimate();
        // Clear internal state (last step estimate)

      inline G4int GetNoCalls(); 
      inline G4int GetNoTrials();        // Total number of trials
      inline G4int GetNoMaxTrials();     // Maximum # of trials for one call
        // Get statistics about number of calls & trials in FindNextChord

      virtual void   PrintStatistics(); 
        // A report with the above -- and possibly other stats
      inline G4int SetVerbose( G4int newvalue=1); 
        // Set verbosity and return old value

      void SetFractions_Last_Next( G4double fractLast= 0.90, 
                                   G4double fractNext= 0.95 ); 
        // Parameters for  performance ... change with great care

      inline void SetFirstFraction(G4double fractFirst);
        // Parameter for  performance ... change with great care

   public:  // without description

      void     TestChordPrint( G4int    noTrials, 
                               G4int    lastStepTrial, 
                               G4double dChordStep, 
                               G4double nextStepTrial );

        //   Printing for monitoring ...
 
      inline   G4double GetFirstFraction();         // Originally 0.999
      inline   G4double GetFractionLast();          // Originally 1.000
      inline   G4double GetFractionNextEstimate();  // Originally 0.980
      inline   G4double GetMultipleRadius();        // No original value
        //  Parameters for adapting performance ... use with great care

   private: 
      G4NewChordFinder();   
        //  Default constructor - only for use in implementer 
        //      other constructors
   
   protected:   // .........................................................

   
      inline  void    AccumulateStatistics( G4int noTrials ); 
        // Accumulate the basic statistics 
        //   - other specialised ones must be kept by derived classes
 
      inline G4bool AcceptableMissDist(G4double dChordStep) const;

      G4double NewStep( G4double stepTrialOld, 
                        G4double dChordStep,     // Current dchord estimate
                        G4double& stepEstimate_Unconstrained ) ;  
      
      virtual G4double FindNextChord( const  G4FieldTrack& yStart,
                              G4double     stepMax,
                              G4FieldTrack& yEnd,
                              G4double&    dyErr,      //  Error of endpoint 
                              G4double     epsStep,
                              G4double*  pNextStepForAccuracy,  // = 0,
                              const G4ThreeVector latestSafetyOrigin,
                              G4double       latestSafetyRadius 
                                      );  

      void     PrintDchordTrial(G4int     noTrials, 
                                G4double  stepTrial, 
                                G4double  oldStepTrial, 
                                G4double  dChordStep);

      inline G4double GetLastStepEstimateUnc(); 
      inline void     SetLastStepEstimateUnc( G4double stepEst ); 

   private:  // ............................................................

      G4NewChordFinder(const G4NewChordFinder&);
      G4NewChordFinder& operator=(const G4NewChordFinder&);
        // Private copy constructor and assignment operator.

   private:  // ............................................................
                                          // G4int    nOK, nBAD;

      // Constants
      const G4double fDefaultDeltaChord;  // SET in G4NewChordFinder.cc = 0.25 mm

      //  PARAMETERS 
      //  ---------------------
      G4double  fDeltaChord;               //  Maximum miss distance 
      //    Internal parameters
      G4double  fFirstFraction, fFractionLast, fFractionNextEstimate;
      G4double  fMultipleRadius; 
      G4int     fStatsVerbose;  // if > 0, print Statistics in destructor

      //  DEPENDENT Objects
      //  ---------------------
      G4MagInt_Driver*        fIntgrDriver;
      G4MagIntegratorStepper* fDriversStepper;

      G4FSALIntegrationDriver*     fFSALIntgrDriver;
      G4VFSALIntegrationStepper* fFSALStepper;
   
      G4bool                  fAllocatedStepper;  // Bookkeeping of dependent object
      G4EquationOfMotion*     fEquation; 

      //  STATE information
      //  --------------------
      G4double    fLastStepEstimate_Unconstrained;
        //  State information for efficiency

      // For Statistics
      // -- G4int   fNoTrials, fNoCalls;
      G4int   fTotalNoTrials_FNC,  fNoCalls_FNC, fmaxTrials_FNC; // fnoTimesMaxTrFNC; 
};

// Inline function implementation:

#include "G4NewChordFinder.icc"

#endif  // G4NEWCHORDFINDER_HH
