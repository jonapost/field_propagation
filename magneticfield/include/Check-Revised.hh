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
// $Id: $
//
// class G4ChordFinder
//
// Class description:
//
// It is a base class for ChordFinder classes
//
// History:
// - Created: D. Sorokin
// --------------------------------------------------------------------
// class G4ChordFinder implementation by Dmitry Sorokin
// Coding is based on G4ChordFinder. Change FindNextChord() and
// AdvanceChordLimited() to apply different algorithms for drivers
// with and without interpolation
// Sponsored by Google in Google Summer of Code 2016
//
// Supervision / code review: John Apostolakis
//
// This code is made available subject to the Geant4 license, a copy of
//  which is available at http://www.geant4.org/geant4/license/
//
//  History:
//   01.08.16 - Dmitry Sorokin - New design & first implementation
///////////////////////////////////////////////////////////////////////////////

#ifndef G4ChordFinder_HH
#define G4ChordFinder_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_EqRhs.hh"
#include "G4EquationOfMotion.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4VIntegrationDriver.hh"


class G4ChordFinder
{
public:

    G4ChordFinder(G4VIntegrationDriver* pIntDriver, G4int VerboseLevel = 1);

    //Constructor that creates defaults for all "children" classes.
    G4ChordFinder(G4MagneticField* magField, G4double stepMinimum = 1.0e-2, // * mm
                         G4MagIntegratorStepper* pStepper = nullptr,
                         G4int VerboseLevel = 1);

    ~G4ChordFinder();

    G4ChordFinder(const G4ChordFinder&) = delete;
    G4ChordFinder& operator = ( G4ChordFinder&) = delete;

    /*  Uses ODE solver's driver to find the endpoint that satisfies
     *  the chord criterion: that d_chord < delta_chord
     *  returns Length of Step taken.
     * */
    G4double  AdvanceChordLimited(G4FieldTrack& trackCurrent, G4double stepLen, G4double eps);

    G4double FindNextChord(const  G4FieldTrack& trackStart, G4double stepMax,
                           G4FieldTrack& trackEnd, G4double& dyErrPos,
                           G4double epsStep, G4double* pStepForAccuracy);


    G4FieldTrack ApproxCurvePointS(const G4FieldTrack&  curveAPointVelocity, const G4FieldTrack&  curveBPointVelocity,
                                   const G4FieldTrack&  ApproxCurveV, const G4ThreeVector& currentEPoint,
                                   const G4ThreeVector& currentFPoint, const G4ThreeVector& PointG,
                                   G4bool first,  G4double epsStep);

    G4FieldTrack ApproxCurvePointV(const G4FieldTrack& trackPointA, const G4FieldTrack& trackPointB,
                                   const G4ThreeVector& pointE, G4double epsStep);




    //inline function definition

    /*
     * Resets internal state of G4ChordFinder
     * this is needed for drivers with interpolation of
     * Adams methods
     * */
    inline void Reset();

    inline G4bool DoStepForIntersection(G4FieldTrack& track, G4double stepLen, G4double eps);

    inline G4EquationOfMotion* GetEquationOfMotion();
    inline void SetEquationOfMotion(G4EquationOfMotion* newEquation);

    inline void SetIntegrationDriver(G4VIntegrationDriver* IntDriver);
    inline G4VIntegrationDriver* GetIntegrationDriver();

    inline G4double InvParabolic( const G4double xa, const G4double ya,
                                  const G4double xb, const G4double yb,
                                  const G4double xc, const G4double yc );

    inline G4double  GetDeltaChord() const;

    inline void SetDeltaChord(G4double newDelta);

    // Clear internal state (last step estimate)
    inline void ResetStepEstimate();

    // Get statistics about number of calls & trials in FindNextChord
    inline G4int GetNoCalls();
    inline G4int GetNoTrials();        // Total number of trials
    inline G4int GetNoMaxTrials();     // Maximum # of trials for one call

    // Set verbosity and return old value
    inline G4int SetVerbose( G4int newvalue=1);

    inline G4int GetVerbose() const;

    // Parameter for  performance ... change with great care
    inline void SetFirstFraction(G4double newVal);
    inline G4double GetLastStepEstimateUnc();
    inline void     SetLastStepEstimateUnc( G4double stepEst );


    //   Printing for monitoring ...
    inline   G4double GetFirstFraction();         // Originally 0.999
    inline   G4double GetFractionLast();          // Originally 1.000
    inline   G4double GetFractionNextEstimate();  // Originally 0.980
    inline   G4double GetMultipleRadius();        // No original value

    inline G4bool AcceptableMissDist(G4double dChordStep) const;

    /*Accumulate the basic statistics
     *other specialised ones must be kept by derived classes
     * */
    inline void AccumulateStatistics( G4int noTrials );


    // A report with the above -- and possibly other stats
    void PrintStatistics();

    // Parameters for  performance ... change with great care
    void SetFractions_Last_Next( G4double fractLast= 0.90,
                                 G4double fractNext= 0.95 );

    void TestChordPrint(G4int noTrials,
                        G4int    lastStepTrial,
                        G4double dChordStep,
                        G4double nextStepTrial );



    G4double NewStep( G4double stepTrialOld,
                      G4double dChordStep,     // Current dchord estimate
                      G4double& stepEstimate_Unconstrained );


    void PrintDchordTrial(G4int noTrials,
                          G4double  stepTrial,
                          G4double  oldStepTrial,
                          G4double  dChordStep);



private:

     // Dependent object
     G4VIntegrationDriver* fpIntDriver;
     G4MagIntegratorStepper* fpStepper;
     G4EquationOfMotion* fpEquation;

     G4bool fAllocatedStepper;
     G4bool fAllocatedEquation;

     // Verbose level
     G4int fVerboseLevel;

     // SET in G4VRevisedChordFinderRefined.cc = 0.25 mm
     const G4double fDefaultDeltaChord;

     //  PARAMETERS
     //  ---------------------
     G4double  fDeltaChord;               //  Maximum miss distance
     //    Internal parameters
     G4double  fFirstFraction, fFractionLast, fFractionNextEstimate;
     G4double  fMultipleRadius;
     G4int     fStatsVerbose;  // if > 0, print Statistics in destructor


     //  STATE information
     G4double    fLastStepEstimate_Unconstrained;
       //  State information for efficiency

     // For Statistics
     // -- G4int   fNoTrials, fNoCalls;
     G4int   fTotalNoTrials_FNC,  fNoCalls_FNC, fmaxTrials_FNC; // fnoTimesMaxTrFNC;
};

//include inline function definition
#include "G4ChordFinder.icc"

#endif
