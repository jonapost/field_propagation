//
// class BSChordFinder
//
// Class description:
//
// BSChordFinder is a copy of G4ChordFinder with
// G4MagInt_Driver replaced by BaseDriver
//
// History:
// - Created: D. Sorokin
// --------------------------------------------------------------------


#ifndef BSChordFinder_HH
#define BSChordFinder_HH

#include "G4MagIntegratorDriver.hh"
#include "G4FieldTrack.hh"

#include "BaseDriver.hh"


class BSChordFinder
{
   public:  // with description

      BSChordFinder(BaseDriver* BSDriver);

      virtual ~BSChordFinder();

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

      inline void SetIntegrationDriver(BaseDriver* IntegrationDriver);
      inline BaseDriver* GetIntegrationDriver();
        // Access and set Driver.

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

      BSChordFinder(const BSChordFinder&);
      BSChordFinder& operator=(const BSChordFinder&);
        // Private copy constructor and assignment operator.

   private:  // ............................................................
                                          // G4int    nOK, nBAD;

      // Constants
      const G4double fDefaultDeltaChord;  // SET in BSChordFinder.cc = 0.25 mm

      //  PARAMETERS
      //  ---------------------
      G4double  fDeltaChord;               //  Maximum miss distance
      //    Internal parameters
      G4double  fFirstFraction, fFractionLast, fFractionNextEstimate;
      G4double  fMultipleRadius;
      G4int     fStatsVerbose;  // if > 0, print Statistics in destructor

      //  DEPENDENT Objects
      //  ---------------------
      BaseDriver*        fIntgrDriver;

      //  STATE information
      //  --------------------
      G4double    fLastStepEstimate_Unconstrained;
        //  State information for efficiency

      // For Statistics
      // -- G4int   fNoTrials, fNoCalls;
      G4int   fTotalNoTrials_FNC,  fNoCalls_FNC, fmaxTrials_FNC; // fnoTimesMaxTrFNC;
};

#include "BSChordFinder.icc"

#endif  // BSChordFinder_HH
