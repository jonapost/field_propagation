//
// class G4BSChordFinder
//
// Class description:
//
// Specialisatoin of G4VChordFinder class
// for Bulirsch-Stoer method with interpolation
// (BulirschStoerDenseOut)
//
// History:
// - Created: D. Sorokin
// --------------------------------------------------------------------

#ifndef G4BSChordFinder_HH
#define G4BSChordFinder_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_EqRhs.hh"
#include "G4EquationOfMotion.hh"
#include "G4MagIntegratorStepper.hh"

#include "BulirschStoerDenseOut.hh"
#include "G4VChordFinder.hh"

class G4BSChordFinder : public G4VChordFinder
{
public:

    G4BSChordFinder(G4double hminimum,
               G4EquationOfMotion* equation,
               G4int numberOfComponents = 6,
               G4int statisticsVerbosity = 1);

    ~G4BSChordFinder();

    G4BSChordFinder(const G4BSChordFinder&);

    G4BSChordFinder& operator = (G4BSChordFinder&);


    virtual G4double AdvanceChordLimited(G4FieldTrack& trackCurrent,
                                         G4double stepLen,
                                         G4double eps) override final;

    virtual void reset() override final;


    virtual G4bool DoStepForIntersection(G4FieldTrack&  track,
                                         G4double stepLen,
                                         G4double eps) override final;

    void SetEquationOfMotion(G4EquationOfMotion* newEquation)
    {GetIntegrationDriver()->SetEquationOfMotion(newEquation);}

    G4EquationOfMotion* GetEquationOfMotion()
    {return GetIntegrationDriver()->GetEquationOfMotion();}



private:

    void DoStepForChord(G4FieldTrack&  track, G4double stepLen);
    void DoInterpolation(G4FieldTrack& track, G4double hstep, G4double eps = 0);

    /*virtual */G4double OneGoodStep(G4FieldTrack&  track,
                                 G4double stepLen,
                                 G4double eps)/* override final*/;

    /*  Runge-Kutta driver with adaptive stepsize control. Integrate starting
     *  values at y_current over hstep x2 with accuracy eps.
     *  On output ystart is replaced by values at the end of the integration
     *  interval. RightHandSide is the right-hand side of ODE system.
     *  The source is similar to odeint routine from NRC p.721-722 .
     */
    G4bool AccurateAdvance(G4FieldTrack& trackInOut, G4double hstep,
                           G4double eps, G4double hinitial);



    /*virtual */G4double FindNextChord(G4FieldTrack& trackCurrent,
                           G4double     stepMax)/* override final*/;


    //inline function definition
    G4int GetNumberOfVariables()
    {return fIntgrDriver.GetNumberOfVariables();}


private:
     BulirschStoerDenseOut fIntgrDriver;
     G4double fMinimumStep;
     G4int fverb;

     G4double tBegin,tEnd;
     G4double yIn[G4FieldTrack::ncompSVEC],
              yOut[G4FieldTrack::ncompSVEC],
              dydx[G4FieldTrack::ncompSVEC];

     G4double nextStepSize;
     G4double lastStepSize;

     G4double eps_prev;
};

#endif
