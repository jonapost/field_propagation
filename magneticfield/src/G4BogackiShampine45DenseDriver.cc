#include "G4BogackiShampine45DenseDriver.hh"
#include "BogackiShampine45.hh"
#include "G4SystemOfUnits.hh"

#define ncomp G4FieldTrack::ncompSVEC

G4BogackiShampine45DenseDriver::G4BogackiShampine45DenseDriver(G4double hminimum,
                                                               G4EquationOfMotion* pEquation,
                                                               G4int numberOfComponents ,
                                                               G4int VerboseLevel):
    G4MagInt_Driver (hminimum,
                     new BogackiShampine45(pEquation, numberOfComponents),
                     numberOfComponents,
                     VerboseLevel),
   eps_prev(0),fNextStepSize(DBL_MAX)
{
}

G4BogackiShampine45DenseDriver::~G4BogackiShampine45DenseDriver()
{
}

void G4BogackiShampine45DenseDriver::DoStep(G4FieldTrack& track, G4double hstep, G4double eps)
{
    //prepare
    track.DumpToArray(yCurrent);
    GetEquationOfMotion()->RightHandSide(yCurrent, dydxCurrent);
    interpolationInterval& interval = GetInterpolationInterval();
    interval.first = interval.second = track.GetCurveLength();
    G4double stepLen = std::min(hstep, fNextStepSize);
    eps_prev = eps;
    G4double hdid = 0;

    //do one step with error control
    OneGoodStep(yCurrent, dydxCurrent, interval.second, stepLen, eps, hdid, fNextStepSize);

    //update track
    track.LoadFromArray(yCurrent, ncomp);
    track.SetCurveLength(interval.second);
}

void G4BogackiShampine45DenseDriver::DoInterpolation(G4FieldTrack& track, G4double hstep, G4double eps)
{
    G4double curveLength = track.GetCurveLength();
    G4double clWant = curveLength + hstep;
    interpolationInterval& interval = GetInterpolationInterval();

    //little upperflow, allow.

    if (clWant > interval.second)
    {
        G4double upperflow = (clWant - interval.second)/clWant;
        if (upperflow < perMillion)
            clWant = interval.second;
    }

    if (curveLength >= interval.first && clWant <= interval.second)
    {
        track.DumpToArray(yCurrent);
        const BogackiShampine45* stepper = static_cast<const BogackiShampine45*>(GetStepper());
        stepper->Interpolate(clWant - interval.first, yCurrent);
        track.LoadFromArray(yCurrent, ncomp);
        track.SetCurveLength(clWant);
        if (eps != 0 && eps != eps_prev)
        {
            char buff[256];
            sprintf(buff,"Accuracy changed. eps: %g, eps_prev: "
                         "%g Interpolation is not accurate!",eps,eps_prev);
            G4Exception("G4BogackiShampine45DenseDriver::DoInterpolation()",
                        "GeomField0001", FatalException, buff);
        }
    }
    else
    {
        char buff[256];
        sprintf(buff,"curveLength = %g is out of the interpolation interval "
                     "(%g,%g)!",clWant, interval.first, interval.second);
        G4Exception("G4BogackiShampine45DenseDriver::DoInterpolation()",
                    "GeomField0001", FatalException, buff);
    }
}

G4bool G4BogackiShampine45DenseDriver::isDense() const
{return true;}
