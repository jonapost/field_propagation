#ifndef G4INTEGRATION_OBSERVER_HH
#define G4INTEGRATION_OBSERVER_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"

class G4IntegrationObserver {
public:
    G4IntegrationObserver(
            G4int verboseLevel,
            G4double minimumStepLength,
            G4double numberOfIntegrationVariables);
    ~G4IntegrationObserver();

    void WarnSmallStepSize(G4double hnext, G4double hstep,
                           G4double h, G4double xDone,
                           G4int stepNumber);

    void WarnTooManyStep(G4double x1start, G4double x2end, G4double xCurrent);

    void WarnEndPointTooFar(G4double endPointDist,
                            G4double hStepSize,
                            G4double epsilonRelative,
                            G4int debugFlag);
    //  Issue warnings for undesirable situations

    void PrintStatus(const G4double* StartArr,
                    G4double xstart,
                    const G4double* CurrentArr,
                    G4double xcurrent,
                    G4double requestStep,
                    G4int subStepNo);

    void PrintStatus(const G4FieldTrack& StartFT,
                     const G4FieldTrack& CurrentFT,
                     G4double requestStep,
                     G4int subStepNo);

    void PrintStat_Aux(const G4FieldTrack& aFieldTrack,
                       G4double requestStep,
                       G4double actualStep,
                       G4int subStepNo,
                       G4double subStepSize,
                       G4double dotVelocities);
     //  Verbose output for debugging

    //void PrintStatisticsReport();
    //  Report on the number of steps, maximum errors etc.

    void onRelativeError(
        G4double H_lg,
        G4double DyerrPos_lg,
        G4double DyerrVel_lg);
private:
    G4int fVerboseLevel;
    G4double fMinimumStep;
    G4double fNoIntegrationVariables;

    G4double fNoTotalSteps;
    G4double fNoBadSteps;
    G4double fNoSmallSteps;

    G4double fSumH_lg;
    G4double fDyerrPos_lgTot;
    G4double fDyerrVel_lgTot;
};

#endif
