#ifndef GUSTAFSSONDRIVER_HH
#define GUSTAFSSONDRIVER_HH

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"

class GustafssonDriver {
public:
    GustafssonDriver(G4double minimumStep, G4MagIntegratorStepper* stepper);

    void OneGoodStep(G4double y[],
                     const G4double dydx[],
                     G4double& trackLength,
                     G4double htry,
                     G4double errorTolerance,
                     G4double& hdid,
                     G4double& hnext);

    G4MagIntegratorStepper* GetStepper() {return fstepper;}
    const G4MagIntegratorStepper* GetStepper() const {return fstepper;}

    void setSafety(G4double safety) {fsafety = safety;}
    void setTrials(G4int trials) {ftrials = trials;}
    void setMinimumStep(G4double minimumStep) {fminimumStep = minimumStep;}
    void setMaxStepIncrease(G4double maxStepIncrease)
        {fmaxStepIncrease = maxStepIncrease;}

private:
    G4double shrinkStep(G4double error, G4double hstep);
    G4double growStep(G4double error, G4double hstep);

    G4MagIntegratorStepper* fstepper;

    G4int ftrials;
    G4int fmaxtrials;
    G4double fsafety;
    G4double fminimumStep;
    G4double KI;
    G4double KP;
    G4double fmaxStepIncrease;
    G4double ferrorPrev;
    G4double fpGrow;
    G4double fErrcon;
};

#endif
