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
// $Id: $
//
// Helper namespace 'magneticfield', 'internal'
//
// Description:
// class G4BulirschStoerDriver implementation by Dmitry Sorokin
// This driver class uses Bulirsch-Stoer method to integrate
// the equation of motion
//
//    Implementation by Dmitry Sorokin - GSoC 2016
//       Work supported by Google as part of Google Summer of Code 2016.
//    Supervision / code review: John Apostolakis
//
///////////////////////////////////////////////////////////////////////////////

#ifndef G4BULIRSCH_STOER_DRIVER_HH
#define G4BULIRSCH_STOER_DRIVER_HH

#include "G4VIntegrationDriver.hh"
#include "G4BulirschStoer.hh"

class G4BulirschStoerDriver: public G4VIntegrationDriver {
public:
    G4BulirschStoerDriver(
        G4double hminimum,
        G4EquationOfMotion* equation,
        G4int integratedComponents = 6,
        G4int verb = 1);

    ~G4BulirschStoerDriver() = default;

    G4BulirschStoerDriver(const G4BulirschStoerDriver&) = delete;
    G4BulirschStoerDriver& operator=(const G4BulirschStoerDriver&) = delete;

    virtual G4bool AccurateAdvance(
        G4FieldTrack& track,
        G4double stepLen,
        G4double eps,
        G4double beginStep = 0) override final;

    virtual G4bool QuickAdvance(
        G4FieldTrack& y_val,
        const G4double dydx[],
        G4double hstep,
        G4double& missDist,
        G4double& dyerr) override final;

    void OneGoodStep(
        G4double y[],
        const G4double dydx[],
        G4double& curveLength,
        G4double htry,
        G4double eps,
        G4double& hdid,
        G4double& hnext);

    virtual void GetDerivatives(
        const G4FieldTrack& track,
        G4double dydx[]) const override final;

    virtual void SetVerboseLevel(G4int level) override final;
    virtual G4int GetVerboseLevel() const override final;

    virtual G4double ComputeNewStepSize(
        G4double  errMaxNorm,    // normalised error
        G4double  hstepCurrent) override final; // current step size

    virtual G4EquationOfMotion* GetEquationOfMotion() override final;
    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override final;

private:
    G4EquationOfMotion* fEquation;
    G4double fMinimumStep;
    G4int fNumberOfComponents;
    G4double fVerbosity;

    internal::G4ModifiedMidpoint fMidpointMethod;
    internal::G4BulirschStoer bulirschStoer;

    G4double yIn[G4FieldTrack::ncompSVEC],
             yMid[G4FieldTrack::ncompSVEC],
             yMid2[G4FieldTrack::ncompSVEC],
             yOut[G4FieldTrack::ncompSVEC],
             yOut2[G4FieldTrack::ncompSVEC],
             yError[G4FieldTrack::ncompSVEC];


    G4double dydxCurrent[G4FieldTrack::ncompSVEC];
    G4double yCurrent[G4FieldTrack::ncompSVEC];

    G4double derivs[2][6][G4FieldTrack::ncompSVEC];

    // G4double diffs[4][2][G4FieldTrack::ncompSVEC];

    const G4int interval_sequence[2];
};


#endif
