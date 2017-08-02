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
// -------------------------------------------------------------------

#include "ARK4.hh"

#include "G4FieldTrack.hh"
#include "BettisFunctions.hh"
#include "Utils.hh"

ARK4::ARK4(
    G4EquationOfMotion* equation,
    G4int numberOfVariables)
    : G4MagErrorStepper(equation, numberOfVariables)
{
}

void ARK4::DumbStepper(
    const G4double yIn[],
    const G4double dydx[],
    G4double h,
    G4double yOut[])
{
    G4double hh = h * 0.5;

    G4double ytemp[G4FieldTrack::ncompSVEC];
    G4double ak2[G4FieldTrack::ncompSVEC],
             ak3[G4FieldTrack::ncompSVEC],
             ak4[G4FieldTrack::ncompSVEC];

    // Initialise time to t0, needed when it is not updated by the integration.
    //        [ Note: Only for time dependent fields (usually electric)
    //                  is it neccessary to integrate the time.]
    ytemp[7] = yIn[7];
    yOut[7] = yIn[7];

    //RightHandSide(yIn, dydx);
    for (G4int i = 0; i < GetNumberOfVariables(); ++i) {
        ytemp[i] = yIn[i] + hh * dydx[i];
    }

    RightHandSide(ytemp, ak2);
    for (G4int i = 0; i < GetNumberOfVariables(); ++i) {
        ytemp[i] = yIn[i] + hh * ak2[i];
    }

    RightHandSide(ytemp, ak3);
    for (G4int i = 0; i < GetNumberOfVariables(); ++i) {
        ytemp[i] = yIn[i] + h * ak3[i];
    }

    // TODO use bettis function!
    //double nu =
    G4double b1 = 1. / 6.,
             b2 = 2. / 6.,
             b3 = 2. / 6.,
             b4 = 1. / 6.;

    RightHandSide(ytemp, ak4);
    for (G4int i = 0; i < GetNumberOfVariables(); ++i) {
        yOut[i] = yIn[i] + h * (b1 * dydx[i] + b2 * ak2[i] + b3 * ak3[i] + b4 * ak4[i]);
    }
}

