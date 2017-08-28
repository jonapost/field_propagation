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
// Helper namespace 'magneticfield'
//
// Description:
//    Simple methods to extract vectors from arrays in conventions of
//     the magnetic field integration
//
//    Implementation by Dmitry Sorokin - GSoC 2017
//       Work supported by Google as part of Google Summer of Code 2017.
//    Supervision / code review: John Apostolakis

#ifndef G4FieldUtils_HH
#define G4FieldUtils_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"

namespace magneticfield {

enum class Value3D {
    Position = 0,
    Momentum = 3,
    Spin = 9
};

enum class Value1D {
    KineticEnergy = 6,
    LabTime = 7,
    ProperTime = 8
};

inline G4double extractValue(const G4double array[], const Value1D& value);
inline G4double extractValue2(const G4double array[], const Value1D& value);

inline G4double extractValue(const G4double array[], const Value3D& value);
inline G4double extractValue2(const G4double array[], const Value3D& value);

inline G4ThreeVector makeVector(const G4double array[], const Value3D& value);

G4double relativeError(
    const G4double y[],
    const G4double yerr[],
    const G4double hstep,
    const G4double errorTolerance = 1);

} //magneticfield

#include "G4FieldUtils.icc"

#endif
