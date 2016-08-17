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
// $Id: NTSTTabulatedField3d.hh,v 1.2 2007-10-26 09:51:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  Class NTSTTabulatedField3d.hh derived from class 
//        PurgMagTabulatedField3D from examples/advanced/PurgMagnet
//  Class NTSTTabulatedField3d.hh reads the FieldTable and
//  gives the value of field in any x,y,z point.
//  FieldTable represents Field Map of Fully symmetric Solenoid with
//  dimentions of Babar Solenoid
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTTabulatedField3d_h
#define NTSTTabulatedField3d_h 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>

class NTSTTabulatedField3d
#ifndef STANDALONE
 : public G4MagneticField
#endif
{
  
  // Storage space for the table
  
  std::vector< std::vector<  G4double  > > rField;
  std::vector< std::vector<  G4double  > > zField;
  // The dimensions of the table
  G4int nr,nz; 
  // The physical limits of the defined region
  G4double  minr, maxr, minz, maxz;
  // The physical extent of the defined region
  G4double  dr, dz;
  G4double fZoffset;//not used in this test
  G4bool  invertR, invertZ;

public:

  NTSTTabulatedField3d(const char* filename, G4double zOffset );
  void  GetFieldValue( const  G4double Point[4],
		       G4double *Bfield          ) const;

  //without description

  G4int GetCount() const {return count;}
  void ClearCount() {count = 0;}

  protected:
  mutable G4int count;

};

#endif
