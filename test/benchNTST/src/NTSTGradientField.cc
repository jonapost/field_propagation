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
// $Id: NTSTGradientField.cc,v 1.2 2007-10-26 09:51:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  NTSTGradientField is derived from G4QuadripoleMagField:
//
//  Bx=fGradient*x; By=-fGradient*y; Bz=1.5*tesla(Babar Solenoid);
//
// -------------------------------------------------------------------

#include "NTSTGradientField.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

NTSTGradientField::NTSTGradientField(G4double pGradient)
{
   fGradient = pGradient ;
   count=0;
}

/////////////////////////////////////////////////////////////////////////

NTSTGradientField::~NTSTGradientField()
{
}

////////////////////////////////////////////////////////////////////////

void NTSTGradientField::GetFieldValue( const G4double y[7],
                                                G4double B[3]  ) const  
{
   //   G4double fGradient = 0.001 ;   // Tesla/mm
   count++;
   B[0] = fGradient*y[1] ;
   B[1] =- fGradient*y[0] ;
  
   B[2]=1.5*tesla;//Babar Solenoid
   // G4cout<<"B="<<B[0]/tesla<<"  "<<B[1]/tesla<<"  "<<B[2]/tesla<<G4endl;
}
