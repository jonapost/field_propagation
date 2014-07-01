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
// $Id: NTSTTabulatedField3d.cc,v 1.2 2007-10-26 09:51:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  NTSTTabulatedField3D is derived from PurgingMagTabulatedField3D
//
//  Field Map of Ideal Solenoid with parameters of Babar Solenoid is used.
//  Ideal Solenoid is fully symmetric: needs 2D simulation:RZplane,no Phi
//  In Z direction only Z<0  is used, positive Z are reconstructed by symetry
//  Still have 3D aproach for possibility of use non-symmetric Solenoid.
//

#include "NTSTTabulatedField3d.hh"
#include "G4SystemOfUnits.hh"

NTSTTabulatedField3d::NTSTTabulatedField3d( const char* filename, G4double zOffset ) 
  :fZoffset(zOffset),invertR(false),invertZ(false)
{    
  count=0; 
  G4double lenUnit= meter;
  G4double fieldUnit= tesla; 
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Magnetic field"
	 << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << G4endl; 
  std::ifstream file( filename ); // Open the file for reading.
  
  // Ignore first blank line
  char buffer[256];
  file.getline(buffer,256);
  
  // Read table dimensions 
  file  >> nz >> nr; // Note dodgy order
  
  G4cout << "  [ Number of values y,z: " 
	 << " " << nr << " " << nz << " ] "
	 << G4endl;

  // Set up storage space for table
    
    rField.resize(nr);
    zField.resize(nr);
    // G4cout<<"Ok for resize nr"<<G4endl;
    G4int ir,iz;
    for (ir=0; ir<nr; ir++) {
      rField[ir].resize(nz);
      zField[ir].resize(nz);
      //  G4cout<<"Ok for resize nz"<<G4endl;
    }
  
 
  // Read in the data
  G4double xval,yval,zval,bx,by,bz,rval,br;
   
    for (ir=0; ir<nr; ir++) {
      for (iz=0; iz<nz; iz++) {
        file >> xval >> yval >> zval >> bx >> by >> bz ;
	// Ideal Solenoid is fully symmetric: needs 2D simulation:RZplane,no Phi
        // In Z direction only Z<0  is used, positive Z are reconstructed by symetry
        // Still have 3D aproach for possibility of use non-symmetric Solenoid
	// -z is used in data file just to remember that only half Z is needed for 
        // reconstruction of Ideal Solenoid Field
        zval=-zval;
        rval=std::sqrt(xval*xval+yval*yval);
        if(rval==0.){br=0.;}
        else{ br=(bx*xval+by*yval)/rval;};        
	// G4cout<<"iz="<<iz<<" ir="<<ir<<" z="<<zval<<" bz="<<bz<<G4endl;
        if (  ir==0 && iz==0 ) {
          
          minr = rval * lenUnit;
          minz = zval * lenUnit;
        }
        rField[ir][iz] = br * fieldUnit;
        zField[ir][iz] = bz * fieldUnit;
      }
    }
   
  file.close();
  
  maxr = yval * lenUnit;
  maxz = zval * lenUnit;

  G4cout << "\n ---> ... done reading " << G4endl;

  // G4cout << " Read values of field from file " << filename << G4endl; 
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values r,z: " 
	 <<  " " << minr/cm << " " << minz/cm << " cm "
	 << "\n ---> Max values r,z: " 
	 <<  " " << maxr/cm << " " << maxz/cm << " cm "
	 << "\n ---> The field will be offset by " << zOffset/cm << " cm " << G4endl;

  // Should really check that the limits are not the wrong way around.
  
  if (maxr < minr) {std::swap(maxr,minr); invertR = true;} 
  if (maxz < minz) {std::swap(maxz,minz); invertZ = true;} 
  G4cout << "\nAfter reordering if neccesary"  
	 << "\n ---> Min values r,z: " 
	 <<  " " << minr/cm << " " << minz/cm << " cm "
	 << " \n ---> Max values r,z: " 
	 <<  " " << maxr/cm << " " << maxz/cm << " cm ";

  dr = maxr - minr;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values r,z (range): " 
	 <<  " " << dr/cm << " " << dz/cm << " cm in z "
	 << "\n-----------------------------------------------------------" << G4endl;
}

void NTSTTabulatedField3d::GetFieldValue(const G4double point[4],
				      G4double *Bfield ) const
{
  count++;
  G4double x = point[0];
  G4double y = point[1];
  G4double z = point[2] + fZoffset;
  G4double absZ=std::fabs(z);
  G4double r =std::sqrt(x*x+y*y);
  G4double phi= std::atan2(x,y);
 
  // Check that the point is within the defined region 
  if ( 
       r>=minr && r<=maxr && 
       absZ>=minz && absZ<=maxz ) {
    
    // Position of given point within region, normalized to the range
    // [0,1]
    
    G4double rfraction = (r - minr) / dr; 
    G4double zfraction = (absZ - minz) / dz;

   
    if (invertR) { rfraction = 1 - rfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    G4double  rdindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    
    G4double rlocal = ( std::modf(rfraction*(nr-1), &rdindex));
    G4double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));

    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point

    G4int rindex = static_cast<G4int>(rdindex);
    G4int zindex = static_cast<G4int>(zdindex);
   
    
#ifdef DEBUG_INTERPOLATING_FIELD
    G4cout << "Local x,y,z: " <<  " " << ylocal << " " << zlocal << G4endl;
    G4cout << "Index x,y,z: " <<  " " << yindex << " " << zindex << G4endl;
    G4double valx0z0, mulx0z0, valx1z0, mulx1z0;
    G4double valx0z1, mulx0z1, valx1z1, mulx1z1;
    valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
    valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
    valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
    valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

        // Full 3-dimensional version
    //Bfield[0] =
    //  xField[yindex  ][zindex  ] *  (1-ylocal) * (1-zlocal) +
    //  xField[yindex  ][zindex+1] *  (1-ylocal) *    zlocal  +
    //  xField[yindex+1][zindex  ] *    ylocal  * (1-zlocal) +
    //  xField[yindex+1][zindex+1] *  ylocal  *    zlocal  ;

    // 2-Dimensional Table for Ideal Solenoid Field   
    Bfield[1] =
      rField[rindex  ][zindex  ] *  (1-rlocal) * (1-zlocal) +
      rField[rindex  ][zindex+1] *  (1-rlocal) *    zlocal  +
      rField[rindex+1][zindex  ] *   rlocal  * (1-zlocal) +
      rField[rindex+1][zindex+1] *     rlocal  *    zlocal  ;
      
    Bfield[2] =
      zField[rindex  ][zindex  ] * (1-rlocal) * (1-zlocal) +
      zField[rindex  ][zindex+1] * (1-rlocal) *    zlocal  +
      zField[rindex+1][zindex  ] *  rlocal  * (1-zlocal) +
      zField[rindex+1][zindex+1] *  rlocal  *    zlocal  ;
    //reconstruction of 3-dimensional field
    if(z<0.0){Bfield[1]=-1.*Bfield[1];};//correction for direction,table is only for positive z
    Bfield[0]=Bfield[1]*std::cos(phi);
    Bfield[1]=Bfield[1]*std::sin(phi);

   
  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;//0.07;//still have small Solenoid Field
  }
}

