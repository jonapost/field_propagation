//
//  Approximation of CMS magnetic field using Bi-linear approximation
//   
//  Implementation by Ananya - Nov/Dec 2015
//    Obtained from vecgeom repository, commit 15fad33
//
//  Adapted for Geant4 by J. Apostolakis, Feb 2016
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>

#include "CMSMagneticFieldG4.hh"

using namespace std;

CMSMagneticFieldG4::CMSMagneticFieldG4()
   : fInitialised(false)
{
}

CMSMagneticFieldG4::
CMSMagneticFieldG4(std::string inputMap)
{
   ReadVectorData(inputMap);   
}


CMSMagneticFieldG4::~CMSMagneticFieldG4()
{
}

void CMSMagneticFieldG4::ReadVectorData(string inputMap)
{
    string line;
    string s1,s2,s3,s4,s5,s0;
    double d1,d2,d3,d4,d5,d0;
    int ind =0;

    ifstream pFile(inputMap);
    if (pFile.is_open()){
        // getline() returns the stream. testing the stream with while returns error such as EOF
        while(getline(pFile,line)){
             // so here we know that the read was a success and that line has valid data
            stringstream ssline (line);
            //parsing all the parts. s0's store the string names which are of no use to us. 
            ssline >> s0 >> d1 >> s1 >> d0 >> s2 >> d2 >> s3 >> d3 >> s4 >> d4 >> s5 >> d5;
            fRadius.push_back(d1);
            fPhi.push_back(d0);
            fZ.push_back(d2);
            d3 *= kAInverse * CLHEP::tesla;
            d4 *= kAInverse * CLHEP::tesla;
            d5 *= kAInverse * CLHEP::tesla;            
            fBz.push_back(d3);
            fBr.push_back(d4);
            fBphi.push_back(d5);
            ind++;
            /***
            if( ind % 10 == 0 ) std::cout << "  Read in line " << ind
                                          << " B-field values= " << d3 << " " << d4 << " "
                                          << d5 << std::endl;
            ***/
        }
        pFile.close();
        const int linesExpected = kNoZValues * kNoRValues;
        G4bool  finished = ( ind == linesExpected );
        if( ! finished ) {
           std::cout<<"CMSMagneticField::ReadVectorData> "
                    << "Did not find all data for CMS Magnetic Field."
                    << std::endl
                    << "     Read  " << ind << " lines vs " << linesExpected << " expected. "
                    << std::endl;
        }
        fInitialised = finished;
    }
    else{
       std::cout<<"Unable to open file" << inputMap << " for CMS Magnetic Field."<< std::endl;
       fInitialised = false;        
    }
    if( !fInitialised ) {
       std::cout<<"FATAL ERROR: unable to read (all) data for CMS Magnetic Field."<< std::endl;
       exit(1);
    }
}


void CMSMagneticFieldG4::GetFieldValueRZ(const double r, const double Z, G4ThreeVector &rzField) const
{
    //Take care that radius and z for out of limit values take values at end points 
    double radius = min(r, kRMax);
    double z = max(min(Z, kZMax), -kZMax); //max(min(Z,Zmax), Zmin ) 

    //to make sense of the indices, consider any particular instance e.g. (25,-200)
    int rFloor = floor(radius*kRDiffInv);
    int rIndLow = rFloor*kNoZValues;
    int rIndHigh = rIndLow + kNoZValues;
    
    //if we use z-z0 in place of two loops for Z<0 and Z>0
    //z-z0 = [0,32000]
    //so indices 0 to 160 : total 161 indices for (z-z0)/200
    //i.e. we are saying:
    int zInd = floor((z-kZ0)*kZDiffInv);
    //need i1,i2,i3,i4 for 4 required indices
    int i1 = rIndLow + zInd;
    int i2 = i1 + 1;
    int i3 = rIndHigh + zInd;
    int i4 = i3 + 1;
    double zLow = (zInd- kHalfZValues)*kZDiff; //80 because it's the middle index in 0 to 160
    double zHigh = zLow + kZDiff;
    double radiusLow = rFloor*kRDiff;
    double radiusHigh = radiusLow + kRDiff;
    //cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<endl;

    //now write function
    double a1 = (radiusHigh - radius)*(zHigh - z); //area to be multiplied with i1
    double a2 = (radiusHigh - radius)*(z - zLow);
    double a3 = (radius - radiusLow)*(zHigh - z);
    double a4 = (radius - radiusLow)*(z- zLow);

    double BR   = (fBr[i1]  *a1 + fBr[i2]  *a2 + fBr[i3]  *a3 + fBr[i4]  *a4); // *kAInverse;
    double BZ   = (fBz[i1]  *a1 + fBz[i2]  *a2 + fBz[i3]  *a3 + fBz[i4]  *a4); // *kAInverse;
    double BPhi = (fBphi[i1]*a1 + fBphi[i2]*a2 + fBphi[i3]*a3 + fBphi[i4]*a4); // *kAInverse;

    //To make it thread safe. Because the previous predicted_B* vectors weren't threadsafe
    rzField.set( BR, BPhi, BZ ); 
}

//Sidenote: For theta =0; xyzField = rzField. 
//theta =0 corresponds to y=0
void CMSMagneticFieldG4::GetFieldValueXYZ(const G4ThreeVector &pos, G4ThreeVector &xyzField) const
{

    double cyl[2];
    CartesianToCylindrical(pos, cyl); 
    G4ThreeVector rzField;
    GetFieldValueRZ(cyl[0], cyl[1], rzField); //cyl[2] =[r,z]
    
    double sinTheta=0.0, cosTheta=1.0; //initialize as theta=0
    //To take care of r =0 case 
    if (cyl[0] != 0.0){
        double rInv = 1.0 / cyl[0]; 
        sinTheta = pos.y() * rInv;
        cosTheta = pos.x() * rInv;
    }

    CylindricalToCartesian(rzField, sinTheta, cosTheta, xyzField);
}


void CMSMagneticFieldG4::GetFieldValueTest(const G4ThreeVector &pos, G4ThreeVector &rzField)
{
    double cyl[2];
    CartesianToCylindrical(pos, cyl); 
    GetFieldValueRZ(cyl[0], cyl[1], rzField); //cyl[] =[r,z]
}

void CMSMagneticFieldG4::GetFieldValueRZ(const std::vector<double>   & radius,
                                         const std::vector<double>   & z,
                                               std::vector<G4ThreeVector> & rzField ) const
{
   for (size_t i = 0; i < radius.size(); ++i)
   {   
      GetFieldValueRZ(radius[i], z[i], rzField[i]);
   }
}

// #include "base/SOA3D.h"

/***
void CMSMagneticFieldG4::GetFieldValues(const vecgeom::SOA3D<double> &posVec,
                                            vecgeom::SOA3D<double> &fieldVec) const
{
for (int i = 0; i < posVec.size(); ++i)
    {
        //fill a vector3D with ith triplet for input to getFieldValue 
        G4ThreeVector pos(posVec.x(i), posVec.y(i) ,posVec.z(i));
        G4ThreeVector xyzField;
        GetFieldValueXYZ(pos, xyzField); // runs for 1 triplet
        //Fill SOA3D field with single field values 
        fieldVec.x(i) = xyzField.x();
        fieldVec.y(i) = xyzField.y();
        fieldVec.z(i) = xyzField.z();

    }
}
 ***/
