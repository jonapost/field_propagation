
#ifndef TQUADRUPOLEMAGFIELD_HH
#define TQUADRUPOLEMAGFIELD_HH

#include "G4QuadrupoleMagField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

static G4RotationMatrix IdentityMatrix; 

//template
//<class T_Field>
class TQuadrupoleMagField : public G4QuadrupoleMagField 
{
public: // with description

	TQuadrupoleMagField(G4double pGradient)
		         :G4QuadrupoleMagField(pGradient)
	{
		fGradient = pGradient ;
		fOrigin   = G4ThreeVector( 0.0, 0.0, 0.0) ;
		fpMatrix  = &IdentityMatrix;
	}

	TQuadrupoleMagField(G4double pGradient, 
						G4ThreeVector pOrigin, 
						G4RotationMatrix* pMatrix)
				:G4QuadrupoleMagField(pGradient,
									   pOrigin, 
									   pMatrix)
	{
		fGradient    = pGradient ;
		fOrigin      = pOrigin ;
		fpMatrix     = pMatrix ;
	}

	~TQuadrupoleMagField() {;}

	inline void GetFieldValue(const G4double y[7],
			                 G4double B[3]     ) const
	{
		G4ThreeVector r_global = G4ThreeVector(
				y[0] - fOrigin.x(), 
				y[1] - fOrigin.y(), 
				y[2] - fOrigin.z());

		G4ThreeVector r_local = G4ThreeVector(
				fpMatrix->colX() * r_global,
				fpMatrix->colY() * r_global,
				fpMatrix->colZ() * r_global);

		G4ThreeVector B_local = G4ThreeVector(
				fGradient * r_local.y(),
				fGradient * r_local.x(),
				0);

		G4ThreeVector B_global = G4ThreeVector(
				fpMatrix->inverse().rowX() * B_local,
				fpMatrix->inverse().rowY() * B_local,
				fpMatrix->inverse().rowZ() * B_local);

		B[0] = B_global.x() ;
		B[1] = B_global.y() ;
		B[2] = B_global.z() ;
	}

	TQuadrupoleMagField* Clone() const
	{
		//TODO: Can the fpMatrix be shared??
		return new TQuadrupoleMagField(this->fGradient,
                                        this->fOrigin,
				                        this->fpMatrix);
	}

private:

	G4double          fGradient;
	G4ThreeVector     fOrigin;
	G4RotationMatrix* fpMatrix;
};
#endif

