
#ifndef TCACHED_MAGNETIC_FIELD_DEF
#define TCACHED_MAGNETIC_FIELD_DEF

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"

template
<class T_Field>
class TCachedMagneticField : public G4MagneticField
{
public:  // with description

TCachedMagneticField(T_Field* pTField, G4double distance)
	    : fLastLocation(DBL_MAX,DBL_MAX,DBL_MAX),
	      fLastValue(DBL_MAX,DBL_MAX,DBL_MAX),
	      fCountCalls(0),  fCountEvaluations(0)
{
	fpMagneticField= pTField;
	fDistanceConst= distance;

	// G4cout << " Cached-B-Field constructor> Distance = " << distance << G4endl;
	this->ClearCounts();
}

TCachedMagneticField(const TCachedMagneticField<T_Field> &rightCMF)
{
	fpMagneticField=      rightCMF.fpMagneticField;
	fDistanceConst = rightCMF.fDistanceConst;
	fLastLocation  = rightCMF.fLastLocation;
	fLastValue     = rightCMF.fLastValue;
	this->ClearCounts(); 
}

TCachedMagneticField* Clone() const
{
	G4cout << "Clone is called" << G4endl;
	//Cannot use copy constructor: I need to clone the associated magnetif field
	T_Field* aF = this->fpMagneticField->T_Field::Clone();
	TCachedMagneticField* cloned = new TCachedMagneticField( aF ,
			this->fDistanceConst );
	cloned->fLastLocation = this->fLastLocation;
	cloned->fLastValue = this->fLastValue;
	return cloned;
}

virtual ~TCachedMagneticField() {;}
// Constructor and destructor. No actions.

void ReportStatistics()
{
  G4cout << " Cached field: " << G4endl
	  << "   Number of calls:        " << fCountCalls << G4endl
	  << "   Number of evaluations : " << fCountEvaluations << G4endl;
}

inline void  GetFieldValue( const G4double Point[4],
		G4double *Bfield ) const
{
	G4ThreeVector newLocation( Point[0], Point[1], Point[2] );

	// G4cout << "Cache-B-field called at " << newLocation << G4endl;

	G4double      distSq= (newLocation-fLastLocation).mag2();
	fCountCalls++;
	if( distSq < fDistanceConst*fDistanceConst ) { 
		Bfield[0] = fLastValue.x();
		Bfield[1] = fLastValue.y();
		Bfield[2] = fLastValue.z();
	}else{
		// G4CachedMagneticField* thisNonC= const_cast<G4CachedMagneticField*>(this);
		fpMagneticField->T_Field::GetFieldValue( Point, Bfield );
		// G4cout << " Evaluating. " << G4endl;
		fCountEvaluations++;
		// thisNonC->
		fLastLocation= G4ThreeVector( Point[0],  Point[1],  Point[2] );
		// thisNonC->
		fLastValue=    G4ThreeVector( Bfield[0], Bfield[1], Bfield[2] );
	}
}

G4double GetConstDistance() const         { return fDistanceConst; } 
void     SetConstDistance( G4double dist ){ fDistanceConst= dist;}

G4int GetCountCalls() const { return fCountCalls; }
G4int GetCountEvaluations() const { return fCountEvaluations; } 
void  ClearCounts() { fCountCalls = 0; fCountEvaluations=0; }


TCachedMagneticField& operator = (const TCachedMagneticField &p)
{
	if (&p == this) return *this; *this = p; return *this;
}


private: 
T_Field *fpMagneticField;
// When the field is evaluated within this distance it will not change
G4double      fDistanceConst;
// Caching state
mutable G4ThreeVector fLastLocation;
mutable G4ThreeVector fLastValue;

protected:
mutable G4int fCountCalls, fCountEvaluations;  
};

#endif /* TCACHED_MAGNETIC_FIELD_DEF */
