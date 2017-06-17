#include "MaterialConverter.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "TGeoManager.h"
#include "TGeoMaterial.h"

// #include "assert.h"
#include <cstdio>

TGeoManager* MaterialConverter::fgTGeomMgr = 0;
MaterialConverter* MaterialConverter::fgInstance = 0;

MaterialConverter::MaterialConverter()
{ 
//  const G4MaterialTable *theG4MaterialTable = G4Material::GetMaterialTable();
//  G4int nmaterials = theG4MaterialTable->size();
//  fRootMatIndices.reserve(defSize); 
//  fG4MatIndices.reserve(defSize); 
  if(!fgTGeomMgr) {
   G4cerr << "MaterialConverter::fgTGeomMgr is NULL. Must be set before calling the constructor."<<G4endl;
   exit(EXIT_FAILURE);
  }
  Initialize();
}

MaterialConverter::~MaterialConverter() { }

MaterialConverter* MaterialConverter::Instance()
{
    if(!fgInstance) fgInstance = new MaterialConverter();
    
    return fgInstance;
}


void MaterialConverter::Initialize()
{
  // Build the list of materials in Root, given the materials in Geant4

  const G4MaterialTable *theG4MaterialTable = G4Material::GetMaterialTable();

  int nmaterials= theG4MaterialTable->size();
  std::cout<<"---- Number of G4Materials found = "<<nmaterials<<std::endl;

  for(G4int imatG4=0; imatG4<nmaterials; ++imatG4) {
     G4Material *g4mat = (*theG4MaterialTable)[imatG4];
     G4int numElements= g4mat->GetNumberOfElements(); 

     // Create corresponding TGeoMaterial
     
     if( numElements == 1 ) {
            new TGeoMaterial(g4mat->GetName(), g4mat->GetA(), g4mat->GetZ(), 
                             g4mat->GetDensity()*CLHEP::cm3/CLHEP::g,   // Units => Root Units ?
                             g4mat->GetRadlen(), 
                             g4mat->GetNuclearInterLength() );
     } else {
        // G4ElementVector* elementVec= g4mat->GetElementVector();
        const G4double*  g4elemFractions= g4mat->GetFractionVector();
        TGeoMixture *tgeoMixture = new TGeoMixture(g4mat->GetName(), numElements, g4mat->GetDensity()*CLHEP::cm3/CLHEP::g );      
        for( int ielem = 0; ielem < numElements ; ielem++ )
        {
           const G4Element* g4elem= g4mat->GetElement(ielem);
           tgeoMixture->AddElement(g4elem->GetA(), g4elem->GetZ(), g4elemFractions[ielem]);
        }

     }

     G4int rtMatIdx = fgTGeomMgr->GetMaterialIndex(g4mat->GetName());
     fgTGeomMgr->GetMaterial(rtMatIdx)->SetUsed();

//     std::cout<<"---   "<<imatG4<< " rtMatIdx =  "<<rtMatIdx<< " rtName  ="<<
//       fgTGeomMgr->GetMaterial(rtMatIdx)->GetName() <<" g4name = "<< g4mat->GetName()<< "  #elements ="<<numElements<< std::endl; 
     
  }
}


