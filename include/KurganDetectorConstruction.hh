#ifndef KurganDetectorConstruction_h
#define KurganDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#define M_PI 3.14159265358979323846

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class KurganDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    KurganDetectorConstruction();
    virtual ~KurganDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    
    const G4LogicalVolume*  GetScoringVolume() const;
    const G4LogicalVolume*  GetSensetiveDetector() const;
    const G4LogicalVolume*  GetSensetiveDetector2() const;
  
	private:
    G4LogicalVolume*  fScoringVolume;
    G4LogicalVolume*  fSensetiveDetector;
    G4LogicalVolume*  fSensetiveDetector2;

};

inline const G4LogicalVolume* KurganDetectorConstruction::GetScoringVolume() const{
return fScoringVolume;
}

inline const G4LogicalVolume* KurganDetectorConstruction::GetSensetiveDetector() const{
return fSensetiveDetector;
}

inline const G4LogicalVolume* KurganDetectorConstruction::GetSensetiveDetector2() const{
return fSensetiveDetector2;
}
#endif
