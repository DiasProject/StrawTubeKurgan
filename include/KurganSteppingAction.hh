#ifndef KurganSteppingAction_h
#define KurganSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

//class G4LogicalVolume;
class KurganEventAction;
class KurganDetectorConstruction;

using namespace std;
/// Stepping action class
/// 

class KurganSteppingAction : public G4UserSteppingAction
{
  public:
    KurganSteppingAction(const KurganDetectorConstruction* detectorConstruction, KurganEventAction* eventAction);
    virtual ~KurganSteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

   
  private:
    KurganEventAction*   fStrawNumber;
    KurganEventAction*   fPositionDet1;
    KurganEventAction*   fPositionDet2;
    KurganEventAction*   fCounter;
	 const KurganDetectorConstruction* fDetConstruction;

};
#endif
