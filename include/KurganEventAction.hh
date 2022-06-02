#ifndef KurganEventAction_h
#define KurganEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4AccumulableManager.hh"
#include "globals.hh"
#include <vector>

class KurganRunAction;
using namespace std; 
using std::vector;

class KurganEventAction : public G4UserEventAction
{
  public:
    KurganEventAction(KurganRunAction* runAction);
    virtual ~KurganEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

	 void AddPositionDet1(G4double X, G4double Y){fPosX1 = X; fPosY1 = Y;}
	 void AddPositionDet2(G4double X, G4double Y){fPosX2 = X; fPosY2 = Y;}
    void AddStrawNumber(G4int StrawNumber){ 
		if( StrawNumber != tmp ){
		  fCounter += 1;
			  if(StrawNumber == 601){fCounter1 = fCounter - 1;}
		  }
		 tmp = StrawNumber;
		}

private:
		KurganRunAction*	fRunAction;
		G4double			fPosX1 = 0;
		G4double			fPosY1 = 0;
		G4double			fPosX2 = 0;
		G4double			fPosY2 = 0;
		G4int				fCounter = 0;
		G4int				fCounter1 = 0;
		G4int				tmp = 0;
};
#endif
