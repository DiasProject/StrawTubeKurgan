#include "KurganSteppingAction.hh"
#include "KurganEventAction.hh"
#include "KurganDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4TrackStatus.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "KurganAnalysis.hh"

#include "fstream"
#include "string"

using namespace std;

KurganSteppingAction::KurganSteppingAction(
			const KurganDetectorConstruction* detectorConstruction,
			KurganEventAction* eventAction) 
: G4UserSteppingAction(),
  fPositionDet1(eventAction),
  fPositionDet2(eventAction),
  fStrawNumber(eventAction),
//  fCounter(eventAction),
//  fCounter2(eventAction),
  fDetConstruction(detectorConstruction)
{}
void KurganSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if(volume != fDetConstruction->GetScoringVolume() 
		&& volume != fDetConstruction->GetSensetiveDetector()
		  && volume != fDetConstruction->GetSensetiveDetector2()){return;}
    
	G4Track* track = step->GetTrack();
	G4ThreeVector position = track->GetPosition();
	G4int Straw_Number = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();//номер строу трубки через которую прошла частица
//	G4int GetPDG =  step->GetTrack()->GetDefinition()->GetPDGEncoding();  
	G4int ParentID =  step->GetTrack()->GetParentID();  
// G4int Idtrack = step->GetTrack()->GetTrackID();//получение ID трека
//	G4int FirstStep = step->IsFirstStepInVolume();
	G4ThreeVector deltaPos = step->GetDeltaPosition();
  
  //cout<<"Volume="<<FirstStep<<" Copy_Number="<<Straw_Number<<endl;

	G4StepPoint* preStepPoint = step->GetPreStepPoint();//индикатор входа частицы в строу камеру
	G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
	G4ThreeVector worldPosition = preStepPoint->GetPosition();
	
if(ParentID==0){
/*  if(Straw_Number != 601){ 
			if(ParentID==0){
	//  	if(FirstStep==1){
#ifdef SHOW_COMMENTS_SteppingAction_cc
	cout<<Straw_Number<<" "<<position.x()<<" "<<position.y()<<" "<<position.z()<<" "<<Z_pos_part<<endl;
#endif
			fStrawNumber->AddStrawNumber(Straw_Number);
		  //fPositionDet->AddPositionDet(position.x(), position.y());	
		}
	}else */
  if(Straw_Number == 601){
    if(preStepPoint){
		 //fCounter->AddStrawNumber(Straw_Number);
		 fPositionDet1->AddPositionDet1(position.x(), position.y());	
     //cout<<"posX1="<<position.x()<<" posY1="<<position.y()<<endl;
    }
	}else if(Straw_Number == 602){
    if(preStepPoint){
		 fPositionDet2->AddPositionDet2(position.x(), position.y());	
     //cout<<"posX2="<<position.x()<<" posY2="<<position.y()<<endl;
    }
  }else{
    return;
  }
}

}
KurganSteppingAction::~KurganSteppingAction()
{}
