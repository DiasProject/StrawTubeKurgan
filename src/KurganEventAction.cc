#include "KurganAnalysis.hh"
#include "KurganEventAction.hh"
#include "KurganRunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AccumulableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VHitsCollection.hh"
#include "fstream"

#include "comments.hh"

using namespace std;
using std::array;
using std::vector;

KurganEventAction::KurganEventAction(KurganRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{} 

KurganEventAction::~KurganEventAction()
{
   delete G4AnalysisManager::Instance();
}

void KurganEventAction::BeginOfEventAction(const G4Event* event)
{    
//  G4int EventID = event->GetEventID(); 
}

void KurganEventAction::EndOfEventAction(const G4Event* event)
{
     // accumulate statistics in run action
     //fRunAction->AddStrawNumber_Run(fCounter);
     //fRunAction->AddCounter1_Run(counter1);
     //fRunAction->AddCounter2_Run(counter2);
     //fRunAction->AddPositionDet_Run(fPosX, fPosY);
     //fRunAction->AddPositionDet2_Run(fPosX2, fPosY2);
 
     // get analysis manager
     auto analysisManager = G4AnalysisManager::Instance();
 
     // fill histograms
     analysisManager->FillH1(0, fCounter);
     //analysisManager->FillH1(1, result1);
     //analysisManager->FillH1(2, result2);
     
     // fill histograms 2D
     analysisManager->FillH2(0, fPosX1, fPosY1);
     analysisManager->FillH2(1, fPosX2, fPosY2);

#ifdef SHOW_COMMENTS_EventAction_cc
     cout<<"fCounter1_srs_event="<<fCounter1<<endl;
     cout<<"fCounter="<<fCounter<<" EventAction"<<endl;
     G4cout<<"X="<<fPosX<<" Y="<<fPosY<<G4endl;
     G4cout<<"X2="<<fPosX2<<" Y2="<<fPosY2<<G4endl;
#endif 

     // fill ntuple
     analysisManager->FillNtupleDColumn(0, fCounter);
     //analysisManager->FillNtupleDColumn(1, counter1);
     //analysisManager->FillNtupleDColumn(2, counter2);
     analysisManager->AddNtupleRow();

}
