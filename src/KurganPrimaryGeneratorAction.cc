#include "KurganPrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

KurganPrimaryGeneratorAction::KurganPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  //G4ParticleDefinition* particle = particleTable->FindParticle(particleName="kaon+");
  //G4ParticleDefinition* particle = particleTable->FindParticle(particleName="proton");
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="mu+");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //fParticleGun->SetParticleEnergy(75.*GeV);
  fParticleGun->SetParticleEnergy(10.*MeV);
}

KurganPrimaryGeneratorAction::~KurganPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void KurganPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//точечный источник
//  G4ThreeVector iso2;
//  G4double costet = G4RandFlat::shoot(0.996194698,1.0);
//  G4double tet = acos(costet);
//  G4double phi = G4RandFlat::shoot(2*CLHEP::pi);
//  iso2.setRThetaPhi(1.,/*(360.0-3.5)*deg*/tet,/*0.*deg*/phi);
//	 static const double pi = 3.14159265358979323846;
//	 G4cout<<"tet:"<<tet*180/pi<<" costet:"<<acos(costet)*180/pi<<" phi:"<<phi<<G4endl;
//  fParticleGun->SetParticleMomentumDirection(iso2.unit());
//  fParticleGun->SetParticlePosition(G4ThreeVector(15.0*cm, 50.0*cm, -10.0*cm));
//  fParticleGun->GeneratePrimaryVertex(anEvent);

//кругляй источник
   //G4ThreeVector iso2;
	//G4double x = (G4RandGauss::shoot(10.0, 0.16985))*mm;
	//G4double y = (G4RandGauss::shoot(0.0, 0.16985))*mm;
	//G4double z = -50*cm;
   //fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
   //fParticleGun->GeneratePrimaryVertex(anEvent);
 
//квадратный источник
//	G4double x = -105.0*cm + 210.0*(1.*G4UniformRand()-0.)*cm;
	G4double y = -100.0*cm + 200.0*(1.*G4UniformRand()-0.)*cm;
	//G4double x = -1.0*cm + 2.0*(1.*G4UniformRand()-0.)*cm;
	//G4double x = -100.0*cm + 90.0*(1.*G4UniformRand()-0.)*cm;
	G4double x = -100.0*cm + 200.0*(1.*G4UniformRand()-0.)*cm;
	//G4double y = -70.0*cm + 15.0*(1.*G4UniformRand()-0.)*cm;
	G4double z = -50*cm;

	G4ThreeVector iso2;
	//G4double costet = G4RandFlat::shoot(0.996194698,1.0);//0-5 degree
	//G4double costet = G4RandFlat::shoot(0.99875026,1.0);//0-2.86 degree
	G4double costet = G4RandFlat::shoot(1.0,1.0);//0 degree
	G4double tet = acos(costet);
	G4double phi = G4RandFlat::shoot(2*CLHEP::pi);
	iso2.setRThetaPhi(1.,/*(360.0-3.5)*deg*/tet,/*0.*deg*/phi);
	fParticleGun->SetParticleMomentumDirection(iso2.unit());
   fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
   fParticleGun->GeneratePrimaryVertex(anEvent);

 
}
