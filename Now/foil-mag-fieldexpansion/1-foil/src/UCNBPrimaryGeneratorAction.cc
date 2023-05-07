#include "UCNBPrimaryGeneratorAction.hh"
#include "UCNBDetectorConstruction.hh"
#include "UCNBAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4VEmModel.hh"
 #include "globals.hh"

UCNBPrimaryGeneratorAction::UCNBPrimaryGeneratorAction(UCNBDetectorConstruction* myDC)
  :UCNBDetector(myDC)
{
}

UCNBPrimaryGeneratorAction::~UCNBPrimaryGeneratorAction()
{
  delete particleGun;
}

void UCNBPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 // G4cout<<" Generating primaries    :PGA"<<G4endl;
/*generating energies for 113Sn
From : https://pdg.lbl.gov/2021/reviews/contents_sports.html 
EC == electron capture 
EC1 = 364 keV intensity 28 % 
EC2 = 388 keV intensit 6 %
Idea :  if we generate random number between 1 to 34 (28 + 6) the chances that EC1 will come out are between 1-28 */

  G4double EC1 = 364*keV;
  G4double EC2 = 388*keV;
  G4double chance = 34.0*G4UniformRand();
  G4double Ee;
  if (chance  <= 28)
  {
    Ee = EC1;
  }
  else Ee = EC2;
  
    G4double cosTHETAE = -1 + 2*G4UniformRand();
    G4double PHIE = 2.*M_PI*G4UniformRand();

  G4double EU = cos(PHIE)*sin(acos(cosTHETAE));
  G4double EV = sin(PHIE)*sin(acos(cosTHETAE));
  G4double EW = cosTHETAE;


  G4double px_hat_e = EU;
  G4double py_hat_e = EV;
  G4double pz_hat_e = EW;
  G4double thetaElectron = acos(pz_hat_e);
 

  // Sample vertex position
   G4double x_vertex, y_vertex, z_vertex;
   G4double x_test, y_test;
   G4double testRadius = 1.e99;
  // while (testRadius > 0.05) {
  //  x_test = (-0.05 + G4UniformRand()*2.*0.05)*m;
  //  y_test = (-0.05 + G4UniformRand()*2.*0.05)*m;
  //  testRadius = sqrt(x_test*x_test + y_test*y_test);
  //}
   x_vertex = 0*m;
   y_vertex = 0*m;
   //z_vertex = (-1.5 + G4UniformRand()*3.0)*m;
  //  z_vertex = 1.5*m + 280e-9*m + 0.5e-6*m; //1.5 m is the end f decay trap right after which is the foil of thickness 280 nm (one foil) and the source holder of thickness 0.5 micron
  z_vertex = -0.25e-06*m; //1.5 m is the end f decay trap right after which is the foil of thickness 280 nm (one foil) and the source holder of thickness 0.5 micron
  
  // Generate electron 
  G4int n_particle = 1;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particleGun = new G4ParticleGun(n_particle);
  G4ParticleDefinition* particle1 = particleTable->FindParticle("e-");
  particleGun->SetParticleDefinition(particle1);
  particleGun->SetParticleEnergy(Ee);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(px_hat_e,py_hat_e,pz_hat_e));

  particleGun->SetParticlePosition(G4ThreeVector(x_vertex,y_vertex,z_vertex));
  particleGun->GeneratePrimaryVertex(anEvent);

  G4double Tp0 = 0;
  G4double Tn0 = 0;
  G4double Tv0 = 0;
  G4double thetaProton = 0;
  // Save initial vertex variables
  UCNBAnalysisManager::getInstance()->saveEventVertex(x_vertex/m, y_vertex/m,z_vertex/m, Ee/keV, Tp0/keV,
						      px_hat_e,py_hat_e,pz_hat_e, thetaElectron, thetaProton, Tn0/keV, Tv0/keV);
 
  // G4cout<<"EE == "<<Ee<<G4endl;
   G4cout<<"z_e - "<<z_vertex<<G4endl;
}
