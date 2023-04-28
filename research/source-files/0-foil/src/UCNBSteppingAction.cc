#include "UCNBSteppingAction.hh"
#include "UCNBEventAction.hh"
#include "UCNBAnalysisManager.hh"
#include "G4SteppingManager.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "TStopwatch.h"
#include "TMath.h"
#include "G4EventManager.hh"
#include "G4SystemOfUnits.hh"
//----------------------------------------------------------------------------
UCNBSteppingAction::UCNBSteppingAction()
{
}
//----------------------------------------------------------------------------
void UCNBSteppingAction::UserSteppingAction(const G4Step* fStep)
{
  G4Track* fTrack = fStep->GetTrack();
  G4double eventTimeSoFar = ((UCNBEventAction*)G4EventManager::GetEventManager()->GetUserEventAction())->getEventCPUTime();
 // G4cout<<"Line 25 Stepping action **********************************"<<G4endl;
  if(eventTimeSoFar > 40.) {
    G4int iFlag = 1;
    UCNBAnalysisManager::getInstance()->killEventFlag(iFlag);
    fTrack->SetTrackStatus(fStopAndKill);
  }
#define COLLECT
#ifdef COLLECT

  // Determine proton total time-of-flight and hit positions in Silicon detectors/DeadLayers
  ///////////////////////////////////////////////////////////////////////////////////////////
  if ( ((fStep->GetTrack()->GetTrackID() == 2)) )
  {
    G4Track* track = fStep -> GetTrack();
    const G4DynamicParticle* dynParticle = track -> GetDynamicParticle();
    G4ParticleDefinition* particle = dynParticle -> GetDefinition();
    G4String particleName = particle -> GetParticleName();
    G4double dTstep = fStep->GetDeltaTime();
    UCNBAnalysisManager::getInstance()->AddUpProtonDriftTime(dTstep/s*1e9);

    if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon1") {
      G4ThreeVector hitSilicon1 = fStep->GetPreStepPoint()->GetPosition();
      G4double xhitSilicon1 = hitSilicon1.x()/m;
      G4double yhitSilicon1 = hitSilicon1.y()/m;
      G4double zhitSilicon1 = hitSilicon1.z()/m;
      UCNBAnalysisManager::getInstance()->recordSilicon1pPosition(xhitSilicon1,yhitSilicon1,zhitSilicon1);
    }

    if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon2") {
      G4ThreeVector hitSilicon2 = fStep->GetPreStepPoint()->GetPosition();
      G4double xhitSilicon2 = hitSilicon2.x()/m;
      G4double yhitSilicon2 = hitSilicon2.y()/m;
      G4double zhitSilicon2 = hitSilicon2.z()/m;
      UCNBAnalysisManager::getInstance()->recordSilicon2pPosition(xhitSilicon2,yhitSilicon2,zhitSilicon2);
    }
//---------------------------Dead Layer-----------------------------------------////////
  /*  if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Dead1") {
      G4ThreeVector hitSilicon1d = fStep->GetPreStepPoint()->GetPosition();
      G4double xhitSilicon1 = hitSilicon1d.x()/m;
      G4double yhitSilicon1 = hitSilicon1d.y()/m;
      G4double zhitSilicon1 = hitSilicon1d.z()/m;
      UCNBAnalysisManager::getInstance()->recordDead1pPosition(xhitSilicon1,yhitSilicon1,zhitSilicon1);
    }

    if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Dead2") {
      G4ThreeVector hitSilicon2d = fStep->GetPreStepPoint()->GetPosition();
      G4double xhitSilicon2 = hitSilicon2d.x()/m;
      G4double yhitSilicon2 = hitSilicon2d.y()/m;
      G4double zhitSilicon2 = hitSilicon2d.z()/m;
      UCNBAnalysisManager::getInstance()->recordDead2pPosition(xhitSilicon2,yhitSilicon2,zhitSilicon2);
    }
  }*/
  }
  //////////////////////////////////////////////////////////////
//---------------Electron Position and Time------------------------------------------
  if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon1")
  {
    G4Track* track = fStep -> GetTrack();
    const G4DynamicParticle* dynParticle = track -> GetDynamicParticle();
    G4ParticleDefinition* particle = dynParticle -> GetDefinition();
    G4String particleName = particle -> GetParticleName();
    if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
      G4ThreeVector hitSilicon11 = fStep->GetPreStepPoint()->GetPosition();
      G4double xhitSilicon1 = hitSilicon11.x()/m;
      G4double yhitSilicon1 = hitSilicon11.y()/m;
      G4double zhitSilicon1 = hitSilicon11.z()/m;
      UCNBAnalysisManager::getInstance()->recordSilicon1ePosition(xhitSilicon1,yhitSilicon1,zhitSilicon1);
    }
  }

  if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon2")
  {
    G4Track* track = fStep -> GetTrack();
    const G4DynamicParticle* dynParticle = track -> GetDynamicParticle();
    G4ParticleDefinition* particle = dynParticle -> GetDefinition();
    G4String particleName = particle -> GetParticleName();
    if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
      G4ThreeVector hitSilicon22 = fStep->GetPreStepPoint()->GetPosition();
      G4double xhitSilicon2 = hitSilicon22.x()/m;
      G4double yhitSilicon2 = hitSilicon22.y()/m;
      G4double zhitSilicon2 = hitSilicon22.z()/m;
      UCNBAnalysisManager::getInstance()->recordSilicon2ePosition(xhitSilicon2,yhitSilicon2,zhitSilicon2);
    }
  }

//----------------------------------------------------------------------------
  ////////   e- / daughter position tracking ///////////////////
   //check if while dereferencing they are null or not. When particle leaves the world there is no valid pointer and there is crash in Get Name.
   //Hence, one could check as in https://geant4-forum.web.cern.ch/t/getting-particles-out-of-the-world/3935/10
   //X == step->GetTrack()->GetNextVolume()  If the volume exists, then X is non-zero; if the volume does not exist, then X is zero.
   if (fStep->GetPostStepPoint()->GetPhysicalVolume()){
   if(fStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Dead1" && fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Drift1"){
        G4cout<<"Entering detector 1 "<<G4endl;G4Track* track = fStep -> GetTrack();
        const G4DynamicParticle* dynParticle = track -> GetDynamicParticle();
        G4ParticleDefinition* particle = dynParticle -> GetDefinition();
        G4String particleName = particle -> GetParticleName();
        if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
    
          G4ThreeVector pIncidentDet1 = fStep->GetPostStepPoint()->GetMomentumDirection();
          G4double pInDet1x = pIncidentDet1.x();
          G4double pInDet1y = pIncidentDet1.y();
          G4double pInDet1z = pIncidentDet1.z();
          G4cout <<" pinx -- befre run ---: "<<pInDet1x<<G4endl;
          UCNBAnalysisManager::getInstance()->p1incident(pInDet1x, pInDet1y, pInDet1z);
        }

      }
      /*finding outoing angle*/
      if(fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Dead1" && fStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Drift1"){
        G4cout<<"Leaving detector 1 "<<G4endl;
        if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
    
          G4ThreeVector pOutgoingDet1 = fStep->GetPreStepPoint()->GetMomentumDirection();
          G4double pOutDet1x = pOutgoingDet1.x();
          G4double pOutDet1y = pOutgoingDet1.y();
          G4double pOutDet1z = pOutgoingDet1.z();
          UCNBAnalysisManager::getInstance()->p1out(pOutDet1x, pOutDet1y, pOutDet1z);
        }
      }
      if(fStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Dead2" && fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Drift2"){
        G4cout<<"Entering detector 2 "<<G4endl;
         if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
    
            G4ThreeVector pIncidentDet2 = fStep->GetPostStepPoint()->GetMomentumDirection();
            G4double pInDet2x = pIncidentDet2.x();
            G4double pInDet2y = pIncidentDet2.y();
            G4double pInDet2z = pIncidentDet2.z();
            UCNBAnalysisManager::getInstance()->p2incident(pInDet2x, pInDet2y, pInDet2z);

        }
      }
/*    finding outoing angle*/
      if(fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Dead2" && fStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Drift2"){
        G4cout<<"Leaving detector 2 "<<G4endl;
        if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
    
          G4ThreeVector pOutgoingDet2 = fStep->GetPreStepPoint()->GetMomentumDirection();
          G4double pOutDet2x = pOutgoingDet2.x();
          G4double pOutDet2y = pOutgoingDet2.y();
          G4double pOutDet2z = pOutgoingDet2.z();
          UCNBAnalysisManager::getInstance()->p2out(pOutDet2x, pOutDet2y, pOutDet2z);

         }  
      }

   }
  if(fStep->GetTrack()->GetTrackID()==1){   // It's the primary electron
    if(fStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()=="e-"){
    // When no detector has been hit (eStop==-1), check for hit on silicon. When it happens, say the electron stopped 
      if((fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Silicon1"||fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Silicon2") && UCNBAnalysisManager::getInstance()->eStop==-1) { UCNBAnalysisManager::getInstance()->eStop=1;}
    // When a detector has been hit (eStop==1), check if the electron leaves the silicon. If so, say it did
      if(!(fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Silicon1"||fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Silicon2") && UCNBAnalysisManager::getInstance()->eStop==1) { UCNBAnalysisManager::getInstance()->eStop=0;}
    // When no detector has been hit (eStop==-1), record the energy at the start of the step.
      if(!(fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Silicon1"||fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Dead1"||fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Silicon2"||fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="Dead2") && UCNBAnalysisManager::getInstance()->eStop==-1) UCNBAnalysisManager::getInstance()->ePreSi=fStep->GetPreStepPoint()->GetKineticEnergy()/keV;
 //       G4cout<<" Breaking hereeee -- ln 129--------- "<<G4endl;
    }
  }


  // Give PID value to each track.
  if(UCNBAnalysisManager::getInstance()->PIDi[fStep->GetTrack()->GetTrackID()]==0) UCNBAnalysisManager::getInstance()->PIDi[fStep->GetTrack()->GetTrackID()]=fStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();
  
  //  Daughter photon tracking
    if(fStep->GetTrack()->GetTrackID()!=1){ 
      if(fStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()=="gamma" && fStep->GetTrack()->GetTrackID()!=UCNBAnalysisManager::getInstance()->gammaID){  
      UCNBAnalysisManager::getInstance()->gammaID=fStep->GetTrack()->GetTrackID();
      UCNBAnalysisManager::getInstance()->numGamma++;
 //G4cout<<" Breaking here -ln 142 photon if loop ---------- "<<G4endl;
      // If it's a gamma, record the initial state of the photon
      G4ThreeVector photonStart = fStep->GetPreStepPoint()->GetPosition();
      G4double xpho0 = photonStart.x()/m;
      G4double ypho0 = photonStart.y()/m;
      G4double zpho0 = photonStart.z()/m;
      G4ThreeVector photonPStart = fStep->GetPreStepPoint()->GetMomentumDirection();
      G4double pXpho0 = photonPStart.x();
      G4double pYpho0 = photonPStart.y();
      G4double pZpho0 = photonPStart.z();
      UCNBAnalysisManager::getInstance()->recordBremPos(xpho0, ypho0, zpho0,pXpho0,pYpho0,pZpho0);
      }
    }
/*==================== ELECTRON HITS EDEP TIME ================*/     


//Separating  out the "hits" by comparing the timme of a given hit to the time of the last hit

  if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon1" || fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon2")
  {
    // G4int TotalNoHits = 0;
	  //use same iterator:
    /*  SELECTING ELECTRON  */
	  if(fStep->GetTrack()->GetTrackID()==1){ 
      G4Track* track = fStep -> GetTrack();
      const G4DynamicParticle* dynParticle33 = track -> GetDynamicParticle();
      G4ParticleDefinition* particle33 = dynParticle33 -> GetDefinition();
      G4String particleName33 = particle33 -> GetParticleName();
     
    /*-----DETECTOR 1 --------------*/                                                                                                                              
      if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon1" ) {    
       // G4cout<<" DEad layer is at : "<<zDeadLayer/m<<G4endl;
        G4ThreeVector hitDetector1 = fStep->GetPreStepPoint()->GetPosition();
        G4double zHitDetector1 = hitDetector1.z()/m; // all of these z's are in mm and not in m 
      //  G4cout<<" position it hits : "<<zHitDetector1 <<G4endl;

      //  if(zHitDetector1 > (zDeadLayer/m)){  
        //  G4cout<<"1"<<G4endl;
          if (UCNBAnalysisManager::getInstance()->dESi1HitTime[0]==0){
            UCNBAnalysisManager::getInstance()->HitNo1 = 0;
            UCNBAnalysisManager::getInstance()->Det1Hits = UCNBAnalysisManager::getInstance()->HitNo1 + 1;
            UCNBAnalysisManager::getInstance()->dESi1HitTime[0] =fStep->GetPreStepPoint()->GetGlobalTime()/nanosecond;
          }
          else if ((fStep->GetPreStepPoint()->GetGlobalTime()/nanosecond) - (UCNBAnalysisManager::getInstance()->dESi1HitTime[UCNBAnalysisManager::getInstance()->HitNo1]) > 10.){//need to hve anther hit statement because in the next nteraction it is coming to the following else if statement. Hence need to have 
            UCNBAnalysisManager::getInstance()->HitNo1++;
            UCNBAnalysisManager::getInstance()->Det1Hits = UCNBAnalysisManager::getInstance()->HitNo1 + 1;
            UCNBAnalysisManager::getInstance()->dESi1HitTime[UCNBAnalysisManager::getInstance()->HitNo1]=fStep->GetPreStepPoint()->GetGlobalTime()/nanosecond;
          }  
      //  } 
      } 
 
    /*------DETECTOR 2 --------*/
      if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon2") {
        G4ThreeVector hitDetector2 = fStep->GetPreStepPoint()->GetPosition();
        G4double zHitDetector2 = -1.*hitDetector2.z()/m; // all of these z's are in mm and not in m 
        G4double dEstep2 = fStep->GetTotalEnergyDeposit();
       // if(zHitDetector2 > zDeadLayer/m){  
        //  G4cout<<"2"<<G4endl;
          if (UCNBAnalysisManager::getInstance()->dESi2HitTime[0]==0){
            UCNBAnalysisManager::getInstance()->HitNo2 = 0;
            UCNBAnalysisManager::getInstance()->Det2Hits = UCNBAnalysisManager::getInstance()->HitNo2 + 1;
            UCNBAnalysisManager::getInstance()->dESi2HitTime[0] =fStep->GetPreStepPoint()->GetGlobalTime()/nanosecond;
          }
          else if ((fStep->GetPreStepPoint()->GetGlobalTime()/nanosecond) - (UCNBAnalysisManager::getInstance()->dESi2HitTime[UCNBAnalysisManager::getInstance()->HitNo2]) > 10.){//need to hve anther hit statement because in the next nteraction it is coming to the following else if statement. Hence need to have 
            UCNBAnalysisManager::getInstance()->HitNo2++;
            UCNBAnalysisManager::getInstance()->Det2Hits = UCNBAnalysisManager::getInstance()->HitNo2 + 1;
            UCNBAnalysisManager::getInstance()->dESi2HitTime[UCNBAnalysisManager::getInstance()->HitNo2]=fStep->GetPreStepPoint()->GetGlobalTime()/nanosecond;
          } 
       // } 
      }
    UCNBAnalysisManager::getInstance()->TotalNoHits =  UCNBAnalysisManager::getInstance()->Det2Hits+ UCNBAnalysisManager::getInstance()->Det1Hits ;
    //G4cout<<"Total number of hits   : "<<UCNBAnalysisManager::getInstance()->TotalNoHits<<G4endl;
  }

//--------------Energy deposition in Silicon detectors and deadlayers-----------------------------------------
  if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Dead1") {
      G4ThreeVector hitSilicon1d = fStep->GetPreStepPoint()->GetPosition();
      if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
        G4double dEstep = fStep->GetTotalEnergyDeposit();
        UCNBAnalysisManager::getInstance()->AddUpElectronDeadLayer1EnergyDeposition(dEstep/keV);
       }
  }

 //G4cout<<" Breaking here ln 225 ----------- "<<G4endl;
  if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Dead2") {
      G4ThreeVector hitSilicon1d = fStep->GetPreStepPoint()->GetPosition();
      if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
        G4double dEstep = fStep->GetTotalEnergyDeposit();
        UCNBAnalysisManager::getInstance()->AddUpElectronDeadLayer2EnergyDeposition(dEstep/keV);
       }
  }

// G4cout<<" Breaking here ln 234 ----------- "<<G4endl;
  // Add up energy deposition in Silicon Detector #1
  if(fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon1")
  {
    G4Track* track = fStep -> GetTrack();
    const G4DynamicParticle* dynParticle1 = track -> GetDynamicParticle();
    G4ParticleDefinition* particle = dynParticle1 -> GetDefinition();
    G4String particleName = particle -> GetParticleName();

 //G4cout<<" Breaking here ln 243 ----------- "<<G4endl;
    if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
      G4double dEstep = fStep->GetTotalEnergyDeposit();
      G4double timeHit1 = track->GetGlobalTime();
      G4double test = UCNBAnalysisManager::getInstance()->globalTimeHit1;
   
// G4cout<<" Breaking here ln 249 ----------- "<<G4endl;  
      if (dEstep/keV > 0. && timeHit1/s > test) { //it seems here that getting the global time since the track is created. and (once done with the event 7 deeper anaysis check what the global timehit1 is . is it initialised at some point because it is just been set equal to somwthign)
      	UCNBAnalysisManager::getInstance()->globalTimeHit1 = timeHit1/s;
       // G4cout<<" Time hit det 1 ============================"<<timeHit1/s<<G4endl;
      }
    }

// G4cout<<" Breaking here ln 256 ----------- "<<G4endl;
    if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ){
      G4Track* track1 = fStep -> GetTrack();
      const G4DynamicParticle* dynParticle = track1 -> GetDynamicParticle();
      G4ParticleDefinition* particle1 = dynParticle -> GetDefinition();
      G4String particleName1 = particle1 -> GetParticleName();
      G4ThreeVector hitDetector1 = fStep->GetPreStepPoint()->GetPosition();
      G4double zHitDetector1 = hitDetector1.z()/m; // all of these z's are in mm and not in m 
     // G4double dEstep = fStep->GetTotalEnergyDeposit();
    //  G4cout<<" position it hits : "<<zHitDetector1<<G4endl;
      G4double dEstep = fStep->GetTotalEnergyDeposit();
 
 //G4cout<<" Breaking here ln 268 ----------- "<<G4endl;       
    //  if(zHitDetector1 > (zDeadLayer/m)){  
        UCNBAnalysisManager::getInstance()->AddUpElectronSilicon1EnergyDeposition(dEstep/keV);
    //  }else {
    //    UCNBAnalysisManager::getInstance()->AddUpElectronDeadLayer1EnergyDeposition(dEstep/keV);
    //  }
      G4double timeHit1 = track1->GetGlobalTime();
      timeHit1 = timeHit1/s * 1.0e9;
      G4int timeInt = (G4int) timeHit1;
      if (timeInt > 499.) timeInt = 499;
    //  UCNBAnalysisManager::getInstance()->EdepTimeBin1[timeInt] += dEstep/keV;
    }

    
 //G4cout<<" Breaking here ln 282 ----------- "<<G4endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
  // Add up energy deposition in Silicon Detector #2
/////////////////////////////////////////////////////////////////////////
  if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Silicon2")
  {
    G4Track* track = fStep -> GetTrack();
    const G4DynamicParticle* dynParticle = track -> GetDynamicParticle();
    G4ParticleDefinition* particle = dynParticle -> GetDefinition();
    G4String particleName = particle -> GetParticleName();

// G4cout<<" Breaking here ln 294 ----------- "<<G4endl;
    if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
      G4double dEstep = fStep->GetTotalEnergyDeposit();
      G4double test = UCNBAnalysisManager::getInstance()->globalTimeHit2;
      G4double timeHit2 = track->GetGlobalTime();
      if (dEstep/keV > 0. && timeHit2/s > test) {
	      UCNBAnalysisManager::getInstance()->globalTimeHit2 = timeHit2/s;
      }
    }
  
    if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) ) {
      G4double dEstep = fStep->GetTotalEnergyDeposit();
      G4ThreeVector hitDetector2 = fStep->GetPreStepPoint()->GetPosition();
      G4double zHitDetector2 = -1.*hitDetector2.z()/m; // all of these z's are in mm and not in m 
    //  if(zHitDetector2 > zDeadLayer/m){  
        UCNBAnalysisManager::getInstance()->AddUpElectronSilicon2EnergyDeposition(dEstep/keV);
    //  }
    //  else{
    //  G4cout<<" hitting the dead layer det 2 at : "<<zHitDetector2<<" with energy : "<< dEstep/keV<<G4endl;
    //    UCNBAnalysisManager::getInstance()->AddUpElectronDeadLayer2EnergyDeposition(dEstep/keV);
    //  }    
      G4double timeHit2 = track->GetGlobalTime();
      timeHit2 = timeHit2/s * 1.0e9;
      G4int timeInt = (G4int) timeHit2;
      if (timeInt > 499.) timeInt = 499;
      UCNBAnalysisManager::getInstance()->EdepTimeBin2[timeInt] += dEstep/keV;
    }

  }

  // Add up energy loss elsewhere (including to secondaries)
// G4cout<<" Breaking here ln 323 ----------- "<<G4endl;
  if ( (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "Silicon1") &&
       (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "Silicon2") 
       )
  {
    G4Track* track = fStep -> GetTrack();
    const G4DynamicParticle* dynParticle = track -> GetDynamicParticle();
    G4ParticleDefinition* particle = dynParticle -> GetDefinition();
    G4String particleName = particle -> GetParticleName();
  //   if ( ((fStep->GetTrack()->GetTrackID() != 2)) && ((fStep->GetTrack()->GetParentID() != 2)) )
    if ( ((fStep->GetTrack()->GetTrackID() == 1)) )
    {
      G4double dEstep = fStep->GetTotalEnergyDeposit();
      UCNBAnalysisManager::getInstance()->AddUpElectronOtherEnergyDeposition(dEstep/keV);
    }
     //if ( ((fStep->GetTrack()->GetTrackID() == 2)) || ((fStep->GetTrack()->GetParentID() == 2)) )
 
  }

#endif
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
