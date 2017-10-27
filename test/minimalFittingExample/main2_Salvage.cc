#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <PlanarMeasurement.h>

#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>
#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <cmath>


#include <TApplication.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TVector3.h>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include <TMath.h>
#include <TString.h>
#include <TEventList.h>

#include <memory>

#include <EventDisplay.h>
#include <AbsFinitePlane.h>
#include "AbsTrackRep.h"
#include <AbsFitterInfo.h>
#include <TGeoMaterialInterface.h>
#include <HelixTrackModel.h>
#include <MeasurementCreator.h>
#include <AbsMeasurement.h>
#include <AbsTrackRep.h>
#include <ConstField.h>
#include <DetPlane.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFittedStateOnPlane.h>
#include <AbsKalmanFitter.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitStatus.h>
#include <DAF.h>
#include <GFGbl.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>
#include <PlanarMeasurement.h>
#include <RectangularFinitePlane.h>
#include <ReferenceStateOnPlane.h>
#include <SharedPlanePtr.h>
#include <SpacepointMeasurement.h>
#include <StateOnPlane.h>
#include <Tools.h>
#include <TrackCand.h>
#include <TrackCandHit.h>
#include <Track.h>
#include <TrackPoint.h>
#include <WireMeasurement.h>
#include <WireMeasurementNew.h>

#include <MaterialEffects.h>
#include <RKTools.h>
#include <RKTrackRep.h>
#include <StepLimits.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

const double noHitMarker = 99999; //Marker given to identify when particle does not hit detector 
const int nDetBeam = 5; //number of detectors before the target. 2 SiPMs, 3 GEMs
const int nDetScatter = 8; //number of detectors after the target. 2 STTs, 2 scintillators on each side

std::vector<double> getHits(std::string branch);
std::vector<double> getCoord(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, int);
std::vector<int> getPositions();
void beamlinestuff();

int main() {
	beamlinestuff();
	/*
  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. ,10., 0.)); // 1 T


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack(); Uncomment for original


  // particle pdg code; pion hypothesis
  const int pdg = 211;

  // start values for the fit, e.g. from pattern recognition
  TVector3 pos(0, 0, 0);
  TVector3 mom(0, 0, 3);


  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

  // create track
  genfit::Track fitTrack(rep, pos, mom);
*/
/* Uncomment for original
  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID

  double detectorResolution(0.001); // resolution of planar detectors
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  hitCov *= detectorResolution*detectorResolution;


  // add some planar hits to track with coordinates I just made up
  TVectorD hitCoords(2);
  hitCoords[0] = 0;
  hitCoords[1] = 0;
  genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
  measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,0), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
  fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

  hitCoords[0] = -0.15;
  hitCoords[1] = 0.70;
  measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
  measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,10), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
  fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

  hitCoords[0] = -0.4;
  hitCoords[1] = 0;
  measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
  measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,20), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
  fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Modifications For GEMs
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/* 
  //Open files to write residuals
  std::ofstream fx;
  std::ofstream fy;
  std::ofstream fx2;
  std::ofstream fy2;
  fx.open("GEM_xres.txt");
  fy.open("GEM_yres.txt");
  fx2.open("GEM2_xres.txt");
  fy2.open("GEM2_yres.txt");

  //Getting hit coordinates from simulation
	std::vector<double> gem1x,gem1y,gem2x,gem2y,gem3x,gem3y,sipm1x,sipm1y,sipm2x,sipm2y,scl1x,scl1y,scl1z,scl2x,scl2y,scl2z,scr1x,scr1y,scr1z,scr2x,scr2y,scr2z,stcl1x,stcl1y,stcl1z,stcl2x,stcl2y,stcl2z,stcr1x,stcr1y,stcr1z,stcr2x,stcr2y,stcr2z; 

	gem1x = getHits("GEM1_PosHitX");
	gem1y = getHits("GEM1_PosHitY");
	gem2x = getHits("GEM2_PosHitX");
	gem2y = getHits("GEM2_PosHitY");
	gem3x = getHits("GEM3_PosHitX");
	gem3y = getHits("GEM3_PosHitY");
	sipm1x = getHits("BCC1_PosHitX");
	sipm1y = getHits("BCC1_PosHitY");
	sipm2x = getHits("BCC2_PosHitX");
	sipm2y = getHits("BCC2_PosHitY");
	scl1x = getHits("SCL1_PosHitX");
	scl1y = getHits("SCL1_PosHitY");
	scl1z = getHits("SCL1_PosHitZ");
	scl2x = getHits("SCL2_PosHitX");
	scl2y = getHits("SCL2_PosHitY");
	scl2z = getHits("SCL2_PosHitZ");
	scr1x = getHits("SCR1_PosHitX");
	scr1y = getHits("SCR1_PosHitY");
	scr1z = getHits("SCR1_PosHitZ");
	scr2x = getHits("SCR2_PosHitX");
	scr2y = getHits("SCR2_PosHitY");
	scr2z = getHits("SCR2_PosHitZ");
	stcl1x = getHits("STCL1_PosHitX");
	stcl1y = getHits("STCL1_PosHitY");
	stcl1z = getHits("STCL1_PosHitZ");
	stcl2x = getHits("STCL2_PosHitX");
	stcl2y = getHits("STCL2_PosHitY");
	stcl2z = getHits("STCL2_PosHitZ");
	stcr1x = getHits("STCR1_PosHitX");
	stcr1y = getHits("STCR1_PosHitY");
	stcr1z = getHits("STCR1_PosHitZ");
	stcr2x = getHits("STCR2_PosHitX");
	stcr2y = getHits("STCR2_PosHitY");
	stcr2z = getHits("STCR2_PosHitZ");

		
	//Create track in GenFit
  genfit::Track* beamlinetrack; //Track before target
  //genfit::Track* scatteredtrack; //Track after hitting target
  //genfit::EventDisplay* display;

	beamlinetrack = beamlinestuff(sipm1x, sipm1y, sipm2x, sipm2y, gem1x, gem1y, gem2x, gem2y, gem3x, gem3y);

  genfit::MaterialEffects::getInstance()->setNoEffects();
  //Turning off material effects (for now)
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,0.));
  //No mag field
  display = genfit::EventDisplay::getInstance();
  genfit::AbsKalmanFitter* fitter = new genfit::DAF(); //KalmanFitter(20,0.001);
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(20,0.001); //KalmanFitter(20,0.001);

	
  TVectorD hitCoords(2); //will hold x and y coordinate of single hit
  const int detId(0); //detector ID
  int planeId(0); //detector plane ID
  int hitId(0); //hit ID
	double smearVal(0.1); //width of applied gaussian smear (mm) centered at 0
  double detectorResolution(1); //resolution of GEMs detectors, mm
	int numtracks = gem3x.size(); //Number of tracks to try to reconstruct
	int GEM2(3), GEM3(4); //For easier access to coordinates in vectors
  
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
	hitCov *= detectorResolution*detectorResolution;

  genfit::AbsKalmanFitter* fitter = new genfit::DAF(); //KalmanFitter(20,0.001);
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(20,0.001); //KalmanFitter(20,0.001);

  //creating z coordinates of GEMs
  std::vector<int> detectorPositions; //hold positions of detectors
	detectorPositions = getPositions();

  //Loop over number of tracks
  for(int k=0; k < numtracks; k++)
  {  
		//Declaring a new track
		//The two TVector3s are posSeed and momSeed. Can not figure out the difference if they are changed,
		//e.g. the two lines below (one of which should remain commented) both produce the same output.
		//More info in GenFit/core/src/track.cc 
    //beamlinetrack = new genfit::Track(new genfit::RKTrackRep(pdg), pos, mom);
    beamlinetrack = new genfit::Track(new genfit::RKTrackRep(11), TVector3(125,125,1700), TVector3(0,0,1));
    std::vector<double> xcoords; //hold positions of hits
    std::vector<double> ycoords;
    
		xcoords = getCoord(sipm1x, sipm2x, gem1x, gem2x, gem3x, k); //X coordinates for 1 particle across all detectors 
		ycoords = getCoord(sipm2y, sipm2y, gem1y, gem2y, gem3y, k); //Y
		
			//put stuff here if things break, add curly bracket after //end loop to find hits, comment below until for loop 
			int noHitPos = 0; //will hold index of first detector marked as having no hit
			for(int it=0; it < nDetBeam; it++)
			{ 
				//This loop is meant to find the first detector that records no hit, store its position in 'noHit'
				//in order to only fit to the detectors that record hits. 
				if(xcoords.at(it) == noHitMarker && it < GEM3){
					noHitPos = it;
					break;
				}
				else if(it == GEM3 && xcoords.at(it) == noHitMarker){
					noHitPos = GEM3;
					break;
				}
				else if(it == GEM3 && xcoords.at(it) != noHitMarker){
					noHitPos = nDetBeam;
					break;
				}
			}//end loop to find noHitPos
			if(noHitPos == 0) continue;

    	for (int i=0; i<noHitPos; i++)
  		{
    		//Add Gaussian smear to the hits
    	  hitCoords[0] = xcoords[i]*(1+gRandom->Gaus(0,smearVal));
    	  hitCoords[1] = ycoords[i]*(1+gRandom->Gaus(0,smearVal));
    	  int z = (detectorPositions[i]);
  
	  	  genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
	  		measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z), TVector3(1,0,0), TVector3(0,1,0))), planeId++);
	  		beamlinetrack->insertPoint(new genfit::TrackPoint(measurement, beamlinetrack));
  		}//end loop to find hits

    if(beamlinetrack->getNumPoints()==0){
    	std::cout<< "No Points in beamlinetrack\n";
    	return 0; //nothing to fit
  	}
		
    assert(beamlinetrack->checkConsistency()); //check
    fitter->processTrack(beamlinetrack); //do the fit
    assert(beamlinetrack->checkConsistency());
    display->addTrack(beamlinetrack);
*/
//......................................................................................................
//.............................Trying to Extrapolate....................................................
//......................................................................................................	
	/*
	for(int poop=0; poop<nDetBeam; poop++){
    unsigned int numhits = beamlinetrack->getNumPointsWithMeasurement();
		//genfit::TrackPoint* tp = beamlinetrack->getPointWithMeasurement(poop); 

    genfit::KalmanFitterInfo* fi;
    genfit::KalmanFitterInfo* prevFi = 0;
    const genfit::MeasuredStateOnPlane* fittedState(NULL);
    const genfit::MeasuredStateOnPlane* prevFittedState(NULL);

    for(unsigned int j = 0; j < numhits; j++) { // loop over all hits in the track

      fittedState = NULL;

      genfit::TrackPoint* tp = beamlinetrack->getPointWithMeasurement(j);
      if (! tp->hasRawMeasurements()) {
        std::cerr<<"trackPoint has no raw measurements"<<std::endl;
        continue;
      }

      const genfit::AbsMeasurement* m = tp->getRawMeasurement();
      //int hit_coords_dim = m->getDim();

      // check if multiple AbsMeasurements are of same type
      if (tp->getNumRawMeasurements() > 1) {
        bool sameTypes(true);
        for (unsigned int iM=1; iM<tp->getNumRawMeasurements(); ++iM) {
          if (typeid(*(tp->getRawMeasurement(iM))) != typeid(*m))
            sameTypes = false;
        }
        if (!sameTypes) {
          std::cerr<<"cannot draw trackpoint containing multiple Measurements of differend types"<<std::endl;
          continue;
        }
      }



      // get the fitter infos ------------------------------------------------------------------
      if (! tp->hasFitterInfo(rep)) {
        std::cerr<<"trackPoint has no fitterInfo for rep"<<std::endl;
        continue;
      } 

      genfit::AbsFitterInfo* fitterInfo = tp->getFitterInfo(rep);

      fi = dynamic_cast<genfit::KalmanFitterInfo*>(fitterInfo);
      if(fi == NULL) {
        std::cerr<<"can only display KalmanFitterInfo"<<std::endl;
        continue;
      }
      if (! fi->hasPredictionsAndUpdates()) {
        std::cerr<<"KalmanFitterInfo does not have all predictions and updates"<<std::endl;
        //continue;
      }
      else {
        try {
          fittedState = &(fi->getFittedState(true));
        }
        catch (genfit::Exception& e) {
          std::cerr << e.what();
          std::cerr<<"can not get fitted state"<<std::endl;
          fittedState = NULL;
          prevFi = fi;
          prevFittedState = fittedState;
          continue;
        }
      }

      if (fittedState == NULL) {
        if (fi->hasForwardUpdate()) {
          fittedState = fi->getForwardUpdate();
        }
        else if (fi->hasBackwardUpdate()) {
          fittedState = fi->getBackwardUpdate();
        }
        else if (fi->hasForwardPrediction()) {
          fittedState = fi->getForwardPrediction();
        }
        else if (fi->hasBackwardPrediction()) {
          fittedState = fi->getBackwardPrediction();
        }
      }

      if (fittedState == NULL) {
        std::cout << "cannot get any state from fitterInfo, continue.\n";
        prevFi = fi;
        prevFittedState = fittedState;
        continue;
      }

      TVector3 track_pos = fittedState->getPos();
			track_pos.Print();
		}
	}
*/	
		/*genfit::StateOnPlane reference(beamlinetrack->getCardinalRep());
		genfit::SharedPlanePtr firstPlane(beamlinetrack->getPointWithMeasurement(0)->getRawMeasurement(0)->constructPlane(reference));
		double residual = beamlinetrack->getFittedState(0).extrapolateToMeasurement(beamlinetrack->getPointWithMeasurement()->getRawMeasurement(),false,false);
    */

/*	
  //----------------------------------------------------------------------------------------------------
  //Finding residuals
  //----------------------------------------------------------------------------------------------------
    for(unsigned i=0; i<beamlinetrack->getNumPoints(); i++)
  	{
			if(beamlinetrack->getNumPoints() != nDetBeam) continue;
    	TVector3 track_at_plane = beamlinetrack->getFittedState(i).getPos();
    	double hit_x = beamlinetrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[0];
    	double hit_y = beamlinetrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[1];

    	// Residuals
    	double res_x = hit_x - track_at_plane.X();
    	double res_y = hit_y - track_at_plane.Y();

    	fx << res_x << "\n";
    	fy << res_y << "\n";
			if(i == GEM2){
				fx2 << res_x << "\n";
				fy2 << res_y << "\n";
			}
  	}//end residuals loop
	}//end loop for number of tracks

  fx.close(); //close the files storing the residuals
  fy.close();
  fx2.close(); 
  fy2.close();
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//End Modifications
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/* //Uncomment for original example
  //check
  assert(fitTrack.checkConsistency());

  // do the fit
  fitter->processTrack(&fitTrack);

  // print fit result
  fitTrack.getFittedState().Print();

  //check
  assert(fitTrack.checkConsistency());


  display->addEvent(&fitTrack);

  for (int i=0; i<nDetBeam; i++){
    detectorPositions.push_back(i*84); //Detectors are 84mm apart
  }

  delete fitter;
*/

  // open event display
  //display->open();
  return 0;
}


std::vector<double> getHits(std::string branch){
	//Function to get hit coordinates from root file
	
	std::vector<double> *gem = 0;
	std::vector<double> hits;
	
	TFile *fin = TFile::Open("default_run_action.root", "READ");
  TTree *t; fin->GetObject("T",t);

  t->SetBranchAddress(branch.c_str(), &gem);
  for (int n = 0; n<t->GetEntries(); ++n){
    t-> GetEntry(n);

		//If there is no hit recorded in the detector, mark that with noHitMarker
		if(gem->size() == 0) hits.push_back(noHitMarker);
		//Else take the first event, just the particle striking the GEM
		else hits.push_back(gem->at(0));
  }
  t->ResetBranchAddresses();

  return hits;
}

std::vector<double> getCoord(std::vector<double> v0, std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, std::vector<double> v4, int loc)
{
	//This function is designed to get the the x or y coordinates from each gem and put them in a vector
	std::vector<double> coords; //hold positions of hits
    
	for(int i=0; i<nDetBeam; i++)
  {
	  switch(i){
	    case 0: {
	      coords.push_back(v0.at(loc));
	    }
	    case 1: {
	      coords.push_back(v1.at(loc));
	    }
	    case 2: {
	      coords.push_back(v2.at(loc));
				}
	    case 3: {
	      coords.push_back(v3.at(loc));
	    }
	    case 4: {
	      coords.push_back(v4.at(loc));
	    }
	  }//end switch statement

  }//end loop 

	return coords;
}

std::vector<int> getPositions(){

	std::vector<int> posVec;

	posVec.push_back(0);
	posVec.push_back(15); //Distance between SiPM detectors is 20mm, I think. Found in g4PSIDetectorConstruction.cc line 1249
	posVec.push_back(55); //Distance between SiPM and GEM1 is 40mm I think, ASK ABOUT THIS
	posVec.push_back(139); //Distance between GEMs is 84mm
	posVec.push_back(223);

	return posVec;
}

void beamlinestuff(){

  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. ,10., 0.)); // 1 T


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack(); Uncomment for original


  // particle pdg code; pion hypothesis
  const int pdg = 211;

  // start values for the fit, e.g. from pattern recognition
  TVector3 pos(0, 0, 0);
  TVector3 mom(0, 0, 3);


  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

  // create track
  genfit::Track fitTrack(rep, pos, mom);

  //Open files to write residuals
  std::ofstream fx;
  std::ofstream fy;
  std::ofstream fx2;
  std::ofstream fy2;
  fx.open("GEM_xres.txt");
  fy.open("GEM_yres.txt");
  fx2.open("GEM2_xres.txt");
  fy2.open("GEM2_yres.txt");

  //Getting hit coordinates from simulation
	std::vector<double> gem1x,gem1y,gem2x,gem2y,gem3x,gem3y,sipm1x,sipm1y,sipm2x,sipm2y,scl1x,scl1y,scl1z,scl2x,scl2y,scl2z,scr1x,scr1y,scr1z,scr2x,scr2y,scr2z,stcl1x,stcl1y,stcl1z,stcl2x,stcl2y,stcl2z,stcr1x,stcr1y,stcr1z,stcr2x,stcr2y,stcr2z; 

	gem1x = getHits("GEM1_PosHitX");
	gem1y = getHits("GEM1_PosHitY");
	gem2x = getHits("GEM2_PosHitX");
	gem2y = getHits("GEM2_PosHitY");
	gem3x = getHits("GEM3_PosHitX");
	gem3y = getHits("GEM3_PosHitY");
	sipm1x = getHits("BCC1_PosHitX");
	sipm1y = getHits("BCC1_PosHitY");
	sipm2x = getHits("BCC2_PosHitX");
	sipm2y = getHits("BCC2_PosHitY");
	scl1x = getHits("SCL1_PosHitX");
	scl1y = getHits("SCL1_PosHitY");
	scl1z = getHits("SCL1_PosHitZ");
	scl2x = getHits("SCL2_PosHitX");
	scl2y = getHits("SCL2_PosHitY");
	scl2z = getHits("SCL2_PosHitZ");
	scr1x = getHits("SCR1_PosHitX");
	scr1y = getHits("SCR1_PosHitY");
	scr1z = getHits("SCR1_PosHitZ");
	scr2x = getHits("SCR2_PosHitX");
	scr2y = getHits("SCR2_PosHitY");
	scr2z = getHits("SCR2_PosHitZ");
	stcl1x = getHits("STCL1_PosHitX");
	stcl1y = getHits("STCL1_PosHitY");
	stcl1z = getHits("STCL1_PosHitZ");
	stcl2x = getHits("STCL2_PosHitX");
	stcl2y = getHits("STCL2_PosHitY");
	stcl2z = getHits("STCL2_PosHitZ");
	stcr1x = getHits("STCR1_PosHitX");
	stcr1y = getHits("STCR1_PosHitY");
	stcr1z = getHits("STCR1_PosHitZ");
	stcr2x = getHits("STCR2_PosHitX");
	stcr2y = getHits("STCR2_PosHitY");
	stcr2z = getHits("STCR2_PosHitZ");

		
	//Create track in GenFit
  genfit::Track* beamlinetrack; //Track before target
  //genfit::Track* scatteredtrack; //Track after hitting target
  //genfit::EventDisplay* display;

  genfit::MaterialEffects::getInstance()->setNoEffects();
  //Turning off material effects (for now)
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,0.));
  //No mag field
  display = genfit::EventDisplay::getInstance();
	
  TVectorD hitCoords(2); //will hold x and y coordinate of single hit
  const int detId(0); //detector ID
  int planeId(0); //detector plane ID
  int hitId(0); //hit ID
	double smearVal(0.1); //width of applied gaussian smear (mm) centered at 0
  double detectorResolution(1); //resolution of GEMs detectors, mm
	int numtracks = gem3x.size(); //Number of tracks to try to reconstruct
	int GEM2(3), GEM3(4); //For easier access to coordinates in vectors
  
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
	hitCov *= detectorResolution*detectorResolution;

  genfit::AbsKalmanFitter* fitter = new genfit::DAF(); //KalmanFitter(20,0.001);
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(20,0.001); //KalmanFitter(20,0.001);

  //creating z coordinates of GEMs
  std::vector<int> detectorPositions; //hold positions of detectors
	detectorPositions = getPositions();

  //Loop over number of tracks
  for(int k=0; k < numtracks; k++)
  {  
		//Declaring a new track
		//The two TVector3s are posSeed and momSeed. Can not figure out the difference if they are changed,
		//e.g. the two lines below (one of which should remain commented) both produce the same output.
		//More info in GenFit/core/src/track.cc 
    //beamlinetrack = new genfit::Track(new genfit::RKTrackRep(pdg), pos, mom);
    beamlinetrack = new genfit::Track(new genfit::RKTrackRep(11), TVector3(125,125,1700), TVector3(0,0,1));
    std::vector<double> xcoords; //hold positions of hits
    std::vector<double> ycoords;
    
		xcoords = getCoord(sipm1x, sipm2x, gem1x, gem2x, gem3x, k); //X coordinates for 1 particle across all detectors 
		ycoords = getCoord(sipm2y, sipm2y, gem1y, gem2y, gem3y, k); //Y
		
			//put stuff here if things break, add curly bracket after //end loop to find hits, comment below until for loop 
			int noHitPos = 0; //will hold index of first detector marked as having no hit
			for(int it=0; it < nDetBeam; it++)
			{ 
				//This loop is meant to find the first detector that records no hit, store its position in 'noHit'
				//in order to only fit to the detectors that record hits. 
				if(xcoords.at(it) == noHitMarker && it < GEM3){
					noHitPos = it;
					break;
				}
				else if(it == GEM3 && xcoords.at(it) == noHitMarker){
					noHitPos = GEM3;
					break;
				}
				else if(it == GEM3 && xcoords.at(it) != noHitMarker){
					noHitPos = nDetBeam;
					break;
				}
			}//end loop to find noHitPos
			if(noHitPos == 0) continue;

    	for (int i=0; i<noHitPos; i++)
  		{
    		//Add Gaussian smear to the hits
    	  hitCoords[0] = xcoords[i]*(1+gRandom->Gaus(0,smearVal));
    	  hitCoords[1] = ycoords[i]*(1+gRandom->Gaus(0,smearVal));
    	  int z = (detectorPositions[i]);
  
	  	  genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
	  		measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z), TVector3(1,0,0), TVector3(0,1,0))), planeId++);
	  		beamlinetrack->insertPoint(new genfit::TrackPoint(measurement, beamlinetrack));
  		}//end loop to find hits

    if(beamlinetrack->getNumPoints()==0){
    	std::cout<< "No Points in beamlinetrack\n";
    	return; //nothing to fit
  	}
		
    assert(beamlinetrack->checkConsistency()); //check
    fitter->processTrack(beamlinetrack); //do the fit
    assert(beamlinetrack->checkConsistency());
    display->addTrack(beamlinetrack);

//......................................................................................................
//.............................Trying to Extrapolate....................................................
//......................................................................................................	
	/*
	for(int poop=0; poop<nDetBeam; poop++){
    unsigned int numhits = beamlinetrack->getNumPointsWithMeasurement();
		//genfit::TrackPoint* tp = beamlinetrack->getPointWithMeasurement(poop); 

    genfit::KalmanFitterInfo* fi;
    genfit::KalmanFitterInfo* prevFi = 0;
    const genfit::MeasuredStateOnPlane* fittedState(NULL);
    const genfit::MeasuredStateOnPlane* prevFittedState(NULL);

    for(unsigned int j = 0; j < numhits; j++) { // loop over all hits in the track

      fittedState = NULL;

      genfit::TrackPoint* tp = beamlinetrack->getPointWithMeasurement(j);
      if (! tp->hasRawMeasurements()) {
        std::cerr<<"trackPoint has no raw measurements"<<std::endl;
        continue;
      }

      const genfit::AbsMeasurement* m = tp->getRawMeasurement();
      //int hit_coords_dim = m->getDim();

      // check if multiple AbsMeasurements are of same type
      if (tp->getNumRawMeasurements() > 1) {
        bool sameTypes(true);
        for (unsigned int iM=1; iM<tp->getNumRawMeasurements(); ++iM) {
          if (typeid(*(tp->getRawMeasurement(iM))) != typeid(*m))
            sameTypes = false;
        }
        if (!sameTypes) {
          std::cerr<<"cannot draw trackpoint containing multiple Measurements of differend types"<<std::endl;
          continue;
        }
      }



      // get the fitter infos ------------------------------------------------------------------
      if (! tp->hasFitterInfo(rep)) {
        std::cerr<<"trackPoint has no fitterInfo for rep"<<std::endl;
        continue;
      } 

      genfit::AbsFitterInfo* fitterInfo = tp->getFitterInfo(rep);

      fi = dynamic_cast<genfit::KalmanFitterInfo*>(fitterInfo);
      if(fi == NULL) {
        std::cerr<<"can only display KalmanFitterInfo"<<std::endl;
        continue;
      }
      if (! fi->hasPredictionsAndUpdates()) {
        std::cerr<<"KalmanFitterInfo does not have all predictions and updates"<<std::endl;
        //continue;
      }
      else {
        try {
          fittedState = &(fi->getFittedState(true));
        }
        catch (genfit::Exception& e) {
          std::cerr << e.what();
          std::cerr<<"can not get fitted state"<<std::endl;
          fittedState = NULL;
          prevFi = fi;
          prevFittedState = fittedState;
          continue;
        }
      }

      if (fittedState == NULL) {
        if (fi->hasForwardUpdate()) {
          fittedState = fi->getForwardUpdate();
        }
        else if (fi->hasBackwardUpdate()) {
          fittedState = fi->getBackwardUpdate();
        }
        else if (fi->hasForwardPrediction()) {
          fittedState = fi->getForwardPrediction();
        }
        else if (fi->hasBackwardPrediction()) {
          fittedState = fi->getBackwardPrediction();
        }
      }

      if (fittedState == NULL) {
        std::cout << "cannot get any state from fitterInfo, continue.\n";
        prevFi = fi;
        prevFittedState = fittedState;
        continue;
      }

      TVector3 track_pos = fittedState->getPos();
			track_pos.Print();
		}
	}
*/	
		/*genfit::StateOnPlane reference(beamlinetrack->getCardinalRep());
		genfit::SharedPlanePtr firstPlane(beamlinetrack->getPointWithMeasurement(0)->getRawMeasurement(0)->constructPlane(reference));
		double residual = beamlinetrack->getFittedState(0).extrapolateToMeasurement(beamlinetrack->getPointWithMeasurement()->getRawMeasurement(),false,false);
    */

	
  //----------------------------------------------------------------------------------------------------
  //Finding residuals
  //----------------------------------------------------------------------------------------------------
    for(unsigned i=0; i<beamlinetrack->getNumPoints(); i++)
  	{
			if(beamlinetrack->getNumPoints() != nDetBeam) continue;
    	TVector3 track_at_plane = beamlinetrack->getFittedState(i).getPos();
    	double hit_x = beamlinetrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[0];
    	double hit_y = beamlinetrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[1];

    	// Residuals
    	double res_x = hit_x - track_at_plane.X();
    	double res_y = hit_y - track_at_plane.Y();

    	fx << res_x << "\n";
    	fy << res_y << "\n";
			if(i == GEM2){
				fx2 << res_x << "\n";
				fy2 << res_y << "\n";
			}
  	}//end residuals loop
	}//end loop for number of tracks

  fx.close(); //close the files storing the residuals
  fy.close();
  fx2.close(); 
  fy2.close();

  display->open();
	return; 
}
