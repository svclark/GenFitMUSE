/*
This code takes information from the g4PSI simulation and attempts to reconstruct tracks using GenFit.
Written by Steven Clark. svclark96@gmail.com
Summer 2017
*/
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
#include <math.h>

const double noHitMarker = 99999; //Marker given to identify when particle does not hit detector 
const int nDetBeam = 5; //number of detectors before the target. 2 SiPMs, 3 GEMs
const int nDetScatter = 4; //number of detectors after the target. 2 STTs, 2 scintillators on each side

//Function Declarations
std::vector<double> getHits(std::string branch);
double getEvent(std::string branch, int number);
//std::vector<double> getBeamCoord(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, int);
//std::vector<double> getScatterCoord(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, int);
std::vector<double> getZPositions();
void beamlinestuff(genfit::EventDisplay*);
void scatterstuff(genfit::EventDisplay*, bool);

int main() {
	
  // init geometry and mag. field
  //new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");

	// init event display
	genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

	//Do fitting on beamline detectors
	beamlinestuff(display);
	
	//Do fitting on scattered particle detectors
	scatterstuff(display, false); //false does right side detectors, true does left side
	scatterstuff(display, true);
	
	display->open();

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
		//Else take the first event
		else hits.push_back(gem->at(0) / 10);
  }
  t->ResetBranchAddresses();
	gem->clear();
  return hits;
}

double getEvent(std::string branch, int number){

	std::vector<double> *gem = 0;
	double coord = 0;
	
	TFile *fin = TFile::Open("default_run_action.root", "READ");
  TTree *t; fin->GetObject("T",t);

  t->SetBranchAddress(branch.c_str(), &gem);
  
	if(t->GetEntry(number) != 0){
		if(gem->size() == 0) coord = 99999;
		else coord = gem->at(0);
	}
	else coord = 11111;

  t->ResetBranchAddresses();

	fin->Close();


  return coord;
}


/*
std::vector<double> getBeamCoord(std::vector<double> v0, std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, std::vector<double> v4, int loc)
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

std::vector<double> getScatterCoord(std::vector<double> v0, std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, int loc)
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
	  }//end switch statement

  }//end loop 

	return coords;
}
*/

std::vector<double> getZPositions(){
	//This function creates a vector with the z coordinates of the beamline detectors. Coordinates are from
	//the g4PSI simulation, which are output in geant4's native units (mm)
	std::vector<double> posVec;
	int offset = 0; //To line up with geometry

	posVec.push_back(-550.038/10 + offset); //Hardcoded in from museGeometry.gdml, line 2870ish
	posVec.push_back(-530.038/10 + offset);
	posVec.push_back(-468.038/10 + offset);
	posVec.push_back(-384.038/10 + offset);
	posVec.push_back(-300.038/10 + offset);

	return posVec;
}


void beamlinestuff(genfit::EventDisplay* display){

  // particle pdg code; pion hypothesis
  const int pdg = 11;

  // start values for the fit, e.g. from pattern recognition
  TVector3 pos(1, 2, 0);
  TVector3 mom(1, 2, 3);


  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

  //Open files to write residuals
  std::ofstream fx, fy, fx2, fy2, fveto, fhmveto;
  fx.open("GEM_xres.txt");
  fy.open("GEM_yres.txt");
  fx2.open("GEM2_xres.txt");
  fy2.open("GEM2_yres.txt");
	fveto.open("Veto_Residuals.txt");
	fhmveto.open("Veto_HitMissResiduals.txt");


	/*
  //Getting hit coordinates from simulation
	std::vector<double> gem1x,gem1y,gem2x,gem2y,gem3x,gem3y,sipm1x,sipm1y,sipm2x,sipm2y, vetox, vetoy; 

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
	vetox = getHits("VSC_PosHitX");
	vetoy = getHits("VSC_PosHitY");

	std::cout << "CHECK: " << gem1x.size() << " " << sipm2y.size() << std::endl;
	
	

	//Find actual number of veto hits
	//int vetoHitCount = 0;
	//for(unsigned it=0; it<vetox.size(); it++){
		//if(vetox.at(it) != noHitMarker && vetoy.at(it) != noHitMarker) ++vetoHitCount;
	//}
	std::vector<double> vrad;
	for(unsigned it=0; it<vetox.size(); it++){
		if(vetox.at(it) != noHitMarker && vetoy.at(it) != noHitMarker){
			vrad.push_back(std::sqrt(vetox.at(it)*vetox.at(it) + vetoy.at(it)*vetoy.at(it)));
		}
		else vrad.push_back(0);
	}
	//Loop below finds minimum r value, used to get an idea of how small it can go
	//double vrmin=5;
	//for(unsigned it=0; it<vrad.size(); it++){
		//if(vrad.at(it) < 3) std::cout << it << " " << vrad.at(it) << std::endl;
		//if(vrad.at(it) < vrmin) vrmin = vrad.at(it);
	//}
	//Now I've run this I few times, I will manually type in r
	//TODO: Find a better way to determine if there are veto hits
	double vrmin = 2.87;
	//std::cout << "Veto Hit Count: " << vetoHitCount << std::endl;
	//int vrCount = 0;
	//for(unsigned it=0; it<vrad.size(); it++){
		//if(vrad.at(it) > vrmin) ++vrCount;
	//}
	//std::cout << "Veto Hits, According to Radius: " << vrCount << std::endl;
	*/
		
	//Create track in GenFit
  genfit::Track* beamlinetrack; //Track before target

  //genfit::MaterialEffects::getInstance()->setNoEffects();
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  //Turning off material effects (for now)
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,0.));
  //No mag field
	
  TVectorD hitCoords(2); //will hold x and y coordinate of single hit
  const int detId(0); //detector ID
  int planeId(0); //detector plane ID
  int hitId(0); //hit ID
	double smearVal(0.1); //width of applied gaussian smear (mm) centered at 0
  double detectorResolution(1); //resolution of GEMs detectors, mm
	//int numtracks = gem3x.size(); //Number of tracks to try to reconstruct
	int GEM3(4); //For easier access to coordinates in vectors
	unsigned GEM2(3);
  
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
	hitCov *= detectorResolution*detectorResolution;

  //genfit::AbsKalmanFitter* fitter = new genfit::DAF(); //KalmanFitter(20,0.001);
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(); //KalmanFitter(20,0.001);
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(20,0.001); //KalmanFitter(20,0.001);
	int predictionCount = 0;
	int hmCount = 0;
	int whereDidYouGo = 0;

  //creating z coordinates of GEMs
  std::vector<double> detectorPositions; //hold positions of detectors
	detectorPositions = getZPositions();
  //Loop over number of tracks
	int k=0;
	std::cout << "Starting big beamline loop\n";
	if(k == -1) return;
  //for(int k=0; k < numtracks; k++){
	while(k != -1){
		//Declaring a new track
    //beamlinetrack = new genfit::Track(new genfit::RKTrackRep(11), TVector3(125,125,1700), TVector3(0,0,1));
    beamlinetrack = new genfit::Track(rep, pos, mom);
    //std::vector<double> xcoords; //hold positions of hits
    //std::vector<double> ycoords;
		/*
		double xcoords [5] = {sipm1x.at(0), sipm2x.at(0), gem1x.at(0), gem2x.at(0), gem3x.at(0)};
		double ycoords [5] = {sipm1y.at(0), sipm2y.at(0), gem1y.at(0), gem2y.at(0), gem3y.at(0)};
		sipm1x.erase(sipm1x.begin());
		sipm2x.erase(sipm2x.begin());
		gem1x.erase(gem1x.begin());
		gem2x.erase(gem2x.begin());
		gem3x.erase(gem3x.begin());
		sipm1y.erase(sipm1y.begin());
		sipm2y.erase(sipm2y.begin());
		gem1y.erase(gem1y.begin());
		gem2y.erase(gem2y.begin());
		gem3y.erase(gem3y.begin());
		*/
		double sipm1x(0), sipm1y(0), sipm2x(0), sipm2y(0), gem1x(0), gem1y(0), gem2x(0), gem2y(0), gem3x(0), gem3y(0);
		sipm1x = getEvent("BCC1_PosHitX", k);
		if(sipm1x == 11111){
			k=-1;
			return;
		}
		sipm1y = getEvent("BCC1_PosHitY", k);
		sipm2x = getEvent("BCC2_PosHitX", k);
		sipm2y = getEvent("BCC2_PosHitY", k);
		gem1x = getEvent("GEM1_PosHitX", k);
		gem1y = getEvent("GEM1_PosHitY", k);
		gem2x = getEvent("GEM2_PosHitX", k);
		gem2y = getEvent("GEM2_PosHitY", k);
		gem3x = getEvent("GEM3_PosHitX", k);
		gem3y = getEvent("GEM3_PosHitY", k);


		double xcoords[5] = {sipm1x, sipm2x, gem1x, gem2x, gem3x};
		double ycoords[5] = {sipm1y, sipm2y, gem1y, gem2y, gem3y};
    
		//xcoords = getBeamCoord(sipm1x, sipm2x, gem1x, gem2x, gem3x, k); //X coordinates for 1 particle across all detectors 
		//ycoords = getBeamCoord(sipm2y, sipm2y, gem1y, gem2y, gem3y, k); //Y
		
			//put stuff here if things break, add curly bracket after //end loop to find hits, comment below until for loop 
			int noHitPos = 0; //will hold index of first detector marked as having no hit
			for(int it=0; it < nDetBeam; it++)
			{ 
				//This loop is meant to find the first detector that records no hit, store its position in 'noHitPos'
				//in order to only fit to the detectors that record hits. 
				if(xcoords[it] == noHitMarker && it < GEM3){
					noHitPos = it;
					break;
				}
				else if(it == GEM3 && xcoords[it] == noHitMarker){
					noHitPos = GEM3;
					break;
				}
				else if(it == GEM3 && xcoords[it] != noHitMarker){
					noHitPos = nDetBeam;
					break;
				}
			}//end loop to find noHitPos
			if(noHitPos == 0 || noHitPos == 1) continue; //Need to have hit on at least 2 detectors to fit

    	for (int i=0; i<noHitPos; i++)
  		{
    		//Add Gaussian smear to the hits
    	  hitCoords[0] = xcoords[i]*(1+gRandom->Gaus(0,smearVal));
    	  hitCoords[1] = ycoords[i]*(1+gRandom->Gaus(0,smearVal));
    	  double z = (detectorPositions[i]);
  
	  	  genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
	  		measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z), TVector3(1,0,0), TVector3(0,1,0))), planeId++);
	  		beamlinetrack->insertPoint(new genfit::TrackPoint(measurement, beamlinetrack));
  		}//end loop to find hits

			try{
    	assert(beamlinetrack->checkConsistency()); //check
    	fitter->processTrack(beamlinetrack); //do the fit
    	assert(beamlinetrack->checkConsistency());
    	display->addTrack(beamlinetrack);
			}
			catch(const std::exception& e){
				std::cout << "Exception: " << e.what() << std::endl;
				continue;
			}
	

  //----------------------------------------------------------------------------------------------------
  //Extrapolation
  //----------------------------------------------------------------------------------------------------

			//only extrapolate if particle strikes GEM3
			/*
			double vetoz = -250.038/10;
			if(beamlinetrack->getNumPoints() == nDetBeam){
				genfit::MeasuredStateOnPlane stateRef(rep);
				genfit::StateOnPlane stateRefOrig(stateRef);
				genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,0,vetoz), TVector3(0,0,1)));
      	genfit::TrackPoint* tp = beamlinetrack->getPointWithMeasurementAndFitterInfo(4, rep);
      	if (tp == NULL) {
        	std::cout << "Track has no TrackPoint with fitterInfo! \n";
        	continue;
      	}
      	genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));

      	// extrapolate back to reference plane.
      	try{
        	//rep->extrapolateToPlane(kfsop, stateRefOrig.getPlane());
        	rep->extrapolateToPlane(kfsop, plane);
      	}
      	catch(genfit::Exception& e){
        	std::cerr<<"Exception, next track"<<std::endl;
        	std::cerr << e.what();
        	continue;
				}

				TVector3 prediction = kfsop.getPos();
				double residx(0), residy(0), prad(0), hmresidx(0), hmresidy(0);
				prad = std::sqrt(prediction.X()*prediction.X() + prediction.Y()*prediction.Y());
				if(vrad.at(k) > vrmin && prad >= vrmin){
					predictionCount++;
					residx = vetox.at(k) - prediction.X();
					residy = vetoy.at(k) - prediction.Y();
					fveto << residx << "\n";
					fveto << residy << "\n";	
				}
				//Find residuals for particles that hit in simulation, miss in extrapolation
				if(vrad.at(k) > vrmin && prad < vrmin){
					hmCount++;
					hmresidx = vetox.at(k) - prediction.X();
					hmresidy = vetoy.at(k) - prediction.Y();
					fhmveto << hmresidx << "\n";
					fhmveto << hmresidy << "\n";
				}
			}//end extrapolation
			else{
				if(vrad.at(k) > vrmin) ++whereDidYouGo;
			}
			*/
	
	
  //----------------------------------------------------------------------------------------------------
  //Finding residuals
  //----------------------------------------------------------------------------------------------------
		
    for(unsigned i=0; i<beamlinetrack->getNumPoints(); i++)
  	{
			//if(beamlinetrack->getNumPoints() != nDetBeam) continue;
			try{
			if(beamlinetrack->getNumPoints() < 2) continue;
      }
      catch(const std::exception& e){
        std::cout << "Exception: " << e.what() << std::endl;
        continue;
      }

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
		std::cout << "END\n";
		++k;
	}//end loop for number of tracks

  fx.close(); //close the files storing the residuals
  fy.close();
  fx2.close(); 
  fy2.close();
	fveto.close();
	fhmveto.close();
	//std::cout << "Predicted? " <<  predictionCount << std::endl;
	//std::cout << "Hit Sim, Miss GenFit: " << hmCount << std::endl;
	//std::cout << "Others? " << whereDidYouGo << std::endl;
	//std::cout << "PLEASE BE THE SAME NUMBER: " << vrCount << " " << predictionCount+hmCount+whereDidYouGo << std::endl;


	return; 
}


void scatterstuff(genfit::EventDisplay* display, bool lor){
  // particle pdg code; pion hypothesis
  const int pdg = 211;

  // start values for the fit, e.g. from pattern recognition
  TVector3 pos(0, 0, 0);
  TVector3 mom(0, 0, 3);

  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

  // create track
  genfit::Track fitTrack(rep, pos, mom);

	//Files that will hold residuals
  std::ofstream fxs;
  std::ofstream fys;
  fxs.open("Scatter_xres.txt");
  fys.open("Scatter_yres.txt");

  //Getting hit coordinates from simulation
	//std::vector<double> sc1x,sc1y,sc1z,sc2x,sc2y,sc2z,stc1x,stc1y,stc1z,stc2x,stc2y,stc2z;

	/*
	if(lor == false){
		sc1x = getHits("SCR1_PosHitX");
		sc1y = getHits("SCR1_PosHitY");
		sc1z = getHits("SCR1_PosHitZ");
		sc2x = getHits("SCR2_PosHitX");
		sc2y = getHits("SCR2_PosHitY");
		sc2z = getHits("SCR2_PosHitZ");
		stc1x = getHits("STCR1_PosHitX");
		stc1y = getHits("STCR1_PosHitY");
		stc1z = getHits("STCR1_PosHitZ");
		stc2x = getHits("STCR2_PosHitX");
		stc2y = getHits("STCR2_PosHitY");
		stc2z = getHits("STCR2_PosHitZ");
	}

	else if(lor == true){
		sc1x = getHits("SCL1_PosHitX");
		sc1y = getHits("SCL1_PosHitY");
		sc1z = getHits("SCL1_PosHitZ");
		sc2x = getHits("SCL2_PosHitX");
		sc2y = getHits("SCL2_PosHitY");
		sc2z = getHits("SCL2_PosHitZ");
		stc1x = getHits("STCL1_PosHitX");
		stc1y = getHits("STCL1_PosHitY");
		stc1z = getHits("STCL1_PosHitZ");
		stc2x = getHits("STCL2_PosHitX");
		stc2y = getHits("STCL2_PosHitY");
		stc2z = getHits("STCL2_PosHitZ");
	}

	else{ 
		std::cout << "WARNING: Need to call scatterstuff() with 0 for right or 1 for left\n";
		return; }
	std::cout << "Scatter Check: " << stc2y.size() << "\n\n\n\n";
	*/

	//Create track in GenFit 
	genfit::Track* scatteredtrack; //Track before target

  //Turning off material effects (for now)
  //genfit::MaterialEffects::getInstance()->setNoEffects();
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,0.));
  //No mag field
  display = genfit::EventDisplay::getInstance();
	
  TVectorD hitCoords(2); //will hold x and y coordinate of single hit
  const int detId(0); //detector ID
  int planeId(0); //detector plane ID
  int hitId(0); //hit ID
	double smearVal(0.1); //width of applied gaussian smear (mm) centered at 0
  double detectorResolution(1); //resolution of GEMs detectors, mm
	//int numtracks = stc1x.size(); //Number of tracks to try to reconstruct
  
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
	hitCov *= detectorResolution*detectorResolution;

  //genfit::AbsKalmanFitter* fitter = new genfit::DAF(); //KalmanFitter(20,0.001);
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(); //KalmanFitter(20,0.001);
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(20,0.001); //KalmanFitter(20,0.001);

  //creating z coordinates of GEMs
  std::vector<int> detectorPositions; //hold positions of detectors

  //Loop over number of tracks
	int k=0;
	std::cout << "Starting big scatter loop\n\n\n\n";
	if(k == -1) return;
  //for(int k=0; k < numtracks; k++){
  while(k != -1){  
		//Declaring a new track
    //scatteredtrack = new genfit::Track(new genfit::RKTrackRep(pdg), pos, mom);
    scatteredtrack = new genfit::Track(new genfit::RKTrackRep(11), TVector3(125,125,1700), TVector3(0,0,1));
    //std::vector<double> xcoords, ycoords, zcoords; //hold positions of hits
		/*
		double xcoords [5] = {stc1x.at(0), stc2x.at(0), sc1x.at(0), sc2x.at(0)};
		double ycoords [5] = {stc1y.at(0), stc2y.at(0), sc1y.at(0), sc2y.at(0)};
		double zcoords [5] = {stc1z.at(0), stc2z.at(0), sc1z.at(0), sc2z.at(0)};
		stc1x.erase(stc1x.begin());
		stc2x.erase(stc2x.begin());
		sc1x.erase(sc1x.begin());
		sc2x.erase(sc2x.begin());
		stc1y.erase(stc1y.begin());
		stc2y.erase(stc2y.begin());
		sc1y.erase(sc1y.begin());
		sc2y.erase(sc2y.begin());
		stc1z.erase(stc1z.begin());
		stc2z.erase(stc2z.begin());
		sc1z.erase(sc1z.begin());
		sc2z.erase(sc2z.begin());
		*/

		double scx1(0), scy1(0), scz1(0), scx2(0), scy2(0), scz2(0), stcx1(0), stcy1(0), stcz1(0), stcx2(0), stcy2(0), stcz2(0);
		if(lor == true){
			scx1 = getEvent("SCL1_PosHitX", k);
			scy1 = getEvent("SCL1_PosHitY", k);
			scz1 = getEvent("SCL1_PosHitZ", k);
			scx2 = getEvent("SCL2_PosHitX", k);
			scy2 = getEvent("SCL2_PosHitY", k);
			scz2 = getEvent("SCL2_PosHitZ", k);
			stcx1 = getEvent("STCL1_PosHitX", k);
			stcy1 = getEvent("STCL1_PosHitY", k);
			stcz1 = getEvent("STCL1_PosHitZ", k);
			stcx2 = getEvent("STCL2_PosHitX", k);
			stcy2 = getEvent("STCL2_PosHitY", k);
			stcz2 = getEvent("STCL2_PosHitZ", k);
		}
		if(lor == false){
			scx1 = getEvent("SCR1_PosHitX", k);
			std::cout << scx1 << std::endl;
			scy1 = getEvent("SCR1_PosHitY", k);
			scz1 = getEvent("SCR1_PosHitZ", k);
			scx2 = getEvent("SCR2_PosHitX", k);
			scy2 = getEvent("SCR2_PosHitY", k);
			scz2 = getEvent("SCR2_PosHitZ", k);
			stcx1 = getEvent("STCR1_PosHitX", k);
			stcy1 = getEvent("STCR1_PosHitY", k);
			stcz1 = getEvent("STCR1_PosHitZ", k);
			stcx2 = getEvent("STCR2_PosHitX", k);
			stcy2 = getEvent("STCR2_PosHitY", k);
			stcz2 = getEvent("STCR2_PosHitZ", k);
		}

		if(stcx1 == 11111) {
			k=-1;
			return;
		}

		double xcoords [4] = {scx1, scx2, stcx1, stcx2};
		double ycoords [4] = {scy1, scy2, stcy1, stcy2};
		double zcoords [4] = {scz1, scz2, stcz1, stcz2};
    
		//xcoords = getScatterCoord(stc1x, stc2x, sc1x, sc2x, k); //x coordinates for 1 particle across all detectors 
		//ycoords = getScatterCoord(stc1y, stc2y, sc1y, sc2y, k); //y coordinates for 1 particle across all detectors 
		//zcoords = getScatterCoord(stc1z, stc2z, sc1z, sc2z, k); //z coordinates for 1 particle across all detectors 
		
			int noHitPos = 0; //will hold index of first detector marked as having no hit
			for(int it=0; it < nDetScatter; it++)
			{ 
				if(xcoords[it] == noHitMarker && it < nDetScatter-1){
					noHitPos = it;
					break;
				}
				else if(it == nDetScatter-1 && xcoords[it] == noHitMarker){
					noHitPos = nDetScatter-1;
					break;
				}
				else if(it == nDetScatter-1 && xcoords[it] != noHitMarker){
					noHitPos = nDetScatter;
					break;
				}
			}//end loop to find noHitPos
			if(noHitPos == 0 || noHitPos == 1) continue; //Need hit on at least 2 detectors to get fit
			//if(noHitPos == 0) continue;

    	for (int i=0; i<noHitPos; i++)
  		{
    		//Add Gaussian smear to the hits
    	  hitCoords[0] = xcoords[i]*(1+gRandom->Gaus(0,smearVal));
    	  hitCoords[1] = ycoords[i]*(1+gRandom->Gaus(0,smearVal));
    	  double z = zcoords[i];
  
	  	  genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
	  		measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z), TVector3(1,0,0), TVector3(0,1,0))), planeId++);
	  		scatteredtrack->insertPoint(new genfit::TrackPoint(measurement, scatteredtrack));

  		}//end loop to find hits

    if(scatteredtrack->getNumPoints()==0){
    	std::cout<< "No Points in scatteredtrack\n";
    	return; //nothing to fit
  	}
		
  	assert(scatteredtrack->checkConsistency()); //check
  	fitter->processTrack(scatteredtrack); //do the fit
  	assert(scatteredtrack->checkConsistency());
  	display->addTrack(scatteredtrack);


  //----------------------------------------------------------------------------------------------------
  //Finding residuals
  //----------------------------------------------------------------------------------------------------
		
    for(unsigned i=0; i<scatteredtrack->getNumPoints(); i++)
  	{
			try{
				if(scatteredtrack->getNumPoints() != nDetScatter) continue;
      }
      catch(const std::exception& e){
        std::cout << "Exception: " << e.what() << std::endl;
        continue;
      }
			
			double hit_x(0), hit_y(0), res_x(0), res_y(0);
			try{
    	TVector3 track_at_plane = scatteredtrack->getFittedState(i).getPos();
    	hit_x = scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[0];
    	hit_y = scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[1];

    	// Residuals
    	res_x = hit_x - track_at_plane.X();
    	res_y = hit_y - track_at_plane.Y();
      }
      catch(const std::exception& e){
        std::cout << "Exception: " << e.what() << std::endl;
        continue;
      }

    	fxs << res_x << "\n";
    	fys << res_y << "\n";
  	}//end residuals loop
		++k;
	}//end loop for number of tracks
	

  fxs.close(); //close the files storing the residuals
  fys.close();

	return;
}

