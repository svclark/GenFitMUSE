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

const double noHitMarker = -99999; //Marker given to identify when particle does not hit detector 
const int nDetBeam = 5; //number of detectors before the target. 2 SiPMs, 3 GEMs
const int nDetScatter = 4; //number of detectors after the target. 2 STTs, 2 scintillators on each side
const int electronPDG = 11; //PDG code to identify electrons
std::vector<int> scatterSpots; //Holds event number of scattered particles for using between functions

//Function Declarations
std::vector<double> getHits(std::string root, std::string extension, bool extra);
std::vector<double> getZPositions();
void beamlinestuff(genfit::EventDisplay*);
void scatterstuff(genfit::EventDisplay*, bool);
bool findInScatterSpots(int match);

int main() {
	
  // init geometry and mag. field
  TGeoManager::Import("genfitGeom.root");

	// init event display
	genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

	//Do fitting on scattered particle detectors
	//scatterstuff(display, false); //false does stage left detectors, true does left side 
	//scatterstuff(display, true); 

	//std::cout << "Number of Scatters: " << scatterSpots.size() << std::endl;

	//Do fitting on beamline detectors
	beamlinestuff(display);
	
	display->open();

	return 0;
}


std::vector<double> getHits(std::string root, std::string extension, bool extra){
	/*
	Function to read in data from root file. called with root of the branch, for example "GEM1" and 
	also the extension you want to read, like "_PosHitX". The bool should be set to false for any branch
	that has particle ID data you wish to read
	*/
	
	std::vector<double> *gem = 0;
	std::vector<double> hits, edep;
	std::vector<int> id;
	std::string ParticleID = "_ParticleID"; 
	std::string branch = (std::string) root + (std::string) extension; 
	std::string checkID = (std::string) root + (std::string) ParticleID;
	std::string Edep = "_Edep";
	std::string checkEdep = (std::string) root + (std::string) Edep;

	
	TFile *fin = TFile::Open("1k_default_run_action.root", "READ"); //Open the root file
  TTree *t; fin->GetObject("T",t); //Name of tree in file

	if(!extra){ //If the the input bool parameter is false, go read in the Particle ID data into vector 'id'
  	t->SetBranchAddress(checkID.c_str(), &gem);
  	for (int n = 0; n<t->GetEntries(); ++n){
    	t-> GetEntry(n);
			if(gem->size() == 0) id.push_back(noHitMarker);
			else id.push_back(gem->at(0));
  	}
  	t->ResetBranchAddresses();
		gem->clear();
	}

	if(root == "SCR1" || root == "SCL1"){ //if it's a scintillator, put energy deposit into vector
  	t->SetBranchAddress(checkEdep.c_str(), &gem);
  	for (int n = 0; n<t->GetEntries(); ++n){
    	t-> GetEntry(n);
			if(gem->size() == 0) id.push_back(noHitMarker);
			else {
				 edep.push_back(gem->at(0));
			}
  	}
  	t->ResetBranchAddresses();
		gem->clear();
		}

  t->SetBranchAddress(branch.c_str(), &gem);
  for (int n = 0; n<t->GetEntries(); ++n){
    t-> GetEntry(n);

		//If there is no hit recorded in the detector, mark that with noHitMarker
		if(gem->size() == 0) hits.push_back(noHitMarker);
		else if(extra == false && id.at(n) != electronPDG){
			hits.push_back(noHitMarker);
		}
		//If it's a scintillator and the energy deposit is less than 5, mark with noHitMarker
		else if((root == "SCR1" || root =="SCL1") && n < edep.size() && edep.at(n) != noHitMarker && edep.at(n) < 5){
			hits.push_back(noHitMarker);
		}
		//Else take the first event
		else hits.push_back(gem->at(0) / 10);
  }
  t->ResetBranchAddresses();
	gem->clear();
	id.clear();
	edep.clear();
  return hits;
}


std::vector<double> getZPositions(){
	//This function creates a vector with the z coordinates of the beamline detectors. Coordinates are from
	//the g4PSI simulation, which are output in geant4's native units (mm)
	std::vector<double> posVec;

	posVec.push_back(-550.038/10); //Hardcoded in from museGeometry.gdml, line 2870ish
	posVec.push_back(-530.038/10);
	posVec.push_back(-468.038/10);
	posVec.push_back(-384.038/10);
	posVec.push_back(-300.038/10);

	return posVec;
}

bool findInScatterSpots(int match){
	/*
	This function meant to find if the value sent to it exists anywhere in global vector scatterSpots
	Returns true if the value is found
	*/
	bool found = false;
	for(unsigned i=0; i<scatterSpots.size(); ++i){
		if(scatterSpots.at(i) == match){
			found = true;
			break;
		}
		else continue;
	}

	return found;
}

//Big function to do all beamline analysis
void beamlinestuff(genfit::EventDisplay* display){

  // particle pdg code 
  const int pdg = electronPDG;

  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
		
	//Create track in GenFit
  genfit::Track* beamlinetrack; //Track before target

  //genfit::MaterialEffects::getInstance()->setNoEffects();
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  //Turning off material effects (for now)
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,0.));
  //No mag field

  //genfit::AbsKalmanFitter* fitter = new genfit::DAF(); //KalmanFitter(20,0.001);
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(); //KalmanFitter(20,0.001);
 

  //Open files to write residuals output data
  std::ofstream fx, fy, fx2, fy2, fveto, fhmveto, bcx, bcy, bcz, bxcoords, bycoords;
  fx.open("./MyOutputFiles/GEM_xres.txt");
  fy.open("./MyOutputFiles/GEM_yres.txt");
  fx2.open("./MyOutputFiles/GEM2_xres.txt");
  fy2.open("./MyOutputFiles/GEM2_yres.txt");
	fveto.open("./MyOutputFiles/Veto_Residuals.txt");
	fhmveto.open("./MyOutputFiles/Veto_HitMissResiduals.txt");
	bcx.open("./MyOutputFiles/BeamCenterX.txt");
	bcy.open("./MyOutputFiles/BeamCenterY.txt");
	bcz.open("./MyOutputFiles/BeamCenterZ.txt");
	bxcoords.open("./MyOutputFiles/BeamlineFitXCoords.txt");
	bycoords.open("./MyOutputFiles/BeamlineFitYCoords.txt");

  //Getting hit coordinates from simulation
	std::vector<double> gem1x,gem1y,gem2x,gem2y,gem3x,gem3y,sipm1x,sipm1y,sipm2x,sipm2y,vetox,vetoy,x0,y0,z0,p0;

	gem1x = getHits("GEM1","_PosHitX",0);
	gem1y = getHits("GEM1","_PosHitY",0);
	gem2x = getHits("GEM2","_PosHitX",0);
	gem2y = getHits("GEM2","_PosHitY",0);
	gem3x = getHits("GEM3","_PosHitX",0);
	gem3y = getHits("GEM3","_PosHitY",0);
	sipm1x = getHits("BCC1","_PosHitX",0);
	sipm1y = getHits("BCC1","_PosHitY",0);
	sipm2x = getHits("BCC2","_PosHitX",0);
	sipm2y = getHits("BCC2","_PosHitY",0);
	vetox = getHits("VSC","_PosHitX",0);
	vetoy = getHits("VSC","_PosHitY",0);
	x0 = getHits("beam","_x0",1); //The last parameter is true because there is no _ParticleID branch
	y0 = getHits("beam","_y0",1);
	z0 = getHits("beam","_z0",1);
	p0 = getHits("beam","_p",1);

  TVectorD hitCoords(2); //will hold x and y coordinate of single hit
  const int detId(0); //detector ID
  int planeId(0); //detector plane ID
  int hitId(0); //hit ID
  double detectorResolution(0.2); //resolution of detectors, mm
	int numtracks = gem3x.size(); //Number of tracks to try to reconstruct
	int GEM3(nDetBeam-1), GEM1(nDetBeam-3); //For easier access to coordinates in vectors
	unsigned GEM2(nDetBeam-2); //unsigned because I didn't like the annoying warn it gave me when int
	unsigned spot = 0; //Used in extrapolation to center
	int vetoHitCount = 0; //Will count how many particles hit the veto detector
	double vrmin = 2.87; //Used to find when a paricle hits the veto. Approximated the center of the veto as a circle with this radius. This was determined by calculating the radius (sqrt(x^2 + y^2)) of particles that hit the veto for a large sample, and then looking for the minimum 
	double vrad = 0; //Will hold a single radius for comparison to vrmin
  
  TMatrixDSym hitCov(2); //Making a covariance matrix
  hitCov.UnitMatrix();
	hitCov *= detectorResolution*detectorResolution;

  //These variables all hold info about veto hits and are printed eventually, mostly for debugging
	int predictionCount = 0;
	int hmCount = 0;
	int whereDidYouGo = 0;
	int cCount = 0;

  //creating z coordinates of GEMs
  std::vector<double> detectorPositions; //hold positions of detectors
	detectorPositions = getZPositions();
  //Loop over number of tracks
	std::cout << "Starting big beamline loop"  << std::endl;
  for(int k=0; k < numtracks; k++)
  {  
		//Getting the first value in vectors, putting that in xcoords, then removing that element from the array. This was an attempt to save some memory. Could change 0 to k in all instances below and remove lines 334-347
		double xcoords [5] = {sipm1x.at(0), sipm2x.at(0), gem1x.at(0), gem2x.at(0), gem3x.at(0)};
		double ycoords [5] = {sipm1y.at(0), sipm2y.at(0), gem1y.at(0), gem2y.at(0), gem3y.at(0)};
		//double xcoords [3] = {gem1x.at(0), gem2x.at(0), gem3x.at(0)}; if only looking at GEMs, use this
		//double ycoords [3] = {gem1y.at(0), gem2y.at(0), gem3y.at(0)};
		double pos [3] = {x0.at(0), y0.at(0), z0.at(0)};
		double mom [3] = {0, 0, p0.at(0)};
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
		x0.erase(x0.begin());
		y0.erase(y0.begin());
		z0.erase(z0.begin());
		p0.erase(p0.begin());

		//Finding radius at veto detector
		if(vetox.at(k) != noHitMarker && vetoy.at(k) != noHitMarker){
			vrad = (std::sqrt(vetox.at(k)*vetox.at(k) + vetoy.at(k)*vetoy.at(k)));
		}
    
		//declaring a new track
		beamlinetrack = new genfit::Track(rep, pos, mom);
    
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
			
			if(noHitPos == GEM1 || noHitPos == GEM2 || noHitPos == GEM3) continue; //Need to have hit on last 3 GEMs to fit 
			/*To change the criteria for fitting to hitting just a certain number of detectors instead of the big loop above, comment that loop and 380, uncomment the line below and 381 
			*/
			//if(xcoords[0] == noHitMarker || xcoords[1] == noHitMarker) continue;

    	for (int i=0; i<noHitPos; i++)
    	//for (int i=0; i<nDetBeam; i++)
  		{
    		//Add Gaussian smear to the hits
				double rand1(0),rand2(0);
				rand1 = gRandom->Gaus(0,0.1);
				rand2 = gRandom->Gaus(0,0.1);
    	  hitCoords[0] = xcoords[i]*(1+rand1);
    	  hitCoords[1] = ycoords[i]*(1+rand2);
    	  double z = (detectorPositions[i]);
  
				//Here is important GenFit stuff
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
  //Finding residuals
  //----------------------------------------------------------------------------------------------------
		
    for(unsigned i=0; i<beamlinetrack->getNumPoints(); i++)
  	{
			try{
			if(beamlinetrack->getNumPoints() < 2) continue; //only find residuals for tracks that hit more than 2 detectors
      }
      catch(const std::exception& e){
        std::cout << "Exception: " << e.what() << std::endl;
        continue;
      }

			try{
				//Get fit coordinates and actual hit coordinates
				TVector3 track_at_plane = beamlinetrack->getFittedState(i).getPos();
    		double hit_x = beamlinetrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[0];
    		double hit_y = beamlinetrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[1];

    		// Residuals
    		double res_x = hit_x - track_at_plane.X();
    		double res_y = hit_y - track_at_plane.Y();

    		fx << res_x << " " << hit_x << " " << track_at_plane.X() << "\n";
    		fy << res_y << "\n";
				bxcoords << track_at_plane.X() << "\n";
				bycoords << track_at_plane.Y() << "\n";
				if(i == GEM2){
					fx2 << res_x << "\n";
					fy2 << res_y << "\n";
				}
      }
      catch(const std::exception& e){
        std::cout << "Exception: " << e.what() << std::endl;
        continue;
      }
  	}//end residuals loop

  //----------------------------------------------------------------------------------------------------
  //Extrapolation to Veto
  //----------------------------------------------------------------------------------------------------

			double vetoz = -250.038/10; //Found from geometry file, plane where the veto is
			if(beamlinetrack->getNumPoints() == nDetBeam){ //Extrapolate if there is a hit in all bl detectors
				genfit::MeasuredStateOnPlane stateRef(rep);
				genfit::StateOnPlane stateRefOrig(stateRef);
				//Plane to extrapolate to
				genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,0,vetoz), TVector3(0,0,1)));
				//below tells the function to get the fitter info at that point at GEM3
      	genfit::TrackPoint* tp = beamlinetrack->getPointWithMeasurementAndFitterInfo(GEM3, rep);
      	if (tp == NULL) {
        	std::cout << "Track has no TrackPoint with fitterInfo! \n";
        	continue;
      	}
      	genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));

      	// extrapolate back to reference plane.
      	try{
        	rep->extrapolateToPlane(kfsop, plane);
      	}
      	catch(genfit::Exception& e){
        	std::cerr<<"Exception, next track"<<std::endl;
        	std::cerr << e.what();
        	continue;
				}

				TVector3 prediction = kfsop.getPos();
				double residx(0), residy(0), prad(0), hmresidx(0), hmresidy(0);
				//Get predicted radius
				prad = std::sqrt(prediction.X()*prediction.X() + prediction.Y()*prediction.Y());
				//If there is a predicted hit and an actual hit, do some naalysis on that
				if(vrad > vrmin && prad >= vrmin && vetox.at(k) != noHitMarker && vetoy.at(k) != noHitMarker){
					++vetoHitCount;
					++predictionCount;
					residx = vetox.at(k) - prediction.X();
					residy = vetoy.at(k) - prediction.Y();
					fveto << residx << "\n";
					fveto << residy << "\n";	

  				TVectorD extrapCoords(2); //will hold x and y coordinate of single hit
					extrapCoords[0] = prediction.X();
					extrapCoords[1] = prediction.Y();
					//Used for display purposes, not essential
	  			genfit::PlanarMeasurement* extrap = new genfit::PlanarMeasurement(extrapCoords, hitCov, detId, ++hitId, NULL);
	  			extrap->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,vetoz), TVector3(1,0,0), TVector3(0,1,0))), planeId++);

					//Display extrapolation, kind of cheating but this works
	  			beamlinetrack->insertPoint(new genfit::TrackPoint(extrap, beamlinetrack));
    			assert(beamlinetrack->checkConsistency()); //check
    			fitter->processTrack(beamlinetrack); //do the fit
    			assert(beamlinetrack->checkConsistency());
    			display->addTrack(beamlinetrack);

				}
				//Find residuals for particles that hit in simulation, miss in extrapolation
				if(vrad > vrmin && prad < vrmin && vetox.at(k) != noHitMarker && vetoy.at(k) != noHitMarker){
					++vetoHitCount;
					++hmCount;
					hmresidx = vetox.at(k) - prediction.X();
					hmresidy = vetoy.at(k) - prediction.Y();
					fhmveto << hmresidx << "\n";
					fhmveto << hmresidy << "\n";
				}
			}//end extrapolation
			else{
				if(vrad > vrmin) ++whereDidYouGo;
			}
				

  //----------------------------------------------------------------------------------------------------
  //Extrapolation to Center
  //----------------------------------------------------------------------------------------------------
			//Very similar to veto extrapolation, but there is less analysis
			if(findInScatterSpots(k)){
					++spot;
					++cCount;

					genfit::MeasuredStateOnPlane stateRef(rep);
					rep->setPosMom(stateRef, pos, mom);
					genfit::StateOnPlane stateRefOrig(stateRef);
					genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,0,0), TVector3(0,0,1)));
      		genfit::TrackPoint* tp = beamlinetrack->getPointWithMeasurementAndFitterInfo(nDetBeam-1, rep);

      		if (tp == NULL) {
        		std::cout << "Track has no TrackPoint with fitterInfo! \n";
        		continue;
      		}
      		genfit::KalmanFittedStateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate())); //Not sure of difference between Backward and Forward

      		//extrapolate back to reference plane.
      		try{
        		rep->extrapolateToPlane(kfsop, plane);
      		}
      		catch(genfit::Exception& e){
        		std::cerr<<"Exception, next track"<<std::endl;
        		std::cerr << e.what();
        		continue;
					}
					TVector3 vertexPrediction = kfsop.getPos();
					TVectorD extrapCoords(2);
					bcx << k << " " << vertexPrediction.X() << "\n";
					bcy << k << " " << vertexPrediction.Y() << "\n";
					bcz << k << " " << vertexPrediction.Z() << "\n";
					extrapCoords[0] = vertexPrediction.X();
					extrapCoords[1] = vertexPrediction.Y();
	  			genfit::PlanarMeasurement* extrap = new genfit::PlanarMeasurement(extrapCoords, hitCov, detId, ++hitId, NULL);
	  			extrap->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,0), TVector3(1,0,0), TVector3(0,1,0))), planeId++);

					//Display extrapolation, kind of cheating but this works
	  			beamlinetrack->insertPoint(new genfit::TrackPoint(extrap, beamlinetrack));
    			assert(beamlinetrack->checkConsistency()); //check
    			fitter->processTrack(beamlinetrack); //do the fit
    			assert(beamlinetrack->checkConsistency());
    			display->addTrack(beamlinetrack);
			}
			
	
	}//end loop for number of tracks

	//You must sing this when this code is executing: http://www.lyricsmode.com/lyrics/b/barney/clean_up.html
	gem1x.clear();
	gem1y.clear();
	gem2x.clear();
	gem2y.clear();
	gem3x.clear();
	gem3y.clear();
	sipm1x.clear();
	sipm1y.clear();
	sipm2x.clear();
	sipm2y.clear();
	vetox.clear();
	vetoy.clear();
	x0.clear();
	y0.clear();
	z0.clear();
	p0.clear();
	detectorPositions.clear();

  fx.close(); //close the files storing the residuals
  fy.close();
  fx2.close(); 
  fy2.close();
	fveto.close();
	fhmveto.close();
	bcx.close();
	bcy.close();
	bcx.close();
	bxcoords.close();
	bycoords.close();
	std::cout << "Beam Center Extrapolation: " << cCount << std::endl;
	std::cout << "Veto Hit Count: " << vetoHitCount << std::endl;
	std::cout << "Correctly Predicted: " <<  predictionCount << std::endl;
	std::cout << "Hit Sim, Miss GenFit: " << hmCount << std::endl;
	std::cout << "Others? " << whereDidYouGo << std::endl;


	return; 
}


//Now do everything again for the scattered detectors. This works very similar, I will try to comment differences
void scatterstuff(genfit::EventDisplay* display, bool lor){
  // particle pdg code; pion hypothesis
  const int pdg = electronPDG; 

  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
	//Create track in GenFit 
	genfit::Track* scatteredtrack; //Track before target

  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,0.));
  //No mag field
  display = genfit::EventDisplay::getInstance();

  //genfit::AbsKalmanFitter* fitter = new genfit::DAF(); //KalmanFitter(20,0.001);
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(); //KalmanFitter(20,0.001);

	//Files that will hold residuals
  std::ofstream fxs, fys, sxcoords, sycoords, szcoords;
 	std::ofstream fscx, fscy, fscz;
	fxs.open("./MyOutputFiles/Scatter_xres.txt");
  fys.open("./MyOutputFiles/Scatter_yres.txt");
	sxcoords.open("./MyOutputFiles/ScatterXCoords.txt");
	sycoords.open("./MyOutputFiles/ScatterYCoords.txt");
	szcoords.open("./MyOutputFiles/ScatterZCoords.txt");
	if(!lor){ //Only open if calling scatterstuff with false. this is why false must be called first
		fscx.open("./MyOutputFiles/ScatterCenterX.txt");
		fscy.open("./MyOutputFiles/ScatterCenterY.txt");
		fscz.open("./MyOutputFiles/ScatterCenterZ.txt");
	}
	if(lor){
		fscx.open("./MyOutputFiles/ScatterCenterX.txt", fstream::out | fstream::app);
		fscy.open("./MyOutputFiles/ScatterCenterY.txt", fstream::out | fstream::app);
		fscz.open("./MyOutputFiles/ScatterCenterZ.txt", fstream::out | fstream::app);
	}

  //Getting hit coordinates from simulation
	std::vector<double> sc1x,sc1y,sc1z,sc2x,sc2y,sc2z,stc1x,stc1y,stc1z,stc2x,stc2y,stc2z,x0,y0,z0,p0;

	if(lor == false){
		sc1x = getHits("SCR1","_PosHitX",0);
		sc1y = getHits("SCR1","_PosHitY",0);
		sc1z = getHits("SCR1","_PosHitZ",0);
		sc2x = getHits("SCR2","_PosHitX",0);
		sc2y = getHits("SCR2","_PosHitY",0);
		sc2z = getHits("SCR2","_PosHitZ",0);
		stc1x = getHits("STCR1","_PosHitX",0);
		stc1y = getHits("STCR1","_PosHitY",0);
		stc1z = getHits("STCR1","_PosHitZ",0);
		stc2x = getHits("STCR2","_PosHitX",0);
		stc2y = getHits("STCR2","_PosHitY",0);
		stc2z = getHits("STCR2","_PosHitZ",0);
	}

	else if(lor == true){
		sc1x = getHits("SCL1","_PosHitX",0);
		sc1y = getHits("SCL1","_PosHitY",0);
		sc1z = getHits("SCL1","_PosHitZ",0);
		sc2x = getHits("SCL2","_PosHitX",0);
		sc2y = getHits("SCL2","_PosHitY",0);
		sc2z = getHits("SCL2","_PosHitZ",0);
		stc1x = getHits("STCL1","_PosHitX",0);
		stc1y = getHits("STCL1","_PosHitY",0);
		stc1z = getHits("STCL1","_PosHitZ",0);
		stc2x = getHits("STCL2","_PosHitX",0);
		stc2y = getHits("STCL2","_PosHitY",0);
		stc2z = getHits("STCL2","_PosHitZ",0);
	}
	
	else{ 
		std::cout << "WARNING: Need to call scatterstuff() with 0 for right or 1 for left\n";
		return; }
	std::cout << "Scatter Check: " << stc2y.size() << "\n";

	x0 = getHits("beam","_x0",1);
	y0 = getHits("beam","_y0",1);
	z0 = getHits("beam","_z0",1);
	p0 = getHits("beam","_p",1);
	
  TVectorD hitCoords(2); //will hold x and y coordinate of single hit
  const int detId(0); //detector ID
  int planeId(0); //detector plane ID
  int hitId(0); //hit ID
  double detectorResolution(0.1); //resolution of GEMs detectors, mm
	int numtracks = stc1x.size(); //Number of tracks to try to reconstruct
  
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
	hitCov *= detectorResolution*detectorResolution;

  //Loop over number of tracks
	std::cout << "Starting big scatter loop\n";
  for(int k=0; k < numtracks; k++)
  {  

		//Declaring a new track
		double xcoords [4] = {stc1x.at(0), stc2x.at(0), sc1x.at(0), sc2x.at(0)};
		double ycoords [4] = {stc1y.at(0), stc2y.at(0), sc1y.at(0), sc2y.at(0)};
		double zcoords [4] = {stc1z.at(0), stc2z.at(0), sc1z.at(0), sc2z.at(0)};
  	// start values for the fit, e.g. from pattern recognition
		double pos[3] = {0,0,0};
		double mom[3] = {0,0,115};
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
		x0.erase(x0.begin());
		y0.erase(y0.begin());
		z0.erase(z0.begin());
		p0.erase(p0.begin());

    scatteredtrack = new genfit::Track(rep, pos, mom);
		
			/*
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
			*/
			//Doing cuts if particle doesn't hit all detectors, or if there is backward scattering or if xcoords are suspiciously close together
			double diff12 = std::abs(std::abs(xcoords[1])-std::abs(xcoords[0]));
			double diff23 = std::abs(std::abs(xcoords[2])-std::abs(xcoords[1]));
			double diff34 = std::abs(std::abs(xcoords[3])-std::abs(xcoords[2]));
			if(xcoords[0] == noHitMarker || xcoords[1] == noHitMarker || xcoords[2] == noHitMarker || xcoords[3] == noHitMarker) continue;
			//else if(std::abs(xcoords[0]) > 250 || std::abs(xcoords[3]) > 50) continue;
			else if(diff12 > 10 || diff23 > 10 || diff34 > 10){
				std::cout << "Cut" << std::endl;
				continue;
			}

    	for (int i=0; i<nDetScatter; i++)
  		{
    		//Add Gaussian smear to the hits
    	  hitCoords[0] = xcoords[i]*(1+gRandom->Gaus(0,0.1));
    	  hitCoords[1] = ycoords[i]*(1+gRandom->Gaus(0,0.1));
				
    	  double x = xcoords[i];
    	  double y = ycoords[i];
    	  double z = zcoords[i];
				double xlim = 1;
				double zlim = 1.7;
				if(lor) zlim *= -1;
				//double norm = std::sqrt(xlim*xlim + ylim*ylim + zlim*zlim);
				double norm = 1;

/*  
	  	  genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
				//Here is the line where the plane is set. It is troublesome
	  		//measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z), TVector3(1,0,0), TVector3(0,1,0))), planeId++);
	  		measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(x,y,z), TVector3(xlim,0,zlim), TVector3(0,1,0))), planeId++);
	  		scatteredtrack->insertPoint(new genfit::TrackPoint(measurement, scatteredtrack));
*/
				
				const genfit::HelixTrackModel* trackModel_ = nullptr; 
				double thetaDetPlane_ = 90; //TODO: Make these right
				double phiDetPlane_ = 0;
				int measurementCounter_ = 0;
				double tracklength = std::sqrt(x*x + y*y + z*z);
				int type = 0; //GenFit measurement type, 0=Pixel detector

				TVector3 dir(0,0,1);
				TVector3 point(x,y,z);
				//trackModel_->getPosDir(tracklength, point, dir);

  			TVector3 planeNorm(dir);
  			planeNorm.SetTheta(thetaDetPlane_*TMath::Pi()/180);
  			planeNorm.SetPhi(planeNorm.Phi()+phiDetPlane_);
  			static const TVector3 zdir(0,0,1);
  			static const TVector3 xdir(1,0,0);

      	genfit::SharedPlanePtr plane(new genfit::DetPlane(point, planeNorm.Cross(zdir), (planeNorm.Cross(zdir)).Cross(planeNorm)));
      genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, int(type), measurementCounter_, nullptr);
      static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, measurementCounter_);
	  	scatteredtrack->insertPoint(new genfit::TrackPoint(measurement, scatteredtrack));

  		}//end loop to find hits

    if(scatteredtrack->getNumPoints()==0){
    	std::cout<< "No Points in scatteredtrack\n";
    	return; //nothing to fit
  	}
		
		try{
  		assert(scatteredtrack->checkConsistency()); //check
  		fitter->processTrack(scatteredtrack); //do the fit
  		assert(scatteredtrack->checkConsistency());
  		display->addTrack(scatteredtrack);
		}
		catch(const std::exception& e){
			std::cout << "Exception: " << e.what() << std::endl;
			continue;
		}

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
			sxcoords << hit_x << "\n";
			sycoords << hit_y << "\n";
  	}//end residuals loop

/*
  //----------------------------------------------------------------------------------------------------
  //Extrapolation to Center
  //----------------------------------------------------------------------------------------------------
			
			if(xcoords[0] != noHitMarker){
				
				genfit::StateOnPlane state(rep);
				rep->setPosMom(state, pos, mom);
				genfit::SharedPlanePtr origPlane = state.getPlane();
				genfit::StateOnPlane origState(state);
				TVector3 point(0,0,0);
				try{
					rep->extrapolateToPoint(state, point);
				}
				catch(const std::exception& e){
					std::cout << "Exception: " << e.what() << std::endl;
					continue;
				}
				TVector3 vertexPrediction = state.getPos();
				
				TVectorD extrapCoords(2);
				fscx << k << " " << vertexPrediction.X() << "\n";
				fscy << k << " " << vertexPrediction.Y() << "\n";
				fscz << k << " " << vertexPrediction.Z() << "\n";
				extrapCoords[0] = vertexPrediction.X();
				extrapCoords[1] = vertexPrediction.Y();
  			genfit::PlanarMeasurement* extrap = new genfit::PlanarMeasurement(extrapCoords, hitCov, detId, ++hitId, NULL);
  			extrap->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,0), TVector3(1,0,0), TVector3(0,1,0))), planeId++);

				
				//Display extrapolation, kind of cheating but this works
  			scatteredtrack->insertPoint(new genfit::TrackPoint(extrap, scatteredtrack));
   			assert(scatteredtrack->checkConsistency()); //check
   			fitter->processTrack(scatteredtrack); //do the fit
   			assert(scatteredtrack->checkConsistency());
   			display->addTrack(scatteredtrack);

				scatterSpots.push_back(k);
			}
		*/
	}//end loop for number of tracks

	//Spring Cleaning
	sc1x.clear();
	sc1y.clear();
	sc1z.clear();
	sc2x.clear();
	sc2y.clear();
	sc2z.clear();
	stc1x.clear();
	stc1y.clear();
	stc1z.clear();
	stc2x.clear();
	stc2y.clear();
	stc2z.clear();
	x0.clear();
	y0.clear();
	z0.clear();
	p0.clear();

  fxs.close(); //close the files storing the residuals
  fys.close();
	sxcoords.close();
	sycoords.close();
	szcoords.close();
	fscx.close();
	fscy.close();
	fscz.close();


	return;
}

