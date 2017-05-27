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

int main()
{
	new TGeoManager("Geometry", "Geane geometry");
  	TGeoManager::Import("genfitGeom.root");

	TFile *file= new TFile("run7099_output.root");
	//TCanvas *c1 = new TCanvas("c1", "Fit1");

    //genfit::MeasurementCreator measurementCreator;
	genfit::MaterialEffects::getInstance()->setNoEffects();                         // Turning off Material effects (for now)
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. ,0., 0.));  // No mag field
    genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

/////////////////////////ROOT TREE STUFF///////////////////////////
    TFile *f = TFile::Open("Event.root","RECREATE");

	std::vector<double> fitdist;
   	std::vector<double> residuals;
   	std::vector<int> plane;
	std::vector<int> straw_in_plane;
	std::vector<double> wireY;
	std::vector<double> wireZ;
	std::vector<double> trackY;
	std::vector<double> trackZ;
	std::vector<int> badevent;
   	Int_t numhits;
   	// Create a ROOT Tree

   	TTree *tree = new TTree("T","An example of a ROOT tree");
	tree->Branch("fitdist", "vector<double>", &fitdist);
	tree->Branch("Residuals","vector<double>",&residuals);
	tree->Branch("Plane","vector<int>",&plane);
	tree->Branch("straw","vector<int>",&straw_in_plane);
	tree->Branch("WireZ","vector<double>",&wireZ);
	tree->Branch("WireY","vector<double>",&wireY);
	tree->Branch("TrackZ","vector<double>",&trackZ);
	tree->Branch("TrackY","vector<double>",&trackY);
	tree->Branch("numhits",&numhits);

 /////////////////////////////////////////////////////////////
	gRandom->SetSeed(0);

  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::DAF();
  //fitter->setDebugLvl(1);
  genfit::Track* scatteredtrack;
  genfit::WireMeasurementNew* measurement;
  
  // main loop
  for (unsigned int iEvent=0; iEvent<1000; ++iEvent)
  {
  	int plane_comp[10] ={0,0,0,0,0,0,0,0,0,0};
  	int uniqueplanes = 0;
	//numhits = gRandom->Integer(6)+2;
  	numhits = 3;
	scatteredtrack = new genfit::Track(new genfit::RKTrackRep(11), TVector3(125,125,1700), TVector3(0,0,1));
    for(int i = 0; i< numhits; i++)
    {
    	int p = gRandom->Integer(10);
    	plane_comp[p] = 1;
    	plane.push_back(p);//Randomly choose which plane to have a hit
    	// plane.push_back(0);
    	// plane.push_back(1);
    	// plane.push_back(2);
    	// plane.push_back(3);
    	// plane.push_back(4);

    	straw_in_plane.push_back((int)gRandom->Gaus(8,1));
    	// straw_in_plane.push_back(7);
    	// straw_in_plane.push_back(6);
    	// straw_in_plane.push_back(7);
    	// straw_in_plane.push_back(6);
    	// straw_in_plane.push_back(7);
	}
	// while(straw_in_plane.size() <numhits)//eliminate multiple hits in same straw
	// {
	// 	int s = (int)gRandom->Gaus(8,0.5);//randomly choose which straw had a hit. Gaus as beam is gaussian
	// 	if(std::find(straw_in_plane.begin(),straw_in_plane.end(),s)==straw_in_plane.end())
	// 		straw_in_plane.push_back(s);
	// }

	for(int j =0;j<10;j++)
	{
		uniqueplanes = uniqueplanes + plane_comp[j];//how many unique planes exist
		//uniqueplanes = 3;
	}

	if(uniqueplanes > 2 && numhits > 2)//We want three unique planes to make a track and at least 3 hits.
	{
		for(int hitID = 0; hitID < numhits; hitID++) 
		{   
		    
		      
			std::pair<TVector3,TVector3> wire;
		  	//double tube_diameter = 10;
		    double distance_between_planes = 8.7;
		    double pitch=10.1; //10.1 mm spacing
		    double zoffset=0;
		      
		    double pos=straw_in_plane.at(hitID)*pitch+(plane.at(hitID) % 2)*pitch/2;
		    wire.first[0]=-125;
		    wire.first[1]=pos;
		    wire.first[2]=zoffset+plane.at(hitID)*distance_between_planes;
		    wire.second[0]=125;
		    wire.second[1]=pos;
		    wire.second[2]=zoffset+plane.at(hitID)*distance_between_planes;
		    wireY.push_back(wire.second[1]);
		    wireZ.push_back(wire.second[2]);

		    double tim = gRandom->Gaus(60,10);//+10;
		 	//double tim = 50;
		 	double distance = (0.05*tim);


		 	fitdist.push_back(distance);
		 	//args are dist, disterr, startpos, stop pos, STT side, hit ID, NULL ptr
			measurement = new genfit::WireMeasurementNew(distance,0.4,wire.first,wire.second,2, hitID+1, NULL);
			measurement->setLeftRightResolution(1);
			measurement->setMaxDistance(5);
			scatteredtrack->insertPoint(new genfit::TrackPoint(measurement, scatteredtrack));	
		}
	}

    if( numhits >2 && uniqueplanes>2 )
    {
    	//check
    	assert(scatteredtrack->checkConsistency());
    	// do the fit
    	fitter->processTrack(scatteredtrack);

    	//check
    	assert(scatteredtrack->checkConsistency());
	}
	if(!scatteredtrack->getFitStatus()->isFitConvergedPartially())
		badevent.push_back(iEvent);
	double temp = 0;
	TVector3 track_at_tube;
	if(numhits>2 && uniqueplanes>2)
	{
		for(int hitID = 0; hitID < numhits; hitID++) 
		{
			//std::cout << " converged fully: " << scatteredtrack->getFitStatus()->isFitConvergedFully() << std::endl;
			//if(!scatteredtrack->getFitStatus()->isFitConvergedPartially())
			//std::cout << " converged partially: " << scatteredtrack->getFitStatus()->isFitConvergedPartially() << std::endl;
			//std::cout << " Residual? " << fitter->MeasuredStateOnPlane()->getResidual(iEvent,true,true) << std::endl;
		 	track_at_tube = scatteredtrack->getFittedState(hitID).getPos();
		 	trackY.push_back(track_at_tube.Y());
		 	trackZ.push_back(track_at_tube.Z());
		 	temp = sqrt(pow(wireY.at(hitID)-track_at_tube.Y(),2)+pow(wireZ.at(hitID)-track_at_tube.Z(),2)) - fabs(fitdist.at(hitID));
		 	if(scatteredtrack->getFitStatus()->isFitConvergedPartially() || scatteredtrack->getFitStatus()->isFitConvergedFully())
		 	residuals.push_back(temp);
		}
	}
	if(iEvent<1000)
		display->addEvent(scatteredtrack);

    // if (!scatteredtrack->getFitStatus()->isFitConvergedPartially()  && uniqueplanes >2) {
    //   // add track to event display
    // 	std::cout << " Bad Event Number: " << iEvent << " Unique Planes: "<< uniqueplanes << std::endl;
    //   display->addEvent(scatteredtrack);

    // }

   

  	tree->Fill();
  	fitdist.clear();
  	plane.clear();
  	straw_in_plane.clear();
  	wireZ.clear();
  	wireY.clear();
  	trackY.clear();
  	trackZ.clear();
  	residuals.clear();
  }// end loop over events


  f->Write();
  delete f;
  delete fitter;


  display->open();
  return 0;
}