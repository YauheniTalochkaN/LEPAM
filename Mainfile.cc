#include "RunManager.hh"

int main(int argc, char *argv[])
{
	gRandom->SetSeed(time(NULL));

	if(argc<2)
	{
    	printf("Not indicated config file! \n");
		return -1;
    }

	RunManager* Run_Manager = RunManager::getInstance();
	Run_Manager->SetDetectorConstruction(new DetectorConstruction());
	Run_Manager->SetRunAction(new RunAction());
	Run_Manager->SetPrimaryGeneratorAction(new PrimaryGeneratorAction());
	Run_Manager->SetTrackPropagation(new TrackPropagation());	
	Run_Manager->SetPhysicsList(new PhysicsList());
	Run_Manager->SetTGeoManager(new TGeoManager("Geometry", "Geometry"));	
	
	Run_Manager->ReadParam(argv[1]);
	Run_Manager->Initialization();

	auto start = std::chrono::high_resolution_clock::now();
	Run_Manager->Launch();
	auto stop = std::chrono::high_resolution_clock::now();

	std::cout << "\nSpent time: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " s\n";

	Run_Manager->SaveGeom("Geometry.root");

	delete Run_Manager;

	return 0;
}