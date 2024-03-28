#include "TrackPropagation.hh"
#include "RunManager.hh"

TrackPropagation::TrackPropagation()
{
   NUM_THREADs = 1;
   omp_set_num_threads(NUM_THREADs);

   ListOfTracks.reserve(10000);
}

TrackPropagation::~TrackPropagation()
{
    ClearTrackList();
}

void TrackPropagation::LockMutex()
{
    MTX.lock();
}

void TrackPropagation::UnLockMutex()
{
    MTX.unlock();
}

void TrackPropagation::SetNumThreads(Int_t num)
{
    NUM_THREADs = num;
    omp_set_num_threads(NUM_THREADs);
    std::cout << "Thread number: " << NUM_THREADs << "\n";
}

Int_t TrackPropagation::GetNumThreads()
{
    return NUM_THREADs;
}

void TrackPropagation::ThreadInitialization()
{
    NewTrackList = new std::vector<Track*>[NUM_THREADs];
}

Track* TrackPropagation::GetTrack(Int_t index)
{
    if((index >= SizeofTrackList()) || (index < 0)) {std::cout << "TrackPropagation: Wrong track number! \n"; return nullptr;}
    auto it = ListOfTracks.begin(); 
    std::advance(it, index);
    return *it;
}

Int_t TrackPropagation::SizeofTrackList()
{
    return ListOfTracks.size();
}

void TrackPropagation::ClearTrackList()
{
    if(ListOfTracks.size() == 0) return;
    
    for (auto &it : ListOfTracks)
    {
        delete it;
    }

    ListOfTracks.clear();

    ClearNewTrackList();
}

void TrackPropagation::ClearNewTrackList()
{
    for (int i = 0; i < NUM_THREADs; i++)
    {
        if(NewTrackList[i].size() > 0)
        {
            for (auto &it : NewTrackList[i])
            {
                delete it;
            }
            NewTrackList[i].clear();
        }
    } 
}

void TrackPropagation::CopyTracksFromNewTrackList()
{
    for (int i = 0; i < NUM_THREADs; i++)
    {
        if(NewTrackList[i].size() > 0)
        {            
            for (auto &it : NewTrackList[i])
            {
                ListOfTracks.push_back(it);
                RunManager::getInstance()->GetPhysicsList()->SetTGeoTrack(it);
            }
            NewTrackList[i].clear();
        }
    }
}

void TrackPropagation::AddNewTrack(Track* track)
{
    int id = omp_get_thread_num();
    if((id < 0) || (id >= NUM_THREADs)) 
    {
        std::cout << "Error: Wrong thread number when a new track is being added! \n";
        return;
    }

    NewTrackList[id].push_back(track);
}

Int_t TrackPropagation::MoveTracks()
{      
    RunManager::getInstance()->GetPhysicsList()->ApplyMultiTrackProcesses(ListOfTracks);

    RunManager::getInstance()->GetPhysicsList()->ApplySingleTrackProcesses(ListOfTracks);

    RunManager::getInstance()->GetPhysicsList()->ParticleStep(ListOfTracks);

    CopyTracksFromNewTrackList();

    int num_tasks = 0;
    for (auto &it : ListOfTracks)
    {
        if(it->GetKilled() == false) num_tasks++;
    }

    return num_tasks;
}
