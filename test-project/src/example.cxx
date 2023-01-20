#include "PROPOSAL/PROPOSAL.h"
#include <iostream>

using namespace PROPOSAL;

int main(){
    Propagator prop(MuMinusDef::Get(), "resources/config_iceaxel.json");
    Particle& mu = prop.GetParticle();
    Particle mu_backup(mu);

    mu_backup.SetEnergy(9e6);
    mu_backup.SetDirection(Vector3D(0, 0, -1));

    std::vector<double> ranges;

    for (int i = 0; i <10; i++)
    {
    mu.InjectState(mu_backup);
    
    std::vector<DynamicData*> dd = prop.Propagate();
    std::cout << "Doing muon " << i << std::endl;
    for(DynamicData* d : dd) {std::cout << d->GetNameFromType(d->GetTypeId()) << std::endl;}
    
    ranges.push_back(mu.GetPropagatedDistance());
    }
    
    // ... Do stuff with ranges, e.g. plot histogram
    // for(auto i:ranges) std::cout << i << std::endl;
}
