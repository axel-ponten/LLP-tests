#include "PROPOSAL/PROPOSAL.h"
#include <iostream>
#include <string>
#include <map>

using namespace PROPOSAL;


void printDynamicData(DynamicData* D);
void printDynamicDataVector(std::vector<DynamicData*> DD);
bool inDetector(Vector3D pos);
void saveDataFrameToCSV(std::map<std::string, std::vector<double>> df, std::string filename);

int main(){

    Output::getInstance().EnableASCIIOutput("test");
    Propagator prop(MuMinusDef::Get(), "/home/axel/i3/testing-PROPOSAL/minimal-example-standalone-cpp/resources/config_iceaxel.json");
    Particle& mu = prop.GetParticle();
    Particle mu_backup(mu); // to reset muon before new propagation

    mu_backup.SetEnergy(9e6); // [MeV]
    mu_backup.SetDirection(Vector3D(0, 0, -1));
    mu_backup.SetPosition(Vector3D(0,0,190000)); // [cm]

    std::map<std::string, std::vector<double>> LLP_frame; // "dataframe" for event information
    std::string filename = "LLP_dataframe_2.csv";

    const int n_LLP = 10;
    int i = 0;
    while (i < n_LLP)
    {
        mu.InjectState(mu_backup);
        std::vector<DynamicData*> dd = prop.Propagate();
        DynamicData* LLP_DynamicData = nullptr;
        for(DynamicData* d : dd) {
            // if LLP interaction
            if(d->GetNameFromType(d->GetTypeId()) == "LLPInt" && inDetector(d->GetPosition())) {
                i++; // number of LLP interactions in detector
                LLP_DynamicData = d; // set LLP DynamicData temp holder

                // add event information to dataframe
                LLP_frame["Muon_propagated_distance"].push_back((d->GetPosition() - mu_backup.GetPosition()).magnitude() / 100.); // [m]
                LLP_frame["LLP_parent_energy"].push_back(d->GetParentParticleEnergy() / 1000.); // [GeV]
                LLP_frame["LLP_energy"].push_back(d->GetEnergy() / 1000.); // [GeV]
                LLP_frame["fractional_energy"].push_back(d->GetEnergy() / d->GetParentParticleEnergy());
                LLP_frame["LLP_propagated_distance"].push_back(LLP_DynamicData->GetPropagatedDistance() / 100.); // [m]
                LLP_frame["LLP_x"].push_back(d->GetPosition().GetX() / 100.); // [m]
                LLP_frame["LLP_y"].push_back(d->GetPosition().GetY() / 100.); // [m]
                LLP_frame["LLP_z"].push_back(d->GetPosition().GetZ() / 100.); // [m]

                std::cout << "\n\n ### NEW EVENT ### \n";
                printDynamicDataVector(dd);
                std::cout << "\nEvent characteristics: Muon traveled = " <<  (d->GetPosition() - mu_backup.GetPosition()).magnitude()
                << " cm and x = " << d->GetEnergy() / d->GetParentParticleEnergy()
                << "and decay length = " << d->GetPropagatedDistance() << "\n";

                // iterate through the daughter particles
                std::vector<Particle> LLP_daughters;
                for(DynamicData* D : dd) {
                    if (D->GetNameFromType(D->GetTypeId()) == "Particle") {
                        if(dynamic_cast<Particle*>(D)->GetName() == "MuPlus" || dynamic_cast<Particle*>(D)->GetName() == "MuMinus") {
                            Particle p(*dynamic_cast<Particle*>(D));
                            LLP_daughters.push_back(p);
                        }
                    }
                }
                for(Particle daughter : LLP_daughters) {
                    mu.InjectState(daughter);
                    std::cout << "Propagating " << daughter.GetName() <<"\n";
                    printDynamicData(&daughter);
                    std::vector<DynamicData*> dd_LLP_daughter = prop.Propagate();
                    printDynamicDataVector(dd_LLP_daughter);
                }
            }
        }

    }
    saveDataFrameToCSV(LLP_frame, filename);
}


// helper functions
void printDynamicData(DynamicData* D) {
    if(D->GetNameFromType(D->GetTypeId()) == "Particle") {
        std::cout << dynamic_cast<Particle*>(D)->GetName() << ", ";
    } else std::cout << D->GetNameFromType(D->GetTypeId()) << ", ";
    std::cout << "Parent Energy: " << D->GetParentParticleEnergy() << ", "
    << "Energy taken: " << D->GetEnergy() << ", "
    << "Position: " << D->GetPosition().GetX() << ", " << D->GetPosition().GetY() << ", " << D->GetPosition().GetZ() << ", "
    << "Time: " << D->GetTime() << ", "
    << "Propagated Distance: " << D->GetPropagatedDistance()
    <<  "\n";
}

void printDynamicDataVector(std::vector<DynamicData*> DD) {
    double previous_energy = DD.at(0)->GetParentParticleEnergy();
    for(DynamicData* D : DD) {
        printDynamicData(D);
        if(D->GetNameFromType(D->GetTypeId()) != "ContinuousEnergyLoss" && std::abs(previous_energy - D->GetParentParticleEnergy()) > 100.)
        {
            //std::cout << "OOPS! previous parentenergy - previous energy loss doesn't add up to current parentenergy\n";
        }
        
        if(D->GetNameFromType(D->GetTypeId()) == "ContinuousEnergyLoss") {
            previous_energy -= D->GetEnergy();
        } else
            previous_energy = D->GetParentParticleEnergy() - D->GetEnergy();
    }
}

bool inDetector(Vector3D pos) {
    double height = 80000.; // 800 m
    double radius = 80000.;

    if(pos.GetZ() < -height) {
        return false;
    } else {
        return true;
    }
}

void saveDataFrameToCSV(std::map<std::string, std::vector<double>> df, std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    // create vector of the frame keys
    std::vector<std::string> df_keys;
    for(const auto& pair : df) {
        df_keys.push_back(pair.first);
    }
    // print header
    for(int j = 0; j < df_keys.size(); j++) {
        if (j == 0) outfile << df_keys.at(j);
        else outfile << "," << df_keys.at(j);
    }
    outfile << "\n";
    // print content of vectors
    for(int i = 0; i < df[df_keys.at(0)].size(); i++) {
        for(int j = 0; j < df_keys.size(); j++) {
            if (j == 0) outfile << df[df_keys.at(j)].at(i);
            else outfile << "," << df[df_keys.at(j)].at(i);
        }
        outfile << "\n";
    }
    outfile.close();
}