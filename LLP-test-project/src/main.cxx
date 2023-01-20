#include "PROPOSAL/PROPOSAL.h"
#include <iostream>
#include <fstream>


#define CHECK(cond) if( !(cond) )\
{printf( "assertion error line %d, file(%s)\n", \
__LINE__, __FILE__ );} \
else {printf( "pass check line %d\n", \
__LINE__);}


using namespace PROPOSAL;


void testLLPModelIndependent() {
    std::cout << "------ start " << __func__  << " ------" << std::endl;
    ParticleDef particle_def = MuMinusDef::Get();
    Ice medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    double mass = 0.3;
    double epsilon = 1.;
    double energy = 100.;

    LLPModelIndependent* LLP_A = new LLPModelIndependent(particle_def, medium, multiplier, mass, epsilon);


    delete LLP_A;
    std::cout << "------" << __func__  << " done ------" << std::endl;
}

void testCrossSections() {
    std::cout << "------ start " << __func__  << " ------" << std::endl;


    ParticleDef particle_def = MuMinusDef::Get();
    Ice medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    double mass = 0.3;
    double epsilon = 1;
    double energy = 1000.;
    
    LLPModelIndependent* LLP_A = new LLPModelIndependent(particle_def, medium, multiplier, mass, epsilon);
    LLPModelIndependent param_int(particle_def, medium, multiplier, mass, epsilon);
    LLPIntegral* Int_A        = new LLPIntegral(param_int);

    // InterpolationDef InterpolDef;
    // LLPInterpolant* Interpol_A = new LLPInterpolant(param_int, InterpolDef);

    // // check cross sections
    // CHECK(Interpol_A->CalculatedEdx(1.) == 0)
    // CHECK(Interpol_A->CalculatedEdxWithoutMultiplier(1.) == 0)
    // CHECK(Interpol_A->CalculatedE2dx(1.) == 0)

    CHECK(Int_A->CalculatedEdx(1.) == 0)
    CHECK(Int_A->CalculatedEdxWithoutMultiplier(1.) == 0)
    CHECK(Int_A->CalculatedE2dx(1.) == 0)


    double chi = LLP_A->chiIWW(energy);
    double cs = LLP_A->DifferentialCrossSection(100., 0.001);
    double cs1 = LLP_A->DifferentialCrossSection(100., 0.2);
    double cs2 = LLP_A->DifferentialCrossSection(100., 0.7);
    double cs3 = LLP_A->DifferentialCrossSection(100., 0.9);
    std::cout << chi << " " << cs << " " << cs1 << " " << cs2 << " " << cs3 << std::endl;


    delete LLP_A;
    delete Int_A;
    // delete Interpol_A;

    std::cout << "------" << __func__  << " done ------" << std::endl;
}

int testComparisonEqualParticle(){
    std::cout << "------ start " << __func__  << " ------" << std::endl;
    int test_val = 1;

    ParticleDef particle_def = MuMinusDef::Get();
    Ice medium;
    EnergyCutSettings ecuts;
    double multiplier = 1.;
    double mass = 1.;
    double epsilon = 1e-4;
    
    LLPInteraction* LLP_A = new LLPModelIndependent(particle_def, medium, multiplier, mass, epsilon);
    Parametrization* LLP_B = new LLPModelIndependent(particle_def, medium, multiplier, mass, epsilon);
    
    if(*LLP_A != *LLP_B) test_val = 0;
    CHECK(*LLP_A == *LLP_B)

    LLPModelIndependent param_int(particle_def, medium, multiplier, mass, epsilon);
    CHECK(param_int == *LLP_A)


    LLPIntegral* Int_A        = new LLPIntegral(param_int);
    CrossSectionIntegral* Int_B = new LLPIntegral(param_int);
    CHECK(*Int_A == *Int_B)

    InterpolationDef InterpolDef;

    LLPInterpolant* Interpol_A        = new LLPInterpolant(param_int, InterpolDef);
    CrossSectionInterpolant* Interpol_B = new LLPInterpolant(param_int, InterpolDef);
    CHECK(*Interpol_A == *Interpol_B)

    delete LLP_A;
    delete LLP_B;
    delete Int_A;
    delete Int_B;
    delete Interpol_A;
    delete Interpol_B;

    std::cout << "------" << __func__  << " done ------" << std::endl;
    return test_val;
}


int main(){
    testComparisonEqualParticle();
    testCrossSections();
    testLLPModelIndependent();
}