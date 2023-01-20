#include "PROPOSAL/PROPOSAL.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

using namespace PROPOSAL;

double mass = 0.3;
double epsilon = 1.;
double mf = 0.000511;
//double mf = 0.105;
#define PI 3.14159265

double chiIWW(double energy)
{
    // double tmax = mass*mass + mf*mf;
    double tmax = energy*energy;
    double tmin = mass*mass*mass*mass / (4*energy*energy);
    double ZNucl = 8; // TODO: charge of oxygen instead
    double ANucl = 16; // TODO: mass of oxygen instead
    double Mel = 0.000511;
    
    double aa = 111.*std::pow(ZNucl,-1./3)/Mel;
    double d = 0.164*std::pow(ANucl,-2./3);
    double ta = std::pow(1./aa,2.);
    double td = d;
    double fluxAnalytical = ZNucl*ZNucl*(-((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(std::log(ta + tmax) - log(td + tmax) - std::log(ta + tmin) + std::log(td + tmin))))/((ta-td)*(ta-td)*(ta-td))));

    return fluxAnalytical;

}


double dsdyIWWamplitude(double t, double v)
{
    double mi = 0.105;
    double deltaM2 = mass*mass - mf*mf - mi*mi;
    double term1 = -1/(2*t);
    double term2 = deltaM2/(2*t*t);
    double term3 = -deltaM2*(mf*mf + v*(mi*mi*v + deltaM2))/(3*v*t*t*t);
    return term1 + term2 + term3;
}


double DifferentialCrossSection(double energy, double v)
{
    // TODO which units? cm^2?
    // TODO implement dsigma/dy from dark leptonic scalar?
    // equation 50. from master thesis?
    
    double mi = 0.105;
    double alpha = 1./137.;
    double betaf = std::sqrt( 1 - mf*mf / (energy*energy) );
    
    double prefactor = alpha*alpha*epsilon*epsilon/(4*PI);
    
    double psimax = 0.1;
    double ttildemax = mass*mass + v*energy*energy*psimax*psimax + mf*mf*(1-v)/v - mi*mi*(1-v);
    double ttildemin = mass*mass + mf*mf*(1-v)/v - mi*mi*(1-v);
    double amplitude_upper = dsdyIWWamplitude(ttildemax, v);
    double amplitude_lower = dsdyIWWamplitude(ttildemin, v);
    double amplitude = betaf*(1 - v)*(1 - v)*(1 - v)/v * (amplitude_upper - amplitude_lower);
    
    return prefactor*chiIWW(energy)*amplitude;
}



int main() {
    // creat file for output
    std::ofstream outfile("cs-output.csv");
    double energy = 100;

    std::cout << chiIWW(energy) << std::endl;
    std::cout << DifferentialCrossSection(energy, 0.5);
    // create vector of y values
    double ymin = mf/energy;
    double ymax = 1 - mass/energy;
    double size = 300;
    double step = (ymax - ymin)/(size-1);
    std::vector<double> y_vals;
    std::vector<double> cs_vals;
    for(int i = 0; i < size; i++) {
        y_vals.push_back(ymin+i*step);
        cs_vals.push_back(DifferentialCrossSection(energy, y_vals.at(i)));
        // outfile << y_vals.at(i) << "," << cs_vals.at(i) << "\n";
        // std::cout << y_vals.at(i) << "," << cs_vals.at(i) << "\n";
    }
    outfile.close();

    return 1;
}