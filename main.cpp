#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    /*
    for(int temp= 0; temp < 1200; temp = temp + 20){

    }
    */

    double numberOfUnitCells = 5;
    double initialTemperature = UnitConverter::temperatureFromSI(200); // measured in Kelvin
    double latticeConstant =    UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms


    /*
    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));
    */
    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.
    /*
    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    */

    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    system.potential().setEpsilon(UnitConverter::temperatureFromSI(119.8));
    system.potential().setSigma(UnitConverter::lengthFromAngstroms(3.405));
    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    //cout << statisticsSampler.diffusionConstant() << endl;
    //exit(1);
    IO movie("movie.xyz"); // To write the state to file
    /*
    cout << setw(17) << "Timestep" <<
            setw(17) << "Time" <<
            setw(17) << "Temperature" <<
            setw(17) << "KineticEnergy" <<
            setw(17) << "PotentialEnergy" <<
            setw(17) << "TotalEnergy" <<
            setw(17) << "DiffConstant" << endl;
    */
    double averageTemp = 0;
    double averageDiff = 0;
    for(int timestep=0; timestep<5000; timestep++) {
        system.step(dt);
        statisticsSampler.sample(system);

        if (timestep > 2999){
            averageTemp += UnitConverter::temperatureToSI(statisticsSampler.temperature());
            averageDiff += UnitConverter::diffusionToSI(statisticsSampler.diffusionConstant());
        }
        if (timestep == 4999){
            cout << "Diffusion_constant_is: " << statisticsSampler.diffusionConstant() << " " <<"Mean_Diffusion_constant:" << " " << averageDiff/2000.0 << " " << "Average_temperature_is: " << averageTemp/2000.0 << endl;
        }
        /*
        if( timestep % 100 == 0 ) {
            // Print the timestep every 100 timesteps
            cout << setw(17) << system.steps() <<
                    setw(17) << system.time() <<
                    setw(17) << UnitConverter::temperatureToSI(statisticsSampler.temperature()) <<
                    setw(17) << statisticsSampler.kineticEnergy() <<
                    setw(17) << statisticsSampler.potentialEnergy() <<
                    setw(17) << statisticsSampler.totalEnergy() <<
                    setw(17) << UnitConverter::diffusionToSI(statisticsSampler.diffusionConstant()) << endl;
        }
        */
        if( timestep % 10 == 0 ) movie.saveState(system);
    }

    movie.close();

    return 0;
}
