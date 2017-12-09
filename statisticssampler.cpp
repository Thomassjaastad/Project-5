#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{


}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    if(!m_file.good()) {
        m_file.open("statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }
    }

    // Print out values here
    m_file <<  m_kineticEnergy << " " << m_potentialEnergy << " " << totalEnergy() << " " << m_temperature << " " << system.time() << " " << m_diffusionConstant << "\n";

}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    sampleDiffusionConstant(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();

    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential().potentialEnergy();
}


void StatisticsSampler::sampleTemperature(System &system)
{
    // Hint: reuse the kinetic energy that we already calculated

    m_temperature = 2.0/3.0*(m_kineticEnergy/(system.atoms().size()));
}

void StatisticsSampler::sampleDensity(System &system)
{
    vec3 systemsize = system.systemSize();
    double volume =  systemsize(0)*systemsize(1)*systemsize(2);
    //std::cout << "Particle density is :" << "" << system.atoms().size()/volume << std::endl;

}

void StatisticsSampler::sampleDiffusionConstant(System &system)
{
    m_diffusionConstant = 0;
    double r = 0;
    for(Atom *atom : system.atoms()) {
       r += (atom->position - atom->initialPosition).lengthSquared();
       //std::cout << "Test" << std::endl;
    }


    m_diffusionConstant = r/(6.0*system.time()*system.atoms().size());

}
