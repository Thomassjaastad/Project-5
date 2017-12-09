#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include <ctime>

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {

    int nrOfAtoms = m_atoms.size();

    for (int i = 0; i < nrOfAtoms; i++){
        Atom*a = m_atoms[i];


        if (a->position(0) < 0){
            a->position(0) += systemSize().x();
            a->initialPosition(0) += systemSize().x();

        }
        else if (a->position(0) >= systemSize().x()){
            a->position(0) -= systemSize().x();
            a->initialPosition(0) -= systemSize().x();
        }

        if (a->position(1) < 0){
            a->position(1) += systemSize().y();
            a->initialPosition(1) += systemSize().y();
        }
        else if (a->position(1) >= systemSize().y()){
            a->position(1) -= systemSize().y();
            a->initialPosition(1) -= systemSize().y();
        }

        if (a->position(2) < 0){
            a->position(2) += systemSize().z();
            a->initialPosition(2) += systemSize().z();
        }
        else if (a->position(2) >= systemSize().z()){
            a->position(2) -= systemSize().z();
            a->initialPosition(2) -= systemSize().z();
        }

    }

    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeTotalMomentum() {
    //p = mv, momentum of one atom
    //Sum up momentum and set to zero, what needs to be zero?
    //Need to define variables in header file.
    int nrOfAtoms = m_atoms.size();

    totalMomentum.zeros();
    totalVelocity.zeros();
    averageVelocity.zeros();

    for (int i = 0; i < nrOfAtoms ; i++){
        Atom*atom = m_atoms[i];
        totalVelocity += atom->velocity;

    }
    averageVelocity = totalVelocity/nrOfAtoms;


    for (int i = 0; i < nrOfAtoms; i++){
        Atom*atom = m_atoms[i];
        atom->velocity -= averageVelocity;
    }


    //std::cout << "Total momentum is:" << " " << totalMomentum << std::endl;

    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    //Position in one unit cell
    double b = latticeConstant;

    int N = numberOfUnitCellsEachDimension;
    m_systemSize[0] = N*b;
    m_systemSize[1] = N*b;
    m_systemSize[2] = N*b;
    //Going over each dimension and creating crystal structure

    for(int i=0; i< N; i++) {
        for(int j=0; j < N; j++){
            for(int k=0; k< N; k++){
                //Creating four atoms for one unit cell
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                //Setting position for one unit cell and adding its relative origin
                vec3 R = {i*b, j*b, k*b};
                atom1->position.set(0,0,0);
                atom1->position += R;
                atom2->position.set(b/2, b/2, 0);
                atom2->position += R;
                atom3->position.set(0, b/2, b/2);
                atom3->position += R;
                atom4->position.set(b/2,0,b/2);
                atom4->position += R;

                //Setting initial position for each atoms in lattice
                atom1->initialPosition = atom1->position;
                atom2->initialPosition = atom2->position;
                atom3->initialPosition = atom3->position;
                atom4->initialPosition = atom4->position;


                //Random::myrandom(500);
                atom1->resetVelocityMaxwellian(temperature);
                atom2->resetVelocityMaxwellian(temperature);
                atom3->resetVelocityMaxwellian(temperature);
                atom4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);
                m_atoms.push_back(atom2);
                m_atoms.push_back(atom3);
                m_atoms.push_back(atom4);

            }
        }
    }

}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
        atom->position;


    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
