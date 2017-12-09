#include "lennardjones.h"
#include "system.h"
#include "math.h"

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop

    //Constants in force
    double sigma6 = sigma()*sigma()*sigma()*sigma()*sigma()*sigma();
    double sigma12 = 2*sigma6*sigma6;

    for (int i = 0; i < system.atoms().size(); i++){
        for(int j = i+1; j < system.atoms().size(); j++){
            //Creating vector r_ij
            Atom*atom_i = system.atoms()[i];
            Atom*atom_j = system.atoms()[j];
            vec3 r_ij = atom_i->position - atom_j->position;

            //Periodic boundary conditions
            //3 - dimensions

            for (int k = 0; k < 3; k++){
                if (r_ij[k] > system.systemSize()[k]*0.5){
                    r_ij[k] = r_ij[k] - system.systemSize()[k]; //new conditions? or in system.cpp

                }
                else if (r_ij[k] <= - system.systemSize()[k]*0.5){
                    r_ij[k] = r_ij[k] + system.systemSize()[k];

                }
            }

            double len_r_ij = r_ij.lengthSquared();
            //vec3 unit_r(r_ij[0]/len_r_ij, r_ij[1]/len_r_ij , r_ij[2]/len_r_ij);

            double len_r_ij6 = len_r_ij*len_r_ij*len_r_ij;  //6 power term
            double len_r_ij12 = len_r_ij6*len_r_ij6;                                   //12 power term

            m_potentialEnergy += 4*epsilon()*(sigma12/(2*len_r_ij12) - sigma6/len_r_ij6);
            //Calculting forces on atom i
            vec3 dF = 24*(epsilon()/r_ij.lengthSquared())*(sigma12/len_r_ij12 - sigma6/len_r_ij6)*r_ij;
            atom_i->force += dF;
            atom_j->force -= dF;
        }
    }

}
