//
//  PenningTrap.hpp
//
//
//  Created by Elena Muñoz Rivas, Alejandro Carballido Mantecón, David Martínez Hernández and Antonio Gómez Garrido.
//
//
// Includes the class PenningTrap with its constructor and the variables needed to run the methods (magnetic field strength, applied potential, the dimension and the the Particle objects in the Penning trap). It contains some member functions as the forces, the electric and magnetic fields, a function to count how many particles are trapped and the Runge-Kutta 4th orther and the forward Eurler methods are defined here.

#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__
#include "Particle.hpp"

class PenningTrap
{
public:
    double B0_,V0_,d_,k_e,f_,wv_,interactions_,extension_;
    std::vector<Particle> Particle_;
    
    
    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, double f, double wv, int interactions, int extension);

    // Add a particle to the trap
    void add_particle(Particle p_in);
    
    void add_random_particle_Ca();
    
    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r, double t);
    
    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i, double t);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i, double t);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order without and with interactions:
    void evolve_RK4(double dt, double t);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt, double t);
    
    //A small function (either as part of the PenningTrap class or outside it) that can count how many of the particles are still inside the trap region
    int counter_particles();
    //The analytical solution which depends only on time and the initial state:
    arma::vec f_analytical(double t, Particle p_in);
       
    //The relative error of RK4 and Forward Euler for each step
    double relative_error(std::vector<int> n, int i);
};


#endif /* PenningTrap_hpp */

