//
//  Particle.hpp
//  
//
//  Created by Elena Muñoz Rivas, Alejandro Carballido Mantecón, David Martínez Hernández and Antonio Gómez Garrido.
//
//
// Includes the class Particle with its constructor and the properties of a single Particle.

#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <stdio.h>
#include <armadillo>
#include <math.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>

class Particle
{
public:
    double q_,m_;
    arma::vec r_,v_;
    
    Particle(double q, double m, arma::vec r, arma::vec v);
    
};


#endif /* Particle_hpp */
