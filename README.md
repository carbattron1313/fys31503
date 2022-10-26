Study of the orbits of particles inside a Penning trap in C++
---------------------------

There are one .cpp file, Particle.cpp and a pair of .hpp files, Particle.hpp and PenningTrap.hpp.

Particle.hpp
------------
Includes the class Particle with its constructor and the properties of a single Particle.

PenningTrap.hpp
---------------
Includes the class PenningTrap with its constructor and the variables needed to run the methods (magnetic field strength, applied potential, the dimension and the the Particle objects in the Penning trap). It contains some member functions as the forces, the electric and magnetic fields, a function to count how many particles are trapped and the Runge-Kutta 4th orther and the forward Eurler methods are difined here.

Particle.cpp
------------
Main program where we use the classes. It simulates the movement of a single particle and two particles in the Penning trap for a total time using the Euler and RK4 method, obtaining the relative errors and the error convergence rate, the fraction of particles that are still trapped after some time when the potential is time-dependant... among other things.


Build: g++ Particle.cpp -o Paricle.exe 
Run: ./Particle.exe

PenningTrap.py
--------------
Python script that reads the datas obtained in the C++ and generates plots of the trajectories of one Particle and then with two with and without particle interactions, the 3D trayectories, the velocity of different directions against its direcctions, the relative errors made with RK4 and Euler method, the fraction of particles that are trapped in the Penning Trap for some frecuencies and the fraction of trapped particles versus frequency for one case with Coulomb interaction and other without. Plots are saved as pdf files.

Run command: python3 PenningTrap.py

