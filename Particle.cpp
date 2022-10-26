//
//  Particle.cpp
//
//  Created by Elena Muñoz Rivas, Alejandro Carballido Mantecón, David Martínez Hernández and Antonio Gómez Garrido.
//
//  This program simulates the movement of a single particle and two particles in the Penning trap for a total time using the Euler and RK4 method, obtaining the relative errors and the error convergence rate, the fraction of particles that are still trapped after some time when the potential is time-dependant... among other things.
//

#include "Particle.hpp"
#include "PenningTrap.hpp"
//Definition of the constructor Particle
Particle::Particle(double q, double m, arma::vec r, arma::vec v)
{
    q_=q;
    m_=m;           // Assignment of the variables that caracterize each particle;
    r_=r;           //them being the charge, mass, possition and velocity
    v_=v;
}

//Definition of the constructor PenningTrap
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double f, double wv, int interactions, int extension)
{
    B0_=B0_in;
    V0_=V0_in;                               //Assignment of the variables that define our trap
    d_=d_in;
    k_e=1.38935333e+5;                       // Definition of the Coulomb's constant in the magnitudes chosen for our problem
    std::vector<Particle> Particle_;         //Declare a vector that stores every object of the type Particle that enters our Penning trap
    f_=f;                                    //Amplitude of the oscilations in the electric potenctial (needed for Problem 9)
    wv_=wv;                                  //Angular frequency of the same oscilations (needed for Problem 9)
    interactions_=interactions;              //Variable to turn off or on the particle interactions. In order to take into account the interactions, we must set this variable to be =1
    extension_=extension;              //Variable to change the potential and make it time dependent or not. To have a time dependant potential it must be =1.
}

//Definition of add_particle. This function introduces an object of the type Particle inside our Penning trap.
void PenningTrap::add_particle(Particle p_in)
{
    Particle_.push_back(p_in);
}

//Definition of a function that introduces a particle of calcuium ion with random position and velocity inside the Penning trap.
void PenningTrap::add_random_particle_Ca()
{
    arma::arma_rng::set_seed_random();
    arma::vec r_rand = arma::vec(3).randn() * 0.1 * d_;  // random initial position
    arma::vec v_rand = arma::vec(3).randn() * 0.1 * d_;  // random initial velocity
    Particle rand_part(1.,40.078,r_rand,v_rand);         // Create a Particle object with the radomly generated position and velocity
    Particle_.push_back(rand_part);                      // Introduce the particle above inside the trap
}


//This function returns the value for the electric field at a certain time in a given point in space.
arma::vec PenningTrap::external_E_field(arma::vec r,double t)
{
    double potential;
    
    if (extension_==1){
        potential=V0_*(1.+f_*std::cos(wv_*t));           // The variable extension is introduced to enable us to choose between a time-dependant potential and a time-independent one.
    }
    else{
        potential=V0_;
    }
    
    arma::vec E;
    
    if (norm(r, 2)>d_){                                 //"turns off" the electric field outside of the trap
        E = {0.,0.,0.};
    }
    else{
        double ratio = potential/(d_*d_);               //Definition of the electric field
        E={r[0]*ratio, r[1]*ratio, -2.*ratio * r[2]};
    }
    return E;
}


//Definition of external_B_field. This gives out the value of the magnetic field at a certain point in space.
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    arma::vec B = {0., 0., B0_};

    if (norm(r, 2)>d_){                                 //"turns off" the magnetic field outside of the trap
        B(2) = 0.;
    }
    return B;                         //Definition of the magnetic field
}


//Definition of the Coulomb force between two particles. The integer variables make refferences to the index that defines the particle
arma::vec PenningTrap::force_particle(int i, int j)
{
    if (i==j){                                          //A particle doesn't exert an electric force on itself
        arma::vec Null = arma::vec(3).fill(0.);
        return Null;
    }
    else{
        arma::vec force_particle = arma::vec(3);
        arma::vec substraction = Particle_[i].r_-Particle_[j].r_;       //Vector that goes from particle j to particle i. This is the direction of the force exerted by particle j on particle i
            force_particle = Particle_[i].q_*Particle_[j].q_*k_e*substraction/(pow(std::sqrt(substraction(0)*substraction(0)+substraction(1)*substraction(1)+substraction(2)*substraction(2)),3));    //Definition of the Coulomb force that particle j exerts on particle i
        
        return force_particle;
    }
}


// Define the Lorentz force acting on a given particle at a given time.
arma::vec PenningTrap::total_force_external(int i, double t)
{
    return Particle_[i].q_*(external_E_field(Particle_[i].r_,t) + arma::cross(Particle_[i].v_,external_B_field(Particle_[i].r_)));
}


// Define the total force acting on particle_i caused by the interaction with the rest of the particles in the trap
arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec total_force_particles = arma::vec(3).fill(0.);        //Create an emptuy vector that will store the force
    for (int j=0; j<Particle_.size(); j++){
        total_force_particles+=force_particle(i,j);                 //Add the contribution of every particle to this force
    }

    return total_force_particles;
}


// Define the total force on particle_i from both external fields and interactions with other particles
arma::vec PenningTrap::total_force(int i, double t)
{
    arma::vec total_force= total_force_external(i,t);
    if (interactions_==1){                                     //The variable "interactions" is equivalent to the variable "extension" mentioned above. In this case, it helps us choose whether we take into account the interactions between particles or not.
        total_force+=total_force_particles(i);
    }
    
    return total_force;
}

//Define the analytical function:
arma::vec PenningTrap::f_analytical(double t, Particle p_in){

double q_=p_in.q_;
double m_=p_in.m_;
arma::vec r_= p_in.r_;
arma::vec v_= p_in.v_;

//Define the f vector to include the values of the directions (x, y, z)
arma::vec f = arma::vec(3);
    
//Define de variables to find the analytical solution:
double A_n, A_s;
double w0,wz2,w_n,w_s;
double x, y, z;
w0= B0_*q_/m_;
wz2= 2*V0_*q_/(m_*d_*d_);
w_n=(w0+std::sqrt(w0*w0-2*wz2))/2;
w_s=(w0-std::sqrt(w0*w0-2*wz2))/2;

A_n= (v_(1)+w_s*r_(0))/(w_s-w_n);
A_s= -(v_(1)+w_n*r_(0))/(w_s-w_n);

//Values of the cartesians axis analytically:
x= A_n*std::cos(w_n*t)+A_s*std::cos(w_s*t);
y= -(A_n*std::sin(w_n*t)+A_s*std::sin(w_s*t));
z= r_(2)*std::cos(std::sqrt(wz2)*t);

f(0)= x; f(1)= y; f(2)= z;

return f;

}


//Define the relative error and print it to a .txt file
double PenningTrap::relative_error(std::vector<int> n, int i){

//Create a vector that contains the .txt filenames for the different timesteps, n_k:
std::vector<std::string> ofiles_txt;

for (int l=0;l<n.size();l++){
ofiles_txt.push_back("relative_error"+std::to_string(n[l])+".txt");
}

double t_re;
double dt_re;
for (int k=0; k<n.size(); k++){

Particle p_in= Particle_[i];
  std::ofstream file;
        file.open(ofiles_txt[k]);
        //Create 2 copies of the particle so they can be updated separately for each method
        Particle p_in_rk4= p_in;
        Particle p_in_eu= p_in;
        
               
        t_re= 0.;
        dt_re=50./n[k];
        
         //Print the initial values into a .txt
        file << std::setw(20) << std::setprecision(5) << 0
        << std::setw(20) << std::setprecision(5) << norm((p_in.r_ - f_analytical(0,p_in)),2)/norm(f_analytical(0., p_in),2)
        << std::setw(20) << std::setprecision(5) << norm((p_in.r_ - f_analytical(0,p_in))/(f_analytical(0.,p_in)),2)
        << std::endl;
      
        //We evolve the system by one time step and store the results in the particle copies

      for (int j=0; j<n[k];j++){
      
      t_re+=dt_re;
      Particle_[i]= p_in_rk4;
      evolve_RK4(dt_re, t_re);
      p_in_rk4=Particle_[i];
      Particle_[i]= p_in_eu;
      evolve_forward_Euler(dt_re, t_re);
      p_in_eu=Particle_[i];
      
      file << std::setw(20) << std::setprecision(5) << t_re
        << std::setw(20) << std::setprecision(5) << norm(p_in_rk4.r_ - f_analytical(t_re,p_in),2)/norm(f_analytical(t_re,p_in),2)
        << std::setw(20) << std::setprecision(5) << norm(p_in_eu.r_ - f_analytical(t_re,p_in),2)/norm(f_analytical(t_re,p_in),2)
        << std::endl;
       }
        //Return the particle to the initial state so that the following loop starts with the original particle
       Particle_[i]=p_in;
        
file.close();
}
return 0;
}

//Function that evolves the system by one time step applying the Runge-Kutta 4 algorithm.
void PenningTrap::evolve_RK4(double dt, double t)
{
    
    std::vector<Particle> Particle_copy= Particle_;     //Make a copy of the collection of particles inside the trap
    
    std::vector<arma::vec> kr1,kr2,kr3,kr4;             //These vectors will store the different kr's and kv's of every particle, which are needed for the RK4 algorithm.
    std::vector<arma::vec> kv1,kv2,kv3,kv4;
    
    arma::vec kr1i = arma::vec(3);
    arma::vec kr2i = arma::vec(3);                      //Declaration of the variables where we will compute the kr's and kv's of each individual particle.
    arma::vec kr3i = arma::vec(3);
    arma::vec kr4i = arma::vec(3);
    arma::vec kv1i = arma::vec(3);
    arma::vec kv2i = arma::vec(3);
    arma::vec kv3i = arma::vec(3);
    arma::vec kv4i = arma::vec(3);
    
    for (int i=0; i<Particle_.size();i++){             //Loop for the computation of kr1 and kv1 of every particle
        kv1i=dt*total_force(i,t)/Particle_[i].m_;
        kr1i=dt*Particle_[i].v_;
        
        kr1.push_back(kr1i);
        kv1.push_back(kv1i);
        
    }
    
    for (int i=0; i<Particle_.size();i++){        //Now use the values computed in the previous loop to update the velocities and positions of the particles. They will be needed in order to compute the next kr's and kv's
        
        
        Particle_[i].r_=Particle_copy[i].r_+0.5*kr1[i];
        Particle_[i].v_=Particle_copy[i].v_+0.5*kv1[i];
        
    }
    for (int i=0; i<Particle_.size();i++){      //Computation of kr2 and kv2
        kv2i=dt*total_force(i,t+0.5*dt)/Particle_[i].m_;
        kr2i=dt*Particle_[i].v_;
        
        kr2.push_back(kr2i);
        kv2.push_back(kv2i);
        
    }
        
    for (int i=0; i<Particle_.size();i++){    //Next updte of the positions and velocities needed for kr3 and kv3
        
        Particle_[i].r_=Particle_copy[i].r_+0.5*kr2[i];
        Particle_[i].v_=Particle_copy[i].v_+0.5*kv2[i];
    }
    for (int i=0; i<Particle_.size();i++){    //Computation of kr3 and kv3
        kv3i=dt*total_force(i,t+0.5*dt)/Particle_[i].m_;
        kr3i=dt*Particle_[i].v_;
        
        kr3.push_back(kr3i);
        kv3.push_back(kv3i);
    }
        
    for (int i=0; i<Particle_.size();i++){      //Next update
        
        Particle_[i].r_=Particle_copy[i].r_+kr3[i];
        Particle_[i].v_=Particle_copy[i].v_+kv3[i];
    }
    for (int i=0; i<Particle_.size();i++){      //Computation of kr4 and kv4
        kv4i=dt*total_force(i,t+dt)/Particle_[i].m_;
        kr4i=dt*Particle_[i].v_;
        
        kr4.push_back(kr4i);
        kv4.push_back(kv4i);
    }
            
    for (int i=0; i<Particle_.size();i++){      //Definitive update of the velocity and position of the particles. For this one, we need the original positions and velocities of the particles, which have remained unchanged inside the copy of the collection of particles (that is why we needed it in the first time).

        Particle_[i].r_=Particle_copy[i].r_+1/6.*(kr1[i]+2*kr2[i]+2*kr3[i]+kr4[i]);
        Particle_[i].v_=Particle_copy[i].v_+1/6.*(kv1[i]+2*kv2[i]+2*kv3[i]+kv4[i]);

    }
}


//Define the forward Euler method
void PenningTrap::evolve_forward_Euler(double dt, double t)
{
    for (int i=0; i<Particle_.size();i++){
        
        Particle Particle_copy = Particle_[i];
        Particle_copy.v_=Particle_[i].v_+dt*total_force(i,t)/Particle_[i].m_;
        Particle_copy.r_=Particle_[i].r_+dt*Particle_[i].v_;
        Particle_[i]=Particle_copy;
    }
    
}

//A small function (either as part of the PenningTrap class or outside it) that can count how many of the particles are still inside the trap region
int PenningTrap::counter_particles()
{
    int num_particles = 0;
    
    for (int i=0; i<Particle_.size(); i++){
        if (norm(Particle_[i].r_, 2)<=d_){      //Only take into account the particles which are inside thhe influence of the Penning trap
            num_particles += 1;
        }
    }
    
    return num_particles;
}



int main(int argc, char *argv[])
{
    int interactions;
    int extension;
    interactions=atoi(argv[1]);
    extension=atoi(argv[2]);
    
//============================================================================================================
    
    //Set a filename for the .txt that stores the information for the first Penning Trap
    std::string filename = "PenningTrap1.txt";
    
    std::ofstream ofile;
    ofile.open(filename);
    
    //Define the arguments for the constructor of the first Penning Trap
    double V0_in=2.41e+6;
    double B0_in=9.65e+1;
    double d_in=500.;
    double f = 0.;
    double wv = 0.;

    PenningTrap PenningTrap1(B0_in,V0_in,d_in,f,wv,interactions,extension);
    
    //Define the arguments for the particle and add it to the Penning Trap
    arma::vec r01 = arma::vec("20. 0. 20.");
    arma::vec v01 = arma::vec("0. 25. 0.");
    Particle Particle1_1(1.,40.078,r01,v01);
    
    PenningTrap1.add_particle(Particle1_1);
    
    // Define the settings for the simulation (50 microseconds, 1 particle)
    int n = 4000;
    double dt = 50./n;
    double z_1 = r01(2);
    double t = 0.;

    //Print the initial values into a .txt
    ofile << std::setw(20) << std::setprecision(8) << t
    << std::setw(20) << std::setprecision(8) << z_1 << std::endl;
    
    //Run the program evolving the system in time (n times) using Runge-Kutta.
    for (int i=0; i<n; i++){
        PenningTrap1.evolve_RK4(dt,t);
        //We could have also use the method Forward Euler instead:
            //PenningTrap1.evolve_forward_Euler(dt);
        
        //Update the actual time of the system and printing it with the position z of the particle
        t+=dt;
        
        ofile << std::setw(20) << std::setprecision(5) << t;
        
        for (int j=0; j<PenningTrap1.Particle_.size();j++){
            ofile << std::setw(20) << std::setprecision(5) << PenningTrap1.Particle_[j].r_(2);
        }
        ofile<<std::endl;
    }
    ofile.close();
    
//============================================================================================================

    
    //Set a filename for the .txt that stores the information for the second Penning Trap
        std::string filename2 = "PenningTrap2.txt";
        
        std::ofstream ofile2;
        ofile2.open(filename2);

        
    //Set the second Penning Trap with the same arguments as the other one.
        PenningTrap PenningTrap2 = PenningTrap(B0_in,V0_in,d_in,f,wv,interactions,extension);
        
    
    //Define the arguments for two particles and add them to the Penning Trap
    
        Particle Particle1_2 = Particle(1.,40.078,r01,v01);

        arma::vec r02 = arma::vec("25. 25. 0.");
        arma::vec v02 = arma::vec("0. 40. 5.");
        Particle Particle2_2 = Particle(1.,40.078,r02,v02);
        
        PenningTrap2.add_particle(Particle1_2);
        PenningTrap2.add_particle(Particle2_2);
    
    // Define the settings for the simulation (50 microseconds, 2 particles)
    double x_2,y_2,z_2,vx_2,vy_2,vz_2;
    double dt2 = 50./n;
    double t2 = 0.;
        
    //Print the initial values into a .txt
    ofile2 << std::setw(20) << std::setprecision(5) << 0
        << std::setw(20) << std::setprecision(5) << r01(0)
        << std::setw(20) << std::setprecision(5) << r01(1)
        << std::setw(20) << std::setprecision(5) << r01(2)
        << std::setw(20) << std::setprecision(5) << v01(0)
        << std::setw(20) << std::setprecision(5) << v01(1)
        << std::setw(20) << std::setprecision(5) << v01(2)
        << std::setw(20) << std::setprecision(5) << r02(0)
        << std::setw(20) << std::setprecision(5) << r02(1)
        << std::setw(20) << std::setprecision(5) << r02(2)
        << std::setw(20) << std::setprecision(5) << v02(0)
        << std::setw(20) << std::setprecision(5) << v02(1)
        << std::setw(20) << std::setprecision(5) << v02(2) << std::endl;
        
    //Same as before but with two particles and printing their x and y components for each step
        for (int i=0; i<n; i++){
            PenningTrap2.evolve_RK4(dt2,t2);
            t2+=dt2;
            ofile2 << std::setw(20) << std::setprecision(5) << t2;
            for (int j=0; j<PenningTrap2.Particle_.size();j++){
                ofile2 << std::setw(20) << std::setprecision(5) << PenningTrap2.Particle_[j].r_(0)
                    << std::setw(20) << std::setprecision(5) << PenningTrap2.Particle_[j].r_(1)
                    << std::setw(20) << std::setprecision(5) << PenningTrap2.Particle_[j].r_(2)
                    << std::setw(20) << std::setprecision(5) << PenningTrap2.Particle_[j].v_(0)
                    << std::setw(20) << std::setprecision(5) << PenningTrap2.Particle_[j].v_(1)
                    << std::setw(20) << std::setprecision(5) << PenningTrap2.Particle_[j].v_(2);
            }
            ofile2<<std::endl;
        
        }
        ofile2.close();

//============================================================================================================

    //For each time step, we compute the relative error fot both, Forward Euler and RK4, of a single particle (Particle 1) and print them in a .TXT file.

    PenningTrap PenningTrap_rel_error(B0_in,V0_in,d_in,f,wv,interactions,extension);
    Particle particle_rel_err(1.,40.078,r01,v01);
    PenningTrap_rel_error.add_particle(particle_rel_err);
    std::vector<int> n_steps = {4000,8000,16000,32000};
    PenningTrap_rel_error.relative_error(n_steps,0);
    
//============================================================================================================

    //Estimate the error convergence rate for RK4 and Forward Euler
            
            Particle particle_copy= particle_rel_err;
            double t_err,dt_err;
            //rk4
            
            
            arma::vec D_err_rk4= arma::vec(n_steps.size());
            for (int i=0; i<n_steps.size(); i++){
            
                t_err=0.;
                dt_err= 50./double(n_steps[i]);
                arma::vec A = arma::vec(n_steps[i]);
                for (int j=0; j<n_steps[i]; j++){
                    //Append the values of the difference between the analytical solution and the numerical approximation
                    A(j) = norm(PenningTrap_rel_error.f_analytical(t_err, particle_copy) - PenningTrap_rel_error.Particle_[0].r_ ,2);
                    
                    PenningTrap_rel_error.evolve_RK4(dt_err, t_err);
                    t_err+=dt_err;
                    
                }
                //For each n, store the maximum error
                PenningTrap_rel_error.Particle_[0]=particle_copy;
                D_err_rk4(i)=A.max();
            }
    
            double r_err_rk4=0.;
            //Compute the error convergence rate
            for (int k=1; k<n_steps.size(); k++){
                
                r_err_rk4 += (1./3.)*(std::log(D_err_rk4(k)/D_err_rk4(k-1)))/(std::log(double(n_steps[k-1])/double(n_steps[k])));
                
            }
            std::cout << "The error convergence rate for RK4 is " << r_err_rk4 <<std::endl;
            
            //forward-euler
            //Do the same for Forward Euler
            PenningTrap_rel_error.Particle_[0]= particle_copy;
            
            arma::vec D_err_eu= arma::vec(4);
            
            for (int i=0; i<n_steps.size(); i++){
            
                t_err=0.;
                dt_err= 50./double(n_steps[i]);
                arma::vec B = arma::vec(n_steps[i]);
                    for (int j=0; j<n_steps[i]; j++){
                        B(j)= norm(PenningTrap_rel_error.f_analytical(t_err,particle_copy) - PenningTrap_rel_error.Particle_[0].r_ ,2);
                        PenningTrap_rel_error.evolve_forward_Euler(dt_err, t_err);
                        t_err+=dt_err;
                    }
                PenningTrap_rel_error.Particle_[0]=particle_copy;
                D_err_eu(i)=B.max();
            }
    
            double r_err_eu=0.;
            for (int k=1; k<4; k++){
            
                r_err_eu += (1./3.)*(std::log(D_err_eu(k)/D_err_eu(k-1)))/(std::log(double(n_steps[k-1])/double(n_steps[k])));
                
            }
            std::cout << "The error convergence rate for Forward Euler is " << r_err_eu<<std::endl;
//============================================================================================================

    //Create random particles:
    PenningTrap PenningTrap100(B0_in,V0_in,d_in,f,wv,0,1);
    
    for (int i=0; i<100; i++){
        PenningTrap100.add_random_particle_Ca();
    }
    std::vector<Particle> Particle_random_copy = PenningTrap100.Particle_;
    
    double t100;
    double dt100 = 500./n;
    

    //============================================================================================================
    
    
    //Set a filename for the .txt that stores the information for the last Penning Trap
    std::string filename100f1 = "PenningTrap100f1.txt";
    
    std::ofstream ofile100f1;
    ofile100f1.open(filename100f1);
    
    //Value of the amplitude, f=0.1:
    PenningTrap100.f_=0.1;
    
    //Print the counter particles and the angular frecuency:
    for (double wv=0.2; wv<=2.5;wv+=0.02){
        PenningTrap100.wv_=wv;
        t100 = 0.;
        ofile100f1<<std::setw(20) << std::setprecision(5) << wv;
        for (int i=0; i<n; i++){
            PenningTrap100.evolve_RK4(dt100,t100);
            t100+=dt100;
        }
        ofile100f1<<std::setw(20) << std::setprecision(5) <<PenningTrap100.counter_particles()<<std::endl;
        PenningTrap100.Particle_=Particle_random_copy;
    }
    ofile100f1.close();
    
    //============================================================================================================
    
    //Set a filename for the .txt that stores the information for the last Penning Trap
    std::string filename100f4 = "PenningTrap100f4.txt";
    
    std::ofstream ofile100f4;
    ofile100f4.open(filename100f4);
    
    //Value of the amplitude, f=0.4:
    PenningTrap100.f_=0.4;
    
    //Print the counter particles and the angular frecuency:
    for (double wv=0.2; wv<=2.5;wv+=0.02){
        PenningTrap100.wv_=wv;
        t100 = 0.;
        ofile100f4<<std::setw(20) << std::setprecision(5) << wv;
        for (int i=0; i<n; i++){
            PenningTrap100.evolve_RK4(dt100,t100);
            t100+=dt100;
        }
        ofile100f4<<std::setw(20) << std::setprecision(5) <<PenningTrap100.counter_particles()<<std::endl;
        PenningTrap100.Particle_=Particle_random_copy;
    }
    ofile100f4.close();
    
    //============================================================================================================
    
    
    //Set a filename for the .txt that stores the information for the last Penning Trap
    std::string filename100f7 = "PenningTrap100f7.txt";
    
    std::ofstream ofile100f7;
    ofile100f7.open(filename100f7);
    
    //Value of the amplitude, f=0.7:
    PenningTrap100.f_=0.7;
    
    //Print the counter particles and the angular frecuency:
    for (double wv=0.2; wv<=2.5;wv+=0.02){
        PenningTrap100.wv_=wv;
        t100 = 0.;
        ofile100f7<<std::setw(20) << std::setprecision(5) << wv;
        for (int i=0; i<n; i++){
            PenningTrap100.evolve_RK4(dt100,t100);
            t100+=dt100;
        }
        ofile100f7<<std::setw(20) << std::setprecision(5) <<PenningTrap100.counter_particles()<<std::endl;
        PenningTrap100.Particle_=Particle_random_copy;
        
    }
    ofile100f7.close();

    
//============================================================================================================

    //Set a filename for the .txt that stores the information for the last Penning Trap
        std::string filenamezoom = "PenningTrap100zoom2.txt";
        
        std::ofstream ofilezoom;
        ofilezoom.open(filenamezoom);
    
    //Value of the amplitude, f=0.7:
    PenningTrap100.f_=0.7;
    
    //Print the counter particles and the angular frecuency when you do "zoom in" og the f=0.7 particle:
    //With and without interactions (with 0 or 1 in the terminal and saving in different .txt)
    for (double wv=0.42; wv<=0.47;wv+=0.0005){
        PenningTrap100.wv_=wv;
        t100 = 0.;
        ofilezoom<<std::setw(20) << std::setprecision(5) << wv;
        for (int i=0; i<n; i++){
            PenningTrap100.evolve_RK4(dt100,t100);
            t100+=dt100;
        }
        ofilezoom<<std::setw(20) << std::setprecision(5) <<PenningTrap100.counter_particles()<<std::endl;
        PenningTrap100.Particle_=Particle_random_copy;
    }
    
    ofilezoom.close();
    
    
    return 0;
    
}
