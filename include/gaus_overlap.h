// Homework 2:  Analytical and Numerical Integration for Overlap Integrals
// Laura Jones
// Chem 279 - Dr. Mayank Agrawal
// Last revisited: 10/1/2024 
// gaus_overlap.h contains all the necessary information about using gaussians to 
//      approximate the overlap between two s-type atomic orbitals

# pragma once

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

//abstract class for gaussian to derive from
class Integrand
{
    public:
        virtual double get_x0() const = 0; //coord for center of atom
        virtual double get_alpha() const = 0;
        virtual double get_l() const = 0;
};

class Gaussian: public Integrand
{
    private:
        double _x0, _alpha;
        int _l;

    public:
    //getter functions
        double get_x0() const {return _x0;};
        double get_alpha() const {return _alpha;};
        double get_l() const {return _l;};
    
    // evaluation functions
        double eval(double, Gaussian);
        double riemann_approx(float, float, int, Gaussian);

    //constructor
        Gaussian(std::string, bool next_atom = false);
        ~Gaussian();

    //misc
        //set parameters after reading them from input file
        void read_input(std::string, bool next_atom = false);

};

class normGaussian: public Integrand 
{
    private:
        std::vector <float> _coords;
        double _x0, _alpha;
        int _l;

    public:
    //getter functions
        double get_x0() const {return _x0;};
        double get_alpha() const {return _alpha;};
        double get_l() const {return _l;};
        std::vector<float> get_coords() const {return _coords;};

    //evaluating functions
        double factorial(double);
        double double_factorial(double);
        double combinations(double, double);
        double eval(normGaussian, int);
        double calc_center(normGaussian);
        std::vector<double> calc_3d(normGaussian);

    //constructor/destructor
        normGaussian(std::string, bool next_atom = false);
        ~normGaussian();

    // misc
        void read_input(std::string, bool next_atom = false);
};
