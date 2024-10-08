#include <iostream>
#include <sstream>
#include <cmath> //so i can use pi

#include "gaus_overlap.h"

//the following is the source for the 1-d orbital overlap approximation class
Gaussian::Gaussian(std::string file, bool next_atom) {
    read_input(file, next_atom);
}
Gaussian::~Gaussian() {};

void Gaussian::read_input(std::string in_file, bool next_atom) {
    std::fstream data_file;
    std::string line;
    std::string s;

    data_file.open(in_file);
    getline(data_file,line, '\n');
    //if this is the 2nd atom, then skip to the next line
    if (next_atom == true) {
        getline(data_file, line, '\n');
    }

    std::stringstream ss(line);
    getline(ss, s, ' ');

    _x0 = std::stod(s);
    getline(ss, s, ' ');
    _alpha = std::stod(s);
    getline(ss, s, ' ');
    _l = std::stoi(s);
}

double Gaussian::eval(double x, Gaussian obj_2) {
    float obj_x0 = obj_2.get_x0();
    float obj_alpha = obj_2.get_alpha();
    int obj_l = obj_2.get_l();

    double part1 = pow((x - _x0), _l) * pow((x - obj_x0), obj_l);
    double part2 = (-_alpha*pow((x - _x0),2)) - (obj_alpha*pow((x - obj_x0),2));
    return part1 * exp(part2);
}

double Gaussian::riemann_approx(float start, float end, int steps, Gaussian obj_2) {
    double area = 0.0;
    double step_size = (end - start) / steps;

    for (double current = start; current < end; current+= step_size) {
        double mid = current + step_size / 2;
        double height = eval(mid, obj_2);

        area += (step_size) * height;
    }

    return area;
}

//below is the source for the 3-d overlap integral

normGaussian::normGaussian(std::string file, bool next_atom) {
    read_input(file, next_atom);
}

normGaussian::~normGaussian(){}

void normGaussian::read_input(std::string file, bool next_atom) {
        std::fstream data_file;
    std::string line;
    std::string s;

    data_file.open(file);
    getline(data_file,line, '\n');
    //if this is the 2nd atom, then skip to the next line
    if (next_atom == true) {
        getline(data_file, line, '\n');
    }

    std::stringstream ss(line);
    getline(ss, s, ' ');

    _coords.push_back(std::stof(s));
    getline(ss, s, ' ');
    _coords.push_back(std::stof(s));
    getline(ss, s, ' ');
    _coords.push_back(std::stof(s));
    getline(ss, s, ' ');
    
    _x0 = std::stod(s);
    getline(ss, s, ' ');
    _alpha = std::stod(s);
    getline(ss, s, ' ');
    _l = std::stoi(s);
}

double normGaussian::factorial(double n) {
    double result = n;
    if(n==0 || n==1) return 1;
    for(int i = n-1; i > 0; i--) {
        result *= i;
    }
    return result;
}

double normGaussian::double_factorial(double n) {
    double result = n;
    while((n-2) > 0 && n > 1) {
        n-=2;
        result*= n;
    }
    return result;
}

double normGaussian::combinations(double n, double m) {
    if (m > n) {
        return 0.0;
    }
    else {
        return factorial(n) / (factorial(m) * factorial(n - m));
    }
}

double normGaussian::calc_center(normGaussian second_atom) {
    double beta = second_atom.get_alpha();
    double x0b = second_atom.get_x0();
    return ((_alpha * _x0) + (beta*x0b)) / _alpha + beta; 
}

double normGaussian::eval(normGaussian second_atom, int dim) {
    ////Xa, XB are 2 x, y, z coords of the 2 subshells
    ///Xp = calc_center
    //l = angular momentum quantum number
    //constants
    double pi = M_PI;

    //values for second atom
    double beta = second_atom.get_alpha();
    //double x0b = second_atom.get_x0();
    double lb = second_atom.get_l();
    std::vector<float> coords_b = second_atom.get_coords();
    double dim_coordb = coords_b[dim];

    double xp = calc_center(second_atom);
    
    double result = 0.0;
//begin calculation in parts
    double pt1 = exp( -(_alpha*beta*pow((_coords[dim] - dim_coordb),2)) / (_alpha + beta));
    double pt2 = sqrt(pi/(_alpha+beta));

    double pt3 = 0.0;
    for(int i = 0; i <= _l; i++) {
        for(int j = 0; j <= lb; j++) {
            if((i+j) % 2 == 0 && (i+j) != 0) {
                double c1 = combinations(_l,i);
                double c2 = combinations(lb,j);

                double fact = double_factorial(i + j - 1);
                pt3 += c1*c2*fact*(pow((xp - _coords[dim]),_l-i)*pow((xp - dim_coordb),lb-j)) \
                          / pow((2*(_alpha+beta)),(i+j/2));
            }
            }
        }
    return pt1*pt2*pt3;
}

std::vector<double> normGaussian::calc_3d(normGaussian second_atom) {
    std::vector<double> vals_3d;
    for(int i = 0; i < 3; i++) {
        vals_3d.push_back(eval(second_atom, i));
    }
    return vals_3d;
}
