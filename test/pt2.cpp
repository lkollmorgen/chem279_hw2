#include <iostream>
#include <fstream>

#include "gaus_overlap.h"

int main(void) {

//test to make sure input file is read correctly
    std::string file_path = "../sample_input/analytical/3.txt";
    // normGaussian ob(file_path);
    // std::vector<float> ob_coords = ob.get_coords();
    // std::cout << ob_coords[0] << " " << ob_coords[1] << " " << ob_coords[2];
    // std::cout << " " << ob.get_x0() << " " << ob.get_alpha();
    // std::cout << " " << ob.get_l() << std::endl;

    // bool second_atom = true;
    // normGaussian ob2(file_path, second_atom);
    // threeD ob2_coords = ob2.get_coords();
    // std::cout << ob2_coords.x << " " << ob2_coords.y << " " << ob2_coords.z;
    // std::cout << " " << ob2.get_x0() << " " << ob2.get_alpha();
    // std::cout << " " << ob2.get_l() << std::endl;
    
//test that the evaluation function works correctly
    // normGaussian ob3(file_path);
    // std::cout << "factorial test: " << ob3.factorial(5) << std::endl;
    // std::cout << "combinations test: " << ob3.combinations(3, 5) << std::endl;
    // std::cout << "double factorial test: " << ob3.double_factorial(5) << std::endl;

//test that the evaluation function works for overlapping orbitals in all 3 dimensions
    bool second_atom = true;
    
    normGaussian ob4(file_path);
    normGaussian ob5(file_path, second_atom);

    std::vector<double> eval_3d = ob4.calc_3d(ob5);

    std::cout << "3d numerical overlap integral between Gaussian functions is ";
    for(int i = 0; i < 3; i++) {
        std::cout << eval_3d[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}