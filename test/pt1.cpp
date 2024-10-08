#include <iostream>
#include <fstream>

#include "gaus_overlap.h"

int main(void) {

//test to make sure input file is read correctly
    //std::string file_path = "../sample_input/numerical/3.txt";
    // Gaussian ob(file_path);
    // std::cout << ob.get_x0() << " " << ob.get_alpha();
    // std::cout << " " << ob.get_l() << std::endl;

    // bool second_atom = true;
    // Gaussian ob2(file_path, second_atom);
    // std::cout << ob2.get_x0() << " " << ob2.get_alpha();
    // std::cout << " " << ob2.get_l() << std::endl;
    
//test that the evaluation function works correctly
    // Gaussian ob3(file_path);
    // std::cout << "Single orbital evaluation: " << ob3.eval(.25) << std::endl;

//test that the riemann summ approximator works for a single orbital
    // Gaussian ob4(file_path);
    // float start = -1.5;
    // float end = 3.5;
    // int steps = 50;
    // std::cout << ob4.riemann_approx(start, end, steps) << std::endl;

//test that the riemann summ approximator works for overlapping orbitals
    // bool second_atom = true;
    
    // Gaussian ob4(file_path);
    // Gaussian ob5(file_path, second_atom);
    
    // float start = -1.5;
    // float end = 3.5;
    // int steps = 500;

    // double approx = ob4.riemann_approx(start, end, steps, ob5);

    // std::cout << "1d numerical overlap integral between Gaussian functions is ";
    // std::cout << approx << std::endl;
    
//comparing given answers with calculated answers
    //yes I know this is a yucky implementation, this is just what I
    // could muster up at the moment
    bool second_atom = true;
    float start = -1.5;
    float end = 3.5;
    int steps = 500;

    std::string case1_in = "../sample_input/numerical/1.txt";
    std::string case1_out = "../sample_output/numerical/1.txt";
    std::string case2_in = "../sample_input/numerical/2.txt";
    std::string case2_out = "../sample_output/numerical/2.txt";
    std::string case3_in = "../sample_input/numerical/3.txt";
    std::string case3_out = "../sample_output/numerical/3.txt";

    Gaussian ob6(case1_in);
    Gaussian ob7(case1_in, second_atom);

    std::ifstream file(case1_out);
    std::string str;
    std::string file_contents;
    while (std::getline(file, str))
    {
        file_contents += str;
        file_contents.push_back('\n');
    }  
    std::cout << "Expected output: " << file_contents;
    std::cout << "Actual calculated output: " << ob6.riemann_approx(start, end, steps, ob7);
    std::cout << std::endl << std::endl;

    Gaussian ob8(case2_in);
    Gaussian ob9(case2_in, second_atom);

    std::ifstream file2(case2_out);
    std::string file_contents2;
    while (std::getline(file2, str))
    {
        file_contents2 += str;
        file_contents2.push_back('\n');
    }  
    std::cout << "Expected output: " << file_contents2;
    std::cout << "Actual calculated output: " << ob8.riemann_approx(start, end, steps, ob9);
    std::cout << std::endl << std::endl;

    Gaussian ob10(case3_in);
    Gaussian ob11(case3_in, second_atom);

    std::ifstream file3(case3_out);
    std::string file_contents3;
    while (std::getline(file3, str))
    {
        file_contents3 += str;
        file_contents3.push_back('\n');
    }  
    std::cout << "Expected output: " << file_contents3;
    std::cout << "Actual calculated output: " << ob10.riemann_approx(start, end, steps, ob11);
    std::cout << std::endl << std::endl;
    

    return 0;
}