#include <iostream>
#include <posit/posit>
#include "jack_settings.hpp"

//
//Sum is -8.10623e-06, which was calculated by taking (-0.00828552 * 0.000999451)
//Sum is -2.68435e+08 after taking sum += quire_mul(-0.00828552, 0.000999451) (quire_mul(-0.00828552, 0.000999451) = -8.10623e-06
//Sum after += is less than -1.0 (as a quire): -1: 111111111111111111111111111111_111111111111111111111111111111111111111111111111111111111.11111111111111111110101010111000100000000000000000000000
//-262144 = (result + (xcoefs[266] * sum)) == (24.0156 + (0.00101471 * -2.68435e+08)


int main() {

    quireX sum = 2.6226e-05;
    positX sumAsPosit;
    positX argA = -0.016571;
    positX argB = 0.000999451;
    sumAsPosit.convert(sum.to_value());
    std::cout << "sum = " << sum << std::endl;
    std::cout << "sumAsPosit = " << sumAsPosit << std::endl;
    sum += sw::unum::quire_mul(argA, argB);
    std::cout << "sum after += == " << sum << std::endl;
    sumAsPosit.convert(sum.to_value());
    std::cout << "sum after +=(as posit) == " << sumAsPosit << std::endl;
    return 0;
}