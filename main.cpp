#include <iostream>
#include <posit/posit>
#include "jack_settings.hpp"

//WRONG! Sum after += is less than -1.0: -2.68435e+08, after taking sum (2.6226e-05) += quire_mul(-0.016571, 0.000999451)
//We essentially took 2.6226e-05 += -1.66893e-05 and somehow ended up with -2.68435e+08
//Incorrect sum as a quire: -1: 111111111111111111111111111111_111111111111111111111111111111111111111111111111111111111.11111111111111110101111111001010000000000000000000000000
//Delta Quire: -1: 000000000000000000000000000000_000000000000000000000000000000000000000000000000000000000.00000000000000010001010111011101000000000000000000000000, (as posit) = -1.66893e-05
//Row: 266, Acoefs[i]: -0.016571, Acols[i]: 398, xcoefs[398]: 0.000999451


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