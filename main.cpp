#include <iostream>
#include <posit/posit>
#include "jack_settings.hpp"

//
//Sum is -8.10623e-06, which was calculated by taking (-0.00828552 * 0.000999451)
//Sum is -2.68435e+08 after taking sum += quire_mul(-0.00828552, 0.000999451) (quire_mul(-0.00828552, 0.000999451) = -8.10623e-06
//Sum after += is less than -1.0 (as a quire): -1: 111111111111111111111111111111_111111111111111111111111111111111111111111111111111111111.11111111111111111110101010111000100000000000000000000000
//-262144 = (result + (xcoefs[266] * sum)) == (24.0156 + (0.00101471 * -2.68435e+08)


int main() {
    quireX testSumQuire1 = -8.10623e-06;
    positX testSumQuire1AsPosit;
    testSumQuire1AsPosit.convert(testSumQuire1.to_value());

    positX acoefsiValue = -0.00828552; //Acoefs[i]
    positX xcoefsAcolsIValue = 0.000999451; //xcoefs[Acols[i]]


    std::cout << "testSumQuire1 == " << testSumQuire1 << std::endl;
    std::cout << "testSumQuire1AsPosit == " << testSumQuire1AsPosit << std::endl;
    for (int i = 0; i < 1000; i++) {
        testSumQuire1 += sw::unum::quire_mul(acoefsiValue, xcoefsAcolsIValue);
        std::cout << "testSumQuire1 after a += sw::unum::quire_mul() has run: " << testSumQuire1 << std::endl;
        testSumQuire1AsPosit.convert(testSumQuire1.to_value());
        std::cout << "testSumQuire1AsPosit == " << testSumQuire1AsPosit << std::endl;
    }




    return 0;
}