#define _INCLUDE_CHEMPROCHELPER_SOLVER
#include "ChemProcHelper.hpp"

int main()
{
    std::string input;
    std::getline(std::cin, input);
    std::cout << chemprochelper::functions::_balRxnEqn(input) << std::endl;
}