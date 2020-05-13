#define _INCLUDE_CHEMPROCHELPER_SOLVER
#include "ChemProcHelper.hpp"

#include "mpl.h"

#include <boost/numeric/odeint.hpp>


int main()
{
    std::vector<chemprochelper::ChemBase> ChemBaseVec(5);
    std::vector<std::string> ChemStrVec = {"NH3", "CH3OH", "H2O", "CH3NH2", "(CH3)2NH"};
    for (auto i = 0; i < ChemStrVec.size(); ++i)
    {
        ChemBaseVec[i] = chemprochelper::ChemBase(ChemStrVec[i]);
    }

    std::vector<chemprochelper::RxnBase> RxnBaseVec;
    std::vector<std::string> RxnStrVec = {"NH3 + CH3OH = CH3NH2 + H2O", "CH3NH2 + CH3OH = (CH3)2NH + H2O"};
    RxnBaseVec.push_back(chemprochelper::RxnBase(RxnStrVec));

    std::vector<chemprochelper::StreamBase> StreamBaseVec;
    std::vector<chemprochelper::ChemBase*> inputChemBaseVec;
    std::vector<float> inputChemMolVec;

    inputChemBaseVec = {&ChemBaseVec[0], &ChemBaseVec[1]};
    inputChemMolVec = {100.0, 100.0};
    StreamBaseVec.push_back(
        chemprochelper::StreamBase(
            inputChemBaseVec, inputChemMolVec
        )
    );

    inputChemBaseVec.resize(ChemBaseVec.size());
    for (auto i = 0; i < inputChemBaseVec.size(); ++i) inputChemBaseVec[i] = &ChemBaseVec[i];
    inputChemMolVec = {20.0 * 200/100, 9.4 * 200/100, 41.1 * 200/100 ,18.6 * 200/100, 10.9 * 200/100};
    StreamBaseVec.push_back(
        chemprochelper::StreamBase(
            inputChemBaseVec, inputChemMolVec
        )
    );

    std::vector<chemprochelper::RxtorBase> RxtorBaseVec;
    RxtorBaseVec.push_back(
        chemprochelper::RxtorBase(
            &StreamBaseVec[0], &StreamBaseVec[1], &RxnBaseVec[0]
        )
    );

    RxtorBaseVec[0].solveConvRateFromStream();

    const auto& ScalarVec = RxtorBaseVec[0].getScalarVec();
    for (auto i = 0; i < ScalarVec.size(); ++i)
    {
        std::cout << ScalarVec[i] << std::endl;
    }


    StreamBaseVec[1].setAllUnknown();
    RxtorBaseVec[0].solveStreamFromConvRate();


    for (auto i = 0; i < StreamBaseVec[1].getChemIdx().size(); ++i)
    {
        std::cout << StreamBaseVec[1].getChemIdx()[i]->getAbb() << '\t' << StreamBaseVec[1].getChemMol()[i] << std::endl;
    }

}