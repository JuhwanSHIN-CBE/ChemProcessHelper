/*
core/RxtorFamily/RxtorBase.hpp
------------------------------
화학 반응기의 기초가 되는 RxtorBase 클래스를 정의함.
*/
#ifndef _CHEMPROCHELPER_RXTORBASE
#define _CHEMPROCHELPER_RXTORBASE

namespace chemprochelper
{
    /*
    화학 반응기를 지정하는 기본 클래스
    --------------------------------
    RxtorBase는 다음과 같은 멤버 변수를 가짐.
    private:
        _RxnPtr : 화학 반응식을 나타내는 RxnBase 객체의 포인터를 저장함.
    */
    class RxtorBase : public ProcObjBase
    {
        /*
        ProcObjBase로부터,

        std::vector<ChemBase*> __ChemIdx;
        Eigen::MatrixXf __MainMat
        std::vector<float> __ScalarVec
        
        세 변수를 상속받음(모두 protected)
        */
        private:

            // 화학 반응식을 나타내는 RxnBase 객체의 포인터를 저장함.
            RxnBase* _RxnPtr;

            /*
            RxnBase 객체로부터 반응기의 행렬(Murphy의 화학공정계산 p.210 참조)을 __MainMat에 구성함.
            행렬을 구성하는데 성공하면 true, 실패하면 false를 반환함.
            */
            void _setMainMat()
            {
                StreamBase* inStreamPtr = getInStreamIdx()[0];
                StreamBase* outStreamPtr = getOutStreamIdx()[0];
                
                auto RxnEffiMat = _RxnPtr->getEffiMat();
                auto RxnChemIdx = _RxnPtr->getChemIdx();
                auto inStrChemIdx = inStreamPtr->getChemIdx();

                // _ChemIdx에 입력 스트림과 출력 스트림의 화합물들을 다 저장함.
                __ChemIdx = inStreamPtr->getChemIdx();
                for (auto idx : outStreamPtr->getChemIdx())
                {
                    if (!functions::inVector(__ChemIdx, idx)) __ChemIdx.push_back(idx);
                }                
                for (auto idx : RxnChemIdx)
                {
                    if (!functions::inVector(__ChemIdx, idx)) throw std::runtime_error("StreamBase object can't cover RxnBase object");
                }

                // __MainMat의 크기를 맞춤.
                int rows = __ChemIdx.size() + inStrChemIdx.size();
                int cols = __ChemIdx.size() + RxnEffiMat.cols() - 1;

                __MainMat.resize(rows, cols);
                __MainMat.setZero();

                for (auto i = 0; i < __ChemIdx.size(); ++i)
                {
                    auto idx = __ChemIdx[i];

                    // 단위행렬 부분의 값을 설정함.
                    __MainMat(i, i) = 1;

                    // 우측 상단 부분의 값을 결정함.
                    if (functions::inVector(RxnChemIdx, idx))   // 반응에 참여하는 화학종의 경우
                    {
                        for (auto j = 0; j < RxnEffiMat.cols() - 1; ++j)
                        {
                            __MainMat(i, __ChemIdx.size() + j) = -1 * RxnEffiMat(functions::getVecPos(RxnChemIdx, idx), j);
                        }
                    }                  
                }

                // 우측 하단 부분의 값을 결정함.
                for (auto i = 0; i < inStrChemIdx.size(); ++i)
                {
                    auto idx = inStrChemIdx[i];

                    if (functions::inVector(RxnChemIdx, idx))   // 반응에 참여하는 경우
                    {
                        for (auto j = 0; j < RxnEffiMat.cols() - 1; ++j)
                        {
                            __MainMat(__ChemIdx.size() + i, __ChemIdx.size() + j) = -1 * RxnEffiMat(functions::getVecPos(RxnChemIdx, idx), j);
                        }
                    }
                }
            }

        public:

            // 생성자 정의부

            // 임시 객체를 위한 생성자
            RxtorBase():
                ProcObjBase() {}
            
            // StreamBase* 포인터를 이용함. 입/출력 스트림이 정의된 경우
            RxtorBase(StreamBase* inStreamPtr, StreamBase* outStreamPtr, RxnBase* RxnPtr):
                ProcObjBase(std::vector<StreamBase*>(1, inStreamPtr), std::vector<StreamBase*>(1, outStreamPtr)),
                _RxnPtr(RxnPtr) {}
            
            // StreamBase* 포인터를 이용함. 입/출력 스트림이 정의된 경우. 코멘트를 포함함.
            RxtorBase(StreamBase* inStreamPtr, StreamBase* outStreamPtr, RxnBase* RxnPtr, std::string& Comment):
                ProcObjBase(std::vector<StreamBase*>(1, inStreamPtr), std::vector<StreamBase*>(1, outStreamPtr), Comment), _RxnPtr(RxnPtr) {}
    };
} // namespace chemprochelper

#endif