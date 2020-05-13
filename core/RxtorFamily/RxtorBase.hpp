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
        std::vector<float> __ChemMol;
        Eigen::MatrixXf __MainMat;
        std::vector<float> __ScalarVec;
        
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
                _RxnPtr(RxnPtr)
            {
                _setMainMat();
            }
            
            // StreamBase* 포인터를 이용함. 입/출력 스트림이 정의된 경우. 코멘트를 포함함.
            RxtorBase(StreamBase* inStreamPtr, StreamBase* outStreamPtr, RxnBase* RxnPtr, std::string& Comment):
                ProcObjBase(std::vector<StreamBase*>(1, inStreamPtr), std::vector<StreamBase*>(1, outStreamPtr), Comment), _RxnPtr(RxnPtr)
            {
                _setMainMat();
            }
            
            // 인스턴스 정의부

            // 반응기 내부의 물질들에 대해 총 반응지수를 계산함.
            float calcTotalQ() const
            {
                float ans = 1;
                int idx;
                const auto& RxnChemIdx = _RxnPtr->getChemIdx();
                auto effiVec = _RxnPtr->getEffiMat().col(_RxnPtr->getEffiMat().cols()-1);

                for (auto i = 0; i < RxnChemIdx.size(); ++i)
                {
                    idx = functions::getVecPos(__ChemIdx, RxnChemIdx[i]);
                    ans *= std::pow(__ChemMol[idx], effiVec[i]);
                }

                return ans;
            }

            // 반응기 내부의 물질들에 대해 반응지수를 계산함.
            std::vector<float> calcQ() const
            {
                int idx;
                const auto& RxnChemIdx = _RxnPtr->getChemIdx();
                const auto& EffiMat = _RxnPtr->getEffiMat();
                std::vector<float> ans(EffiMat.cols()-1, 1);
                Eigen::VectorXf effiVec;
                
                for (auto i = 0; i < EffiMat.cols()-1; ++i)
                {
                    effiVec = EffiMat.col(i);
                    for (auto j = 0; j < RxnChemIdx.size(); ++j)
                    {
                        idx = functions::getVecPos(__ChemIdx, RxnChemIdx[j]);
                        ans[i] *= std::pow(__ChemMol[idx], effiVec[j]);
                    }
                }

                return ans;
            }

            #ifdef _INCLUDE_CHEMPROCHELPER_SOLVER

            /*
            반응기에 대한 입력 스트림으로부터, 출력 스트림이 평형 상태가 되기 위한 전화율의 값을 각각 계산한다.
            */

           #endif

            #ifdef _INCLUDE_CHEMPROCHELPER_SOLVER

            /*
            값이 알려진 한쪽 스트림과 전화율을 바탕으로, 반대편 스트림의 값을 재설정함.
            단, 반응기가 정상 상태에서 작동한다고 가정한다.
            */
            void solveStreamFromConvRate()
            {
                bool direction;
                
                auto inStreamPtr = getInStreamIdx()[0];
                auto outStreamPtr = getOutStreamIdx()[0];
                const auto rxnChemIdx = _RxnPtr->getChemIdx();
                const auto rxnEffiMat = _RxnPtr->getEffiMat();
                const auto inChemIdx = inStreamPtr->getChemIdx();
                const auto outChemIdx = outStreamPtr->getChemIdx();

                if (inStreamPtr->chemMolIsAllKnown())
                {
                    if (!outStreamPtr->chemMolIsAllKnown()) direction = true;
                    else throw std::runtime_error("All stream has known.");
                }
                else
                {
                    if (outStreamPtr->chemMolIsAllKnown()) direction = false;
                    else throw std::runtime_error("All stream are unknown.");
                }

                Eigen::MatrixXf mat;
                mat.resize(rxnEffiMat.rows(), rxnEffiMat.cols()-1);
                mat.block(0,0,mat.rows(), mat.cols()) = rxnEffiMat.block(0,0,mat.rows(), mat.cols());

                Eigen::VectorXf convVec = Eigen::Map<Eigen::VectorXf>(__ScalarVec.data(), __ScalarVec.size());
                Eigen::VectorXf deltaVec = mat * convVec;
                std::unordered_map<ChemBase*, float> ChemMol;

                // 입력 스트림이 알려진 경우
                if (direction)
                {
                    const auto inChemMol = inStreamPtr->getChemMol();                    

                    for (auto chem : outChemIdx)
                    {
                        if (functions::inVector(rxnChemIdx, chem))
                        {
                            ChemMol[chem] = inStreamPtr->getChemMol(chem) + deltaVec[functions::getVecPos(rxnChemIdx, chem)];
                        }
                        else
                        {
                            ChemMol[chem] = inStreamPtr->getChemMol(chem);
                        }
                    }
                    
                    outStreamPtr->updateChem(ChemMol);
                }

                // 출력 스트림이 알려진 경우
                else
                {
                    const auto outChemMol = outStreamPtr->getChemMol();

                    for (auto chem : inChemIdx)
                    {
                        if (functions::inVector(rxnChemIdx, chem))
                        {
                            ChemMol[chem] = outStreamPtr->getChemMol(chem) - deltaVec[functions::getVecPos(rxnChemIdx, chem)];
                        }
                        else
                        {
                            ChemMol[chem] = outStreamPtr->getChemMol(chem);
                        }
                    }

                    inStreamPtr->updateChem(ChemMol);
                }
            }

            #endif

            #ifdef _INCLUDE_CHEMPROCHELPER_SOLVER

            /*
            입력 스트림과 출력 스트림의 값이 모두 알려진 경우 __ScalarVec의 값을 설정함.
            단, 반응기가 정상 상태에서 작동한다고 가정한다.
            */
            void solveConvRateFromStream()
            {
                const auto& rxnChemIdx = _RxnPtr->getChemIdx();
                
                Eigen::VectorXf deltaVec;
                deltaVec.resize(rxnChemIdx.size());
                deltaVec.setZero();                
                
                auto inStreamPtr = getInStreamIdx()[0];
                auto outStreamPtr = getOutStreamIdx()[0];

                int idx_rci, idx_str;

                for (auto i = 0; i < __ChemIdx.size(); ++i)
                {
                    auto ptr = __ChemIdx[i];
                    idx_rci = functions::getVecPos(rxnChemIdx, ptr);

                    // 반응물의 경우 빼야 함.
                    if (functions::inVector(inStreamPtr->getChemIdx(), ptr))
                    {
                        idx_str = functions::getVecPos(inStreamPtr->getChemIdx(), ptr);
                        deltaVec[idx_rci] -= inStreamPtr->getChemMol()[idx_str];
                    }

                    // 생성물의 경우 더해야 함.
                    if (functions::inVector(outStreamPtr->getChemIdx(), ptr))
                    {
                        idx_str = functions::getVecPos(outStreamPtr->getChemIdx(), ptr);
                        deltaVec[idx_rci] += outStreamPtr->getChemMol()[idx_str];
                    }
                }

                const auto& rxnEffiMat = _RxnPtr->getEffiMat();
                Eigen::MatrixXf mat;
                mat.resize(rxnEffiMat.rows(), rxnEffiMat.cols()-1);

                mat.block(0,0,rxnEffiMat.rows(), rxnEffiMat.cols()-1) = rxnEffiMat.block(0,0,rxnEffiMat.rows(), rxnEffiMat.cols()-1);
                Eigen::VectorXf res = mat.colPivHouseholderQr().solve(deltaVec);

                __ScalarVec.clear();
                __ScalarVec.resize(res.size());
                for (auto i = 0; i < res.size(); ++i)
                {
                    __ScalarVec[i] = res[i];
                }
            }

            #endif
    };
} // namespace chemprochelper

#endif