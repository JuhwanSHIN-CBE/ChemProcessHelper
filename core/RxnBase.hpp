namespace chemprochelper
{
    /*
    화학 반응식을 구성하는 기본 클래스.
    RxnBase::_EffiMat의 마지막 열은 현재 반응의 총 v(nu) 값과 동일하다.
    */
    class RxnBase
    {
        private:

            // 해당 반응식에 대한 간단한 메모를 할 수 있음.
            std::string _Comment = "";

            // 반응식에 포함된 화합물(ChemBase 객체)의 포인터를 저장함.
            std::vector<ChemBase*> _ChemIdx;

            // 반응식의 v(nu) 값을 저장하는 행렬.
            Eigen::MatrixXf _EffiMat;

            /*
            전달받은 화학식을 계수와 화합물의 std::vector로 분리함. effi에는 계수를,
            chem에는 화합물(ChemBase 객체)의 인덱스(ChemBase::_PtrVec 상 인덱스)를 저장함.
            direction = true이면 생성물, false이면 반응물로 생각함.
            */
            void _parseTerm(const std::string& term, const bool& direction,
                std::vector<float>& effiVec, std::vector<ChemBase*>& chemVec)
            {
                int sgn = direction ? 1 : -1;

                ChemBase* chemIdx;

                for (auto m : functions::_RegexIter(term, const_variables::pat_big))
                {
                    (m[1] == "") ? effiVec.push_back(sgn) : effiVec.push_back(sgn*std::stof(m[1]));
                    chemIdx = ChemBase::getChemPtr(m[2]);

                    if (!functions::inVector(_ChemIdx, chemIdx)) _ChemIdx.push_back(chemIdx);
                    chemVec.push_back(chemIdx);
                }
            }

            // 전달받은 eqnVec을 바탕으로 this->_EffiMat을 구성함.
            void _setMat(const std::vector<std::string>& eqnVec)
            {
                std::string reac, prod;
                std::vector<std::vector<float>> effiVec;
                std::vector<std::vector<ChemBase*>> chemVec;
                int curEqnIdx = 0;
                int curChemIdx;

                for (auto& eqn : eqnVec)
                {
                    auto strIdx = eqn.find("=");
                    if (strIdx == -1) throw std::runtime_error("Invalid chemical reaction has entered.");

                    reac = eqn.substr(0, strIdx);
                    prod = eqn.substr(strIdx+1);

                    effiVec.push_back(std::vector<float>());
                    chemVec.push_back(std::vector<ChemBase*>());

                    // 반응물 부분
                    _parseTerm(reac, true, effiVec[curEqnIdx], chemVec[curEqnIdx]);

                    // 생성물 부분
                    _parseTerm(prod, false, effiVec[curEqnIdx], chemVec[curEqnIdx]);

                    ++curEqnIdx;
                }

                // 마지막 열은 v(nu)의 합으로 구성함.
                for (auto i = 0; i < _EffiMat.rows(); ++i)
                {
                    _EffiMat(i, curEqnIdx) = _EffiMat.row(i).sum();
                }
            }

        public:

            // 생성자 정의부

            // 디폴트 생성자
            RxnBase() = default;

            RxnBase(const std::string& eqn)
            {
                _setMat({eqn});
            }

            RxnBase(const std::string& eqn, const std::string& Comment):
                _Comment(Comment)
            {
                _setMat({eqn});
            }

            RxnBase(const std::vector<std::string>& eqnVec)
            {
                _setMat(eqnVec);
            }
            
            RxnBase(const std::vector<std::string>& eqnVec, const std::string& Comment):
                _Comment(Comment)
            {
                _setMat(eqnVec);
            }

            // getter 정의부
            auto getComment() {return _Comment;}
            auto getChemIdx() {return _ChemIdx;}
            auto getEffiMat() {return _EffiMat;}
    };
} // namespace chemprochelper