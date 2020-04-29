/*
ChemicalProcessHelper.hpp
----------------------------------------
복잡한 화학 공정에서의 계산을 물질 흐름을 중심으로 빠르고 편리하게 계산함.

본 헤더 파일은 thread-safe 하지 않음.

주요 클래스: ChemBase, RxnBase, StreamBase, ProcObjBase
----------------------------------------
ProcObjBase
    <= MixerBase, RxtorBase, SpliterBase
*/

#ifndef _CHEMPROCHELPER
#define _CHEMPROCHELPER

// 헤더 정의부
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <regex>
#include <algorithm>

// 이 헤더는 Eigen 3 라이브러리를 필수적으로 요구함.
#include <Eigen/Dense>

namespace chemprochelper
{
    // 주요 상수들 혹은 변수들을 포함함.
    namespace const_variables
    {
        const std::map<std::string, float> atomicMass = 
        {
            {"H", 1.007975},
            {"He", 4.002602},
            {"Li", 6.967499999999999},
            {"Be", 9.0121831},
            {"B", 10.8135},
            {"C", 12.0106},
            {"N", 14.006855},
            {"O", 15.9994},
            {"F", 18.998403163},
            {"Ne", 20.1797},
            {"Na", 22.98976928},
            {"Mg", 24.3055},
            {"Al", 26.9815385},
            {"Si", 28.085},
            {"P", 30.973761998},
            {"S", 32.067499999999995},
            {"Cl", 35.451499999999996},
            {"Ar", 39.948},
            {"K", 39.0983},
            {"Ca", 40.078},
            {"Sc", 44.955908},
            {"Ti", 47.867},
            {"V", 50.9415},
            {"Cr", 51.9961},
            {"Mn", 54.938044},
            {"Fe", 55.845},
            {"Co", 58.933194},
            {"Ni", 58.6934},
            {"Cu", 63.546},
            {"Zn", 65.38},
            {"Ga", 69.723},
            {"Ge", 72.63},
            {"As", 74.921595},
            {"Se", 78.971},
            {"Br", 79.904},
            {"Kr", 83.798},
            {"Rb", 85.4678},
            {"Sr", 87.62},
            {"Y", 88.90584},
            {"Zr", 91.224},
            {"Nb", 92.90637},
            {"Mo", 95.95},
            {"Tc", 98.0},
            {"Ru", 101.07},
            {"Rh", 102.9055},
            {"Pd", 106.42},
            {"Ag", 107.8682},
            {"Cd", 112.414},
            {"In", 114.818},
            {"Sn", 118.71},
            {"Sb", 121.76},
            {"Te", 127.6},
            {"I", 126.90447},
            {"Xe", 131.293},
            {"Cs", 132.90545196},
            {"Ba", 137.327},
            {"La", 138.90547},
            {"Ce", 140.116},
            {"Pr", 140.90766},
            {"Nd", 144.242},
            {"Pm", 145.0},
            {"Sm", 150.36},
            {"Eu", 151.964},
            {"Gd", 157.25},
            {"Tb", 158.92535},
            {"Dy", 162.5},
            {"Ho", 164.93033},
            {"Er", 167.259},
            {"Tm", 168.93422},
            {"Yb", 173.054},
            {"Lu", 174.9668},
            {"Hf", 178.49},
            {"Ta", 180.94788},
            {"W", 183.84},
            {"Re", 186.207},
            {"Os", 190.23},
            {"Ir", 192.217},
            {"Pt", 195.084},
            {"Au", 196.966569},
            {"Hg", 200.592},
            {"Tl", 204.3835},
            {"Pb", 207.2},
            {"Bi", 208.9804},
            {"Po", 209.0},
            {"At", 210.0},
            {"Rn", 222.0},
            {"Fr", 223.0},
            {"Ra", 226.0},
            {"Ac", 227.0},
            {"Th", 232.0377},
            {"Pa", 231.03588},
            {"U", 238.02891},
            {"Np", 237.0},
            {"Pu", 244.0},
            {"Am", 243.0},
            {"Cm", 247.0},
            {"Bk", 247.0},
            {"Cf", 251.0},
            {"Es", 252.0},
            {"Fm", 257.0},
            {"Md", 258.0},
            {"No", 259.0},
            {"Lr", 262.0},
            {"Rf", 267.0},
            {"Db", 268.0},
            {"Sg", 271.0},
            {"Bh", 272.0},
            {"Hs", 270.0},
            {"Mt", 276.0},
            {"Ds", 281.0},
            {"Rg", 280.0},
            {"Cn", 285.0},
            {"Nh", 284.0},
            {"Fl", 289.0},
            {"Mc", 288.0},
            {"Lv", 293.0},
            {"Ts", 292.0},
            {"Og", 295.0}
        };
        const std::regex pat_big("([0-9|.]{0,})([A-Z|a-z|0-9]{1,})");
        const std::regex pat_sml("([A-Z][a-z]?)(\\d{0,})");
    }
    
    // 주요 서브루틴들을 포함함.
    namespace functions
    {
        // std::sregex_iterator의 range-based for 구문 지원을 위한 클래스
        class _RegexIter
        {
            public:
                _RegexIter(const std::string& str, const std::regex& pat):
                    _pat(pat), _begin(str.begin(), str.end(), _pat), _end()
                {};
                auto begin()
                {
                    return _begin;
                }
                auto end()
                {
                    return _end;
                }
            private:
                std::regex _pat;
                std::sregex_iterator _begin;
                std::sregex_iterator _end;
        };

        // target의 vector 상의 인덱스를 반환함. 만약 target이 vector에 없는 경우 runtime error 발생.
        template<typename T>
        int getVecPos(const std::vector<T>& vector, const T& target)
        {
            auto it = std::find(vector.begin(), vector.end(), target);
            if (it == vector.end()) throw std::runtime_error("target is not in vector");
            return it - vector.begin();
        }

        // target이 vector에 있으면 true, 아니면 false를 반환한다. 위의 getVecPos와 동일한 구조를 가지고 있음.
        template<typename T>
        bool inVector(const std::vector<T>& vector, const T& target)
        {
            auto it = std::find(vector.begin(), vector.end(), target);
            if (it == vector.end()) return false;
            else return true;
        }

        // 분자량을 계산하여 반환함. 현재는 괄호는 처리할 수 없음.
        auto calMw(const std::string& eqn)
        {
            float mw = 0;
            std::string effi_str;
            std::string elem;
            int effi;
            for (auto m : _RegexIter(eqn, const_variables::pat_sml))
            {
                effi_str = m[2];
                elem = m[1];
                if (effi_str == "") effi = 1;
                else effi = std::stoi(effi_str);
                auto it = const_variables::atomicMass.find(elem);
                mw += effi * it->second;
            }
            return mw;
        }

        // n개 중에서 r개를 뽑는 모든 경우를 반환함.
        auto _combination(int n, int r)
        {
            assert(n >= r && r > 0);
            std::vector<std::vector<int>> output;
            std::vector<int> buff;
            if (r == 1)
            {
                for (int i = 0; i < n; ++i)
                {
                    buff.clear();
                    buff.push_back(i);
                    output.push_back(buff);
                }
                return output;
            }
            else if (n == r)
            {
                for (int i = 0; i < n; ++i) buff.push_back(i);
                output.push_back(buff);
                return output;
            }
            else
            {
                auto befList = _combination(n-1, r);
                for (auto& v : befList)
                {
                    output.push_back(v);
                }
                befList = _combination(n-1, r-1);
                for (auto& v : befList)
                {
                    buff.clear();
                    buff = v;
                    buff.push_back(n-1);
                    output.push_back(buff);
                }
                return output;
            }
        }
    }

    /*
    화합물을 구성하는 기본 클래스
    이 클래스를 이용해 코드를 구성할 때, 좌측값과 우측값의 특성에 반드시 주의해야 함.

    주의) 이 클래스는 몰질량 등의 정보를 포함하지 않은 클래스이다.
    */
    class ChemBase
    {
        private:

            // 화합물의 명칭을 저장함.
            std::string _Name;

            //화합물의 축약형을 저장함.
            std::string _Abb;

            // 화합물의 축약형과 그에 해당하는 ChemBase 객체의 포인터를 저장함.
            static std::unordered_map<std::string, ChemBase*> _AbbMap;

            // ChemBase 객체들의 포인터를 저장함(화합물을 숫자로 대응시키기 위함).
            static std::vector<ChemBase*> _PtrVec;

        public:

            // 정적 함수 정의부
            
            static auto getChemPtr(const int& i) {return _PtrVec[i];}
            static auto getChemPtr(const std::string Abb)
            {
                auto it = _AbbMap.find(Abb);
                if (it == _AbbMap.end()) throw std::runtime_error("Abb "+Abb+" isn't in ChemBase::_AbbMap");
                return it->second;
            }
            // Abb가 _AbbMap의 key로 있는 경우 true를 반환.
            static bool inAbbMap(const std::string& Abb)
            {
                auto it = _AbbMap.find(Abb);
                if (it == _AbbMap.end()) return false;
                else return true;
            }
            // ChemBasePtr이 ChemBase::_PtrVec의 몇 번째에 위치하는지 반환함.
            static auto getChemIdx(ChemBase* ChemBasePtr)            
            {
                return functions::getVecPos(_PtrVec, ChemBasePtr);
            }
            // Abb을 _Abb으로 가지는 ChemBase 객체의 포인터가 ChemBase::_PtrVec의 몇 번째에 위치하는지 반환함.
            static auto getChemIdx(const std::string Abb)           
            {
                auto it = _AbbMap.find(Abb);
                if (it == _AbbMap.end()) throw std::runtime_error("Abb "+Abb+" isn't in ChemBase::_AbbMap");
                return functions::getVecPos(_PtrVec, it->second);
            }
            
            // 생성자 정의부

            ChemBase() = default;
            ChemBase(const std::string& Abb):
                _Name(Abb), _Abb(_Abb)
            {
                auto it = _AbbMap.find(Abb);
                if (it != _AbbMap.end()) throw std::runtime_error("Abb "+Abb+" is already in ChemBase::_ChemList");
                _AbbMap[Abb] = this;
                _PtrVec.push_back(this);
            }
            ChemBase(const std::string& Name, const std::string& Abb):
                _Name(Name), _Abb(Abb)
            {
                auto it = _AbbMap.find(Abb);
                if (it != _AbbMap.end()) throw std::runtime_error("Abb "+Abb+" is already in ChemBase::_ChemList");
                _AbbMap[Abb] = this;
                _PtrVec.push_back(this);
            }
            
            // 소멸자 정의부

            ~ChemBase()
            {
                auto it = std::find(_PtrVec.begin(), _PtrVec.end(), this);
                if (it != _PtrVec.end())   // this가 ChemBase::_AbbMap이나 ChemBase::_PtrVec에 남은 경우 지움.
                {
                    _PtrVec.erase(it);
                    _AbbMap.erase(_Abb);
                }
            }

            // 대입 연산자 정의부

            // other이 우측값인 경우 대입 직후 other이 소멸하므로, ChemBase::_AbbMap과 ChemBase::_PtrVec에 등록된 포인터 주소를 변경함.
            ChemBase& operator=(ChemBase&& other)
            {
                _Name = other._Name;
                _Abb = other._Abb;

                _AbbMap[_Abb] = this;
                auto it = std::find(_PtrVec.begin(), _PtrVec.end(), &other);
                *it = this;

                return *this;
            }

            // other이 좌측값인 경우 대입 이후에도 other이 잔존하므로, ChemBase::_AbbMap과 ChemBase::_PtrVec을 변경하지 않음.
            ChemBase& operator=(ChemBase& other)            
            {
                _Name = other._Name;
                _Abb = other._Abb;

                return *this;
            }

            // getter 정의부

            auto getName() {return _Name;}
            auto getAbb() {return _Abb;}
    };
    std::unordered_map<std::string, ChemBase*> ChemBase::_AbbMap;
    std::vector<ChemBase*> ChemBase::_PtrVec;

    /*
    화학 반응식을 구성하는 기본 클래스.
    RxnBase::_EffiMat의 마지막 열은 현재 반응의 총 v(nu) 값과 동일하다.
    */
    class RxnBase
    {
        private:

            // 해당 반응식에 대한 간단한 메모를 할 수 있음.
            std::string _Comment = "";

            // 반응식에 포함된 화합물(ChemBase 객체)의 인덱스를 저장함.
            std::vector<int> _ChemIdx;

            // 반응식의 v(nu) 값을 저장하는 행렬.
            Eigen::MatrixXf _EffiMat;

            // RxnBase 객체들의 포인터를 저장함(반응식을 숫자로 대응시키기 위함).
            static std::vector<RxnBase*> _PtrVec;

            /*
            전달받은 화학식을 계수와 화합물의 std::vector로 분리함. effi에는 계수를,
            chem에는 화합물(ChemBase 객체)의 인덱스(ChemBase::_PtrVec 상 인덱스)를 저장함.
            direction = true이면 생성물, false이면 반응물로 생각함.
            */
            inline void _parseTerm(const std::string& term, const bool& direction,
                std::vector<float>& effiVec, std::vector<int>& chemVec)
            {
                int sgn;
                if (direction) sgn = -1;
                else sgn = 1;

                int chemIdx;

                for (auto m : functions::_RegexIter(term, const_variables::pat_big))
                {
                    if (m[1] == "") effiVec.push_back(sgn);
                    else effiVec.push_back(sgn*std::stof(m[1]));
                    chemIdx = ChemBase::getChemIdx(m[2]);

                    if (!functions::inVector(_ChemIdx, chemIdx)) _ChemIdx.push_back(chemIdx);
                    chemVec.push_back(chemIdx);
                }
            }

            // 전달받은 eqnVec을 바탕으로 this->_EffiMat을 구성함.
            void _setMat(const std::vector<std::string>& eqnVec)            
            {
                std::string reac, prod;
                std::vector<std::vector<float>> effiVec;
                std::vector<std::vector<int>> chemVec;
                int curEqnIdx = 0;
                int curChemIdx;

                for (auto eqn : eqnVec)
                {
                    auto strIdx = eqn.find("=");
                    if (strIdx == -1) throw std::runtime_error("Invalid chemical Reaction has entered.");

                    reac = eqn.substr(0, strIdx);
                    prod = eqn.substr(strIdx+1);

                    effiVec.push_back(std::vector<float>());
                    chemVec.push_back(std::vector<int>());

                    // 반응물 부분
                    _parseTerm(reac, true, effiVec[curEqnIdx], chemVec[curEqnIdx]);

                    // 생성물 부분
                    _parseTerm(prod, false, effiVec[curEqnIdx], chemVec[curEqnIdx]);

                    ++curEqnIdx;
                }

                // _EffiMat 초기화. 마지막 열에는 v(nu) 값의 총합을 입력해야 함.
                _EffiMat.resize(_ChemIdx.size(), curEqnIdx+1);
                _EffiMat.setZero();

                // ChemBase::_PtrVec상 인덱스와 _ChemIdx의 인덱스 간 매핑.
                std::unordered_map<int, int> chemIdxMap;
                for (auto i = 0; i < _ChemIdx.size(); ++i)
                {
                    chemIdxMap[_ChemIdx[i]] = i;
                }

                for (auto j = 0; j < curEqnIdx; ++j)
                {
                    auto curChemVec = chemVec[j];
                    auto curEffiVec = effiVec[j];

                    for (auto i = 0; i < curChemVec.size(); ++i)
                    {
                        _EffiMat(chemIdxMap[curChemVec[i]], j) = curEffiVec[i];
                    }
                }

                // 마지막 열은 v(nu)의 합으로 구성함.
                for (auto i = 0; i < _EffiMat.cols(); ++i)
                {
                    _EffiMat(i, curEqnIdx) = _EffiMat.row(i).sum();
                }
            }

        public:

            // 정적 함수 정의부

            static auto getRxnPtr(const int& i) {return _PtrVec[i];}

            // 생성자 정의부

            RxnBase() = default;
            RxnBase(const std::string& eqn)
            {
                _setMat({eqn});
                _PtrVec.push_back(this);
            }
            RxnBase(const std::string& eqn, const std::string& Comment):
                _Comment(Comment)
            {
                _setMat({eqn});
                _PtrVec.push_back(this);
            }
            RxnBase(const std::vector<std::string>& eqnVec)
            {
                _setMat(eqnVec);
                _PtrVec.push_back(this);
            }
            RxnBase(const std::vector<std::string>& eqnVec, const std::string& Comment):
                _Comment(Comment)
            {
                _setMat(eqnVec);
                _PtrVec.push_back(this);
            }

            // 소멸자 정의부

            ~RxnBase()
            {
                auto it = std::find(_PtrVec.begin(), _PtrVec.end(), this);
                if (it != _PtrVec.end())    // this가 RxnBase::_PtrVec에 남은 경우 삭제함.
                {
                    _PtrVec.erase(it);
                }
            }

            // 대입 연산자 정의부

            // other이 우측값인 경우 대입 직후 other이 소멸하므로, RxnBase::_PtrVec에 등록된 포인터 주소를 변경함.
            RxnBase& operator=(RxnBase&& other)            
            {
                _Comment = other._Comment;
                _ChemIdx = other._ChemIdx;
                _EffiMat = other._EffiMat;

                auto it = std::find(_PtrVec.begin(), _PtrVec.end(), &other);
                *it = this;

                return *this;
            }

            // other이 좌측값인 경우 대입 이후에도 other이 잔존하므로, RxnBase::_PtrVec을 변경하지 않음.
            RxnBase& operator=(RxnBase& other)            
            {
                _Comment = other._Comment;
                _ChemIdx = other._ChemIdx;
                _EffiMat = other._EffiMat;

                return *this;
            }

            // getter/setter 정의부

            auto getComment() {return _Comment;}
            auto getChemIdx() {return _ChemIdx;}
            auto getEffiMat() {return _EffiMat;}
            auto setComment(const std::string& Comment) {_Comment = Comment;}
    };
    std::vector<RxnBase*> RxnBase::_PtrVec;

    class StreamBase
    /*
    화학공정흐름도에서 물질의 흐름(flow stream)을 표현하는 클래스.
    */
    {
        private:

            // 흐름을 구성하는 화학종의 ChemBase::_PtrVec 상의 인덱스 값을 저장함.
            std::vector<int> _ChemIdx;

            // 화학종의 몰 유량이 알려져 있으면 true, 아닌 경우 false를 부여함.
            std::vector<bool> _ChemMask;

            // 화학종의 몰 유량을 저장함. 몰 유량을 알 수 없는 경우 0을 저장함.
            std::vector<float> _ChemMol;

            // StreamBase 객체들의 포인터를 저장함(흐름을 숫자로 대응시키기 위함).
            static std::vector<StreamBase*> _PtrVec;

            // StreamBase 객체에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            inline bool _updateChem(const int& ChemIdx, const bool& ChemMask, const float& ChemMol)            
            {
                if (functions::inVector(_ChemIdx, ChemIdx))
                {
                    auto idx = functions::getVecPos(_ChemIdx, ChemIdx);
                    _ChemMask[idx] = ChemMask;
                    _ChemMol[idx] = ChemMol;

                    return false;
                }
                else
                {
                    _ChemIdx.push_back(ChemIdx);
                    _ChemMask.push_back(ChemMask);
                    _ChemMol.push_back(ChemMol);

                    return true;
                }
            }

            // StreamBase 객체에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            inline bool _updateChem(const std::vector<int>& ChemIdx,
                const std::vector<bool>& ChemMask, const std::vector<float>& ChemMol)            
            {
                bool res = true;

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    if (functions::inVector(_ChemIdx, ChemIdx[i]))
                    {
                        res = false;

                        _ChemMask[i] = ChemMask[i];
                        _ChemMol[i] = ChemMol[i];
                    }
                    else
                    {
                        _ChemIdx.push_back(ChemIdx[i]);
                        _ChemMask.push_back(ChemMask[i]);
                        _ChemMol.push_back(ChemMol[i]);
                    }
                }

                return res;
            }

            // StreamBase 객체에서 화학종을 제거함. 성공한 경우 true를 반환함.
            inline bool _delChem(const int& ChemIdx)
            {
                auto it = std::find(_ChemIdx.begin(), _ChemIdx.end(), ChemIdx);
                if (it == _ChemIdx.end()) return false;

                auto idx = it - _ChemIdx.begin();
                _ChemIdx.erase(_ChemIdx.begin()+idx);
                _ChemMask.erase(_ChemMask.begin()+idx);
                _ChemMol.erase(_ChemMol.begin()+idx);

                return true;
            }

            // StreamBase 객체에서 특정 화학종을 미지수로 변경함. 성공한 경우 true를 반환함.
            inline bool _setChemUnkown(const int& ChemIdx)
            {
                auto it = std::find(_ChemIdx.begin(), _ChemIdx.end(), ChemIdx);
                if (it == _ChemIdx.end()) return false;

                auto idx = it - _ChemIdx.begin();
                _ChemMask[idx] = false;
                _ChemMol[idx] = 0;

                return true;
            }

        public:

            // 정적 함수 정의부

            static auto getStreamPtr(const int& i) {return _PtrVec[i];}
            
            // StreamasePtr이 StreamBase::_PtrVec의 몇 번째에 위치하는지 반환함.
            static auto getStreamIdx(StreamBase* StreamBasePtr)
            {
                return functions::getVecPos(_PtrVec, StreamBasePtr);
            }

            // 생성자 정의부

            StreamBase() = default;

            // ChemBase::_PtrVec 상의 인덱스를 이용함. 모든 물질의 몰 유량을 모르는 경우.
            StreamBase(const std::vector<int>& ChemIdx)            
            {
                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    _ChemIdx.push_back(ChemIdx[i]);
                    _ChemMask.push_back(false);
                    _ChemMol.push_back(0);
                }

                _PtrVec.push_back(this);
            }

            // ChemBase 객체를 직접 이용함. 모든 물질의 몰 유량을 모르는 경우.
            StreamBase(const std::vector<ChemBase*>& ChemVec)            
            {
                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    _ChemIdx.push_back(ChemBase::getChemIdx(ChemVec[i]));
                    _ChemMask.push_back(false);
                    _ChemMol.push_back(0);
                }

                _PtrVec.push_back(this);
            }

            // ChemBase::_PtrVec 상의 인덱스를 이용함. 모든 물질의 몰 유량을 아는 경우.
            StreamBase(const std::vector<int>& ChemIdx, const std::vector<float>& ChemMol)            
            {
                assert(ChemIdx.size() == ChemMol.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    _ChemIdx.push_back(ChemIdx[i]);
                    _ChemMask.push_back(true);
                    _ChemMol.push_back(ChemMol[i]);
                }

                _PtrVec.push_back(this);
            }
            // ChemBase 객체를 직접 이용함. 모든 물질의 몰 유량을 아는 경우.
            StreamBase(const std::vector<ChemBase*>& ChemVec, const std::vector<float>& ChemMol)            
            {
                assert(ChemVec.size() == ChemMol.size());

                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    _ChemIdx.push_back(ChemBase::getChemIdx(ChemVec[i]));
                    _ChemMask.push_back(true);
                    _ChemMol.push_back(ChemMol[i]);
                }

                _PtrVec.push_back(this);
            }
            
            // ChemBase::_PtrVec 상의 인덱스를 이용함. 일부 물질의 몰 유량만을 아는 경우.
            StreamBase(const std::vector<int>& ChemIdx, const std::vector<bool>& ChemMask,
                const std::vector<float>& ChemMol)
            {
                assert(ChemIdx.size() == ChemMask.size() && ChemMask.size() == ChemMol.size());

                _ChemIdx = ChemIdx;
                _ChemMask = ChemMask;
                _ChemMol = ChemMol;

                _PtrVec.push_back(this);
            }
            
            // ChemBase 객체를 직접 이용함. 일부 물질의 몰 유량만을 아는 경우.
            StreamBase(const std::vector<ChemBase*>& ChemVec, const std::vector<bool>& ChemMask,
                const std::vector<float>& ChemMol)
            {
                assert(ChemVec.size() == ChemMask.size() && ChemMask.size() == ChemMol.size());
                
                _ChemIdx.resize(ChemVec.size());
                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    _ChemIdx[i] = ChemBase::getChemIdx(ChemVec[i]);
                }

                _ChemMask = ChemMask;
                _ChemMol = ChemMol;

                _PtrVec.push_back(this);
            }
            
            // ChemBase::_PtrVec 상의 인덱스를 이용함. 일부 물질의 몰 유량만을 아는 경우.
            StreamBase(const std::vector<int>& ChemIdx,
                const std::unordered_map<int, float>& ChemMol)
            {
                _ChemIdx = ChemIdx;
                _ChemMask.resize(ChemIdx.size());
                _ChemMol.resize(ChemIdx.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    auto it = ChemMol.find(ChemIdx[i]);
                    if (it == ChemMol.end())
                    {
                        _ChemMask[i] = false;
                        _ChemMol[i] = 0;
                    }
                    else
                    {
                        _ChemMask[i] = true;
                        _ChemMol[i] = it->second;
                    }
                }

                _PtrVec.push_back(this);
            }
            
            // ChemBase 객체를 직접 이용함. 일부 물질의 몰 유량만을 아는 경우.
            StreamBase(const std::vector<ChemBase*>& ChemVec,
                const std::unordered_map<ChemBase*, float>& ChemMol)
            {
                _ChemIdx.resize(ChemVec.size());
                _ChemMask.resize(ChemVec.size());
                _ChemMol.resize(ChemVec.size());

                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    _ChemIdx[i] = ChemBase::getChemIdx(ChemVec[i]);

                    auto it = ChemMol.find(ChemVec[i]);
                    if (it == ChemMol.end())
                    {
                        _ChemMask[i] = false;
                        _ChemMol[i] = 0;
                    }
                    else
                    {
                        _ChemMask[i] = true;
                        _ChemMol[i] = it->second;
                    }
                }

                _PtrVec.push_back(this);
            }

            // 소멸자 정의부

            ~StreamBase()
            {
                auto it = std::find(_PtrVec.begin(), _PtrVec.end(), this);
                if (it != _PtrVec.end())    // this가 StreamBase::_PtrVec에 남은 경우 삭제함.
                {
                    _PtrVec.erase(it);
                }
            }

            // 대입 연산자 정의부
            
            // other이 우측값인 경우 대입 직후 other이 소멸하므로, StreamBase::_PtrVec에 등록된 포인터 주소를 변경함.
            StreamBase& operator=(StreamBase&& other)
            {
                _ChemIdx = other._ChemIdx;
                _ChemMask = other._ChemMask;
                _ChemMol = other._ChemMol;

                auto it = std::find(_PtrVec.begin(), _PtrVec.end(), &other);
                *it = this;

                return *this;
            }

            // other이 좌측값인 경우 대입 이후에도 other이 잔존하므로, StreamBase::_PtrVec을 변경하지 않음.
            StreamBase& operator=(StreamBase& other)
            {
                _ChemIdx = other._ChemIdx;
                _ChemMask = other._ChemMask;
                _ChemMol = other._ChemMol;

                return *this;
            }

            // setter/getter 정의부

            auto getChemIdx() {return _ChemIdx;}
            auto getChemMask() {return _ChemMask;}
            auto getChemMol() {return _ChemMol;}

            // 인스턴스 정의부

            // ChemBase::_PtrVec의 인덱스를 통해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const int& ChemIdx)            
            {
                return _updateChem(ChemIdx, false, 0);
            }

            // ChemBase* 포인터를 이용해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(ChemBase* ChemBasePtr)            
            {
                return _updateChem(ChemBase::getChemIdx(ChemBasePtr), false, 0);
            }

            // ChemBase::_PtrVec의 인덱스를 통해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<int>& ChemIdx)            
            {
                std::vector<bool> ChemMask(ChemIdx.size());
                std::vector<float> ChemMol(ChemIdx.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    ChemMask[i] = false;
                    ChemMol[i] = 0;
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase* 포인터를 이용해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<ChemBase*>& ChemVec)
            {
                std::vector<int> ChemIdx(ChemVec.size());
                std::vector<bool> ChemMask(ChemVec.size());
                std::vector<float> ChemMol(ChemVec.size());

                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    ChemIdx[i] = ChemBase::getChemIdx(ChemVec[i]);
                    ChemMask[i] = false;
                    ChemMol[i] = 0;
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase::_PtrVec의 인덱스를 통해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<int>& ChemIdx, const std::vector<float>& ChemMol)
            {
                assert(ChemIdx.size() == ChemMol.size());

                std::vector<bool> ChemMask(ChemIdx.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    ChemMask[i] = true;
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase* 포인터를 이용해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<ChemBase*>& ChemVec, const std::vector<float>& ChemMol)
            {
                assert(ChemVec.size() == ChemMol.size());

                std::vector<int> ChemIdx(ChemVec.size());
                std::vector<bool> ChemMask(ChemVec.size());

                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    ChemIdx[i] = ChemBase::getChemIdx(ChemVec[i]);
                    ChemMask[i] = true;
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase::_PtrVec의 인덱스를 통해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<int>& ChemIdx, const std::vector<bool>& ChemMask,
                const std::vector<float>& ChemMol)
            {
                assert(ChemIdx.size() == ChemMask.size() && ChemMask.size() == ChemMol.size());

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase* 포인터를 이용해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<ChemBase*>& ChemVec, const std::vector<bool>& ChemMask,
                const std::vector<float>& ChemMol)
            {
                assert(ChemVec.size() == ChemMask.size() && ChemMask.size() == ChemMol.size());

                std::vector<int> ChemIdx(ChemVec.size());

                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    ChemIdx[i] = ChemBase::getChemIdx(ChemVec[i]);
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase::_PtrVec의 인덱스를 통해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<int>& ChemIdx, std::unordered_map<int, float>& ChemMolMap)
            {
                std::vector<bool> ChemMask(ChemIdx.size());
                std::vector<float> ChemMol(ChemIdx.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    auto it = ChemMolMap.find(ChemIdx[i]);
                    if (it == ChemMolMap.end())
                    {
                        ChemMask[i] = false;
                        ChemMol[i] = 0;
                    }
                    else
                    {
                        ChemMask[i] = true;
                        ChemMol[i] = it->second;
                    }
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase* 포인터를 이용해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<ChemBase*>& ChemVec, std::unordered_map<int, float>& ChemMolMap)
            {
                std::vector<int> ChemIdx(ChemVec.size());
                std::vector<bool> ChemMask(ChemVec.size());
                std::vector<float> ChemMol(ChemVec.size());

                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    ChemIdx[i] = ChemBase::getChemIdx(ChemVec[i]);

                    auto it = ChemMolMap.find(ChemIdx[i]);
                    if (it == ChemMolMap.end())
                    {
                        ChemMask[i] = false;
                        ChemMol[i] = 0;
                    }
                    else
                    {
                        ChemMask[i] = true;
                        ChemMol[i] = it->second;
                    }
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase::_PtrVec의 인덱스를 통해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::unordered_map<int, float>& ChemMolMap)
            {
                std::vector<int> ChemIdx(ChemMolMap.size());
                std::vector<bool> ChemMask(ChemMolMap.size());
                std::vector<float> ChemMol(ChemMolMap.size());

                int idx = 0;
                for (auto it : ChemMolMap)
                {
                    ChemIdx[idx] = it.first;
                    ChemMask[idx] = true;
                    ChemMol[idx] = it.second;
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // ChemBase* 포인터를 이용해 스트림에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::unordered_map<ChemBase*, float>& ChemMolMap)
            {
                std::vector<int> ChemIdx(ChemMolMap.size());
                std::vector<bool> ChemMask(ChemMolMap.size());
                std::vector<float> ChemMol(ChemMolMap.size());

                int idx = 0;
                for (auto it : ChemMolMap)
                {
                    ChemIdx[idx] = ChemBase::getChemIdx(it.first);
                    ChemMask[idx] = true;
                    ChemMol[idx] = it.second;
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // StreamBase 객체에서 화학종을 제거함. 성공한 경우 true를 반환함.
            bool delChem(const int& ChemIdx)
            {
                return _delChem(ChemIdx);
            }

            // ChemBase* 포인터를 이용해 스트림에 화학종을 제거함. 성공한 경우 true를 반환함.
            bool delChem(ChemBase* ChemBasePtr)
            {
                return _delChem(ChemBase::getChemIdx(ChemBasePtr));
            }

            // StreamBase 객체에서 화학종을 제거함. 성공한 경우 true를 반환함.
            bool delChem(const std::vector<int>& ChemIdx)
            {
                bool res = true;

                for (auto idx : ChemIdx)
                {
                    if (!_delChem(idx)) res = false;
                }

                return res;
            }

            // ChemBase* 포인터를 이용해 스트림에 화학종을 제거함. 성공한 경우 true를 반환함.
            bool delChem(const std::vector<ChemBase*>& ChemVec)
            {
                bool res = true;

                for (auto ptr : ChemVec)
                {
                    if (!_delChem(ChemBase::getChemIdx(ptr))) res = false;
                }

                return res;
            }

            // StreamBase 객체에서 특정 화학종을 미지수로 변경함. 성공한 경우 true를 반환함.
            bool setChemUnknown(const int& ChemIdx)
            {
                return _setChemUnkown(ChemIdx);
            }

            // ChemBase* 포인터를 이용해 특정 화학종을 미지수로 변경함. 성공한 경우 true를 반환함.
            bool setChemUnknown(ChemBase* ChemBasePtr)
            {
                return _setChemUnkown(ChemBase::getChemIdx(ChemBasePtr));
            }

            // StreamBase 객체에서 특정 화학종을 미지수로 변경함. 성공한 경우 true를 반환함.
            bool setChemUnknown(const std::vector<int>& ChemIdx)
            {
                bool res = true;

                for (auto idx : ChemIdx)
                {
                    if (!_setChemUnkown(idx)) res = false;
                }

                return res;
            }

            // ChemBase* 포인터를 이용해 특정 화학종을 미지수로 변경함. 성공한 경우 true를 반환함.
            bool setChemUnknown(const std::vector<ChemBase*>& ChemVec)
            {
                bool res = true;

                for (auto ptr : ChemVec)
                {
                    if (!_setChemUnkown(ChemBase::getChemIdx(ptr))) res = false;
                }

                return res;
            }

            // StreamBase 객체에 특정 화학종이 스트림에 포함된 경우 true를 반환한다.
            bool inChemList(const int& ChemIdx)
            {
                auto it = std::find(_ChemIdx.begin(), _ChemIdx.end(), ChemIdx);
                if (it == _ChemIdx.end()) return false;
                else return true;
            }

            // ChemBase* 포인터를 이용해 특정 화학종이 스트림에 포함된 경우 true를 반환한다.
            bool inChemList(ChemBase* ChemBasePtr)
            {
                auto it = std::find(_ChemIdx.begin(), _ChemIdx.end(), ChemBase::getChemIdx(ChemBasePtr));
                if (it == _ChemIdx.end()) return false;
                else return true;
            }

            // StreamBase 객체에 특정 화학종이 스트림에 포함된 경우 true를 반환함.
            std::vector<bool> inChemList(const std::vector<int>& ChemIdx)
            {
                std::vector<bool> ChemMask(ChemIdx.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    ChemMask[i] = functions::inVector(_ChemIdx, ChemIdx[i]);
                }

                return ChemMask;
            }

            // ChemBase* 포인터를 이용해 특정 화학종이 스트림에 포함된 경우 true를 반환한다.
            std::vector<bool> inChemList(const std::vector<ChemBase*>& ChemVec)
            {
                std::vector<bool> ChemMask(ChemVec.size());

                for (auto i = 0; i < ChemVec.size(); ++i)
                {
                    ChemMask[i] = functions::inVector(_ChemIdx, ChemBase::getChemIdx(ChemVec[i]));
                }

                return ChemMask;
            }
    };
    std::vector<StreamBase*> StreamBase::_PtrVec;

    class ProcObjBase
    /*
    MixerBase, RxtorBase, SpliterBase의 상위 클래스.
    */
    {
        private:

            typedef std::vector<ChemBase*> ChemPtrVec;
            typedef std::vector<StreamBase*> StreamPtrVec;
            typedef std::vector<bool> BoolVec;
            typedef std::vector<float> FloatVec;
            typedef std::vector <int> IntVec;

            // _ProcObjNum은 0부터 시작함.
            int _ProcObjNum;

            // 다음 생성자 호출시 부여받는 _ProcObjNum
            static int _nextProcObjNum;

            // ProcObjBase 객체들의 포인터를 저장함.
            static std::vector<ProcObjBase*> _ProcObjList;
            
            StreamPtrVec _inStream;
            StreamPtrVec _outStream;

            inline void _setProcObjNum()
            // _ProcObjNum을 할당함.
            {
                _ProcObjNum = _nextProcObjNum++;
                _ProcObjList.push_back(this);
            }

            inline bool _addStream(StreamBase* Stream, const bool& direction)
            // Stream을 추가함. direction = true이면 in, false이면 out이다. 성공시 true 반환.
            {
                if (inStreamList(Stream, direction)) return false;

                if (direction) _inStream.push_back(Stream);
                else _outStream.push_back(Stream);

                return true;
            }

            inline bool _delStream(StreamBase* Stream, const bool& direction)
            // Stream을 제거함. direction = true이면 in, false이면 out이다. 성공시 true 반환.
            {
                if (!inStreamList(Stream, direction)) return false;

                int idx;
                if (direction)
                {
                    idx = functions::getPosition(_inStream, Stream);
                    _inStream.erase(_inStream.begin()+idx);
                }
                else
                {
                    idx = functions::getPosition(_outStream, Stream);
                    _outStream.erase(_outStream.begin()+idx);
                }
                
                return true;
            }
        
        protected:

            // protected 항목에 대해서는 자녀 클래스에서 초기화함.

            std::string __Name;
            Eigen::MatrixXf __MainMat;
            FloatVec __ScalarVec;

            ChemPtrVec __ChemList;
            BoolVec __ChemIsKnown;
            FloatVec __ChemMol;

        public:

            static auto getProcObjList() {return _ProcObjList;}
            static auto getProcObjPtr(const int& i) {return _ProcObjList[i];}

            // 생성자 선언부

            ProcObjBase()
            {
                _setProcObjNum();
            }
            ProcObjBase(const StreamPtrVec& inStream, const StreamPtrVec& outStream):
                _inStream(inStream), _outStream(outStream)
            {
                _setProcObjNum();
            }
            ProcObjBase(const IntVec& inStreamNum, const IntVec& outStreamNum)
            {
                StreamPtrVec inStream(inStreamNum.size());
                StreamPtrVec outStream(outStreamNum.size());
                for (auto i : inStreamNum) inStream[i] = StreamBase::getStreamPtr(i);
                for (auto o : outStreamNum) outStream[o] = StreamBase::getStreamPtr(o);

                _inStream = inStream;
                _outStream = outStream;
            }

            // setter/getter 정의부

            auto getName() {return __Name;}
            auto getInStream() {return _inStream;}
            auto getOutStream() {return _outStream;}
            auto getMainMat() {return __MainMat;}
            auto getScalarVec() {return __ScalarVec;}
            
            bool addStream(StreamBase* Stream, const bool& direction)
            // direction = true이면 in, false이면 out. 추가에 성공하면 true, 아닌 경우 false 반환.
            {
                return _addStream(Stream, direction);
            }
            bool addStream(const int& StreamNum, const bool& direction)
            // StreamBase::_StreamNum을 사용함. direction = true이면 in, false이면 out. 추가에 성공하면 true, 아닌 경우 false 반환.
            {
                return _addStream(StreamBase::getStreamPtr(StreamNum), direction);             
            }
            bool addStream(const StreamPtrVec& StreamList, const BoolVec& direction)
            // direction = true이면 in, false이면 out. 추가에 성공하면 true, 아닌 경우 false 반환.
            {
                assert(StreamList.size() == direction.size());

                bool res = true;
                for (auto i = 0; i < StreamList.size(); ++i)
                {
                    if (!_addStream(StreamList[i], direction[i])) res = false;
                }

                return res;
            }
            bool addStream(const IntVec& StreamNumList, const BoolVec& direction)
            // StreamBase::_StreamNum을 사용함. direction = true이면 in, false이면 out. 추가에 성공하면 true, 아닌 경우 false 반환.
            {
                assert(StreamNumList.size() == direction.size());

                bool res = true;
                for (auto i = 0; i < StreamNumList.size(); ++i)
                {
                    if (!_addStream(StreamBase::getStreamPtr(StreamNumList[i]), direction[i])) res = false;
                }

                return res;
            }

            bool delStream(StreamBase* Stream, const bool& direction)
            // direction = true이면 in, false이면 out. 제거에 성공하면 ture, 아닌 경우 false 반환.
            {   
                return _delStream(Stream, direction);
            }
            bool delStream(const int& StreamNum, const bool& direction)
            // StreamBase::_StreamNum을 사용함. direction = true이면 in, false이면 out. 제거에 성공하면 ture, 아닌 경우 false 반환.
            {
                return _delStream(StreamBase::getStreamPtr(StreamNum), direction);
            }
            bool delStream(const StreamPtrVec& StreamList, const BoolVec& direction)
            // direction = true이면 in, false이면 out. 제거에 성공하면 ture, 아닌 경우 false 반환.
            {
                assert(StreamList.size() == direction.size());

                bool res = true;
                for (auto i = 0; i < StreamList.size(); ++i)
                {
                    if (!_delStream(StreamList[i], direction[i])) res = false;
                }

                return false;
            }
            bool delStream(const IntVec& StreamNumList, const BoolVec& direction)
            // StreamBase::_StreamNum을 사용함. direction = true이면 in, false이면 out. 제거에 성공하면 ture, 아닌 경우 false 반환.
            {
                assert(StreamNumList.size() == direction.size());

                bool res = true;
                for (auto i = 0; i < StreamNumList.size(); ++i)
                {
                    if (!_delStream(StreamBase::getStreamPtr(StreamNumList[i]), direction[i])) res = false;
                }

                return res;
            }
            // 인스턴스 정의부

            bool inStreamList(StreamBase* Stream, bool direction)
            // Stream을 확인함. direction = true이면 in, false이면 out이다. 이미 존재하는 경우 true 반환.
            {   
                StreamPtrVec::const_iterator it;
                if (direction)
                {
                    it = std::find(_inStream.begin(), _inStream.end(), Stream);
                    if (it == _inStream.end()) return false;
                }
                else
                {
                    it = std::find(_outStream.begin(), _outStream.end(), Stream);
                    if (it == _outStream.end()) return false;
                }

                return true;
            }
    };
    int ProcObjBase::_nextProcObjNum = 0;
    std::vector<ProcObjBase*> ProcObjBase::_ProcObjList;

    class RxtorBase : public ProcObjBase
    /*
    화학 반응기를 나타내는 클래스. ProcObjBase로부터 상속됨.
    ProcObjBase에서 화학 반응식 및 몰수지 부분을 추가함.
    이 클래스는 Eigen 3 라이브러리의 Solver 기능을 사용함.
    INTEL(R) Math Kernel Library, OpenBLAS 등이 있는 경우 CMake를 이용해 속도 향상이 가능함.
    */
    {   
        private:

            typedef std::vector<ChemBase*> ChemPtrVec;
            typedef std::vector<StreamBase*> StreamPtrVec;
            typedef std::vector<bool> BoolVec;
            typedef std::vector<float> FloatVec;
            typedef std::vector <int> IntVec;

            /*
            ProcObjBase로부터 다음과 같은 protected 변수들을 상속받음.
            std::string __Name;
            Eigen::MatrixXf __MainMat;
            FloatVec __ScalarVec;
            ------------------------------
            ChemPtrVec __ChemList;
            BoolVec __ChemIsKnown;
            FloatVec __ChemMol;
            */

            // _RxtorNum은 0부터 시작함.
            int _RxtorNum;

            // 다음 생성자 호출시 부여받는 _RxtorNum
            static int _nextRxtorNum;

            // RxtorBase 객체들의 포인터를 저장함.
            static std::vector<RxtorBase*> _RxtorList;
            

            inline void _setRxtorNum()
            {
                _RxtorNum = _nextRxtorNum++;
                _RxtorList.push_back(this);
            }

        public:

            // 다음 생성자 호출시 부여받는  _rxtorNum.
            static int nextRxtorNum;

            // Rxtorase 객체들의 _rxtorNum과 포인터를 저장함.
            static std::map<int, RxtorBase*> RxtorMap;

            // 생성자 선언부

            
            //인스턴스 선언부

            bool checkStreamNum()
            // 입력 스트림과 출력 스트림이 모두 1개뿐인지 확인함.
            {

            }

            bool makeChemList()
            // 입력 스트림과 출력 스트림의 물질들을 모두 저장한 _chemList를 생성함. 성공한 경우 true, 실패한 경우 false 반환.
            {

                // to be continued.
            }

            
    };
}

#endif