/*
ChemicalProcessHelper.hpp
복잡한 화학 공정에서의 계산을 물질 흐름을 중심으로 빠르고 편리하게 계산함.

본 헤더 파일은 멀티 쓰레딩을 염두에 두고 설계한 것이 아니므로 thread-safe한지 알 수 없음.

주요 클래스: ChemBase, ProcObjBase, MixerBase, RxtorBase, SpliterBase, StreamBase
ProcObjBase
    <= MixerBase, RxtorBase, SpliterBase
*/

#ifndef _CHEMPROCHELPER
#define _CHEMPROCHELPER

// 헤더 정의부
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <regex>
#include <algorithm>

// 이 헤더는 Eigen 3 라이브러리를 필수적으로 요구함.
#include <Eigen/Dense>

namespace chemprochelper
{
    namespace const_variables
    // 주요 상수들 혹은 변수들을 포함함.
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
        const std::regex pat_big("~([0-9|.]{0,})([A-Z|a-z|0-9]{1,})");
        const std::regex pat_sml("([A-Z][a-z]?)(\\d{0,})");
    }

    namespace functions
    // 주요 서브루틴들을 포함함.
    {
        class _RegexIter
        // std::sregex_iterator의 range-based for 구문 지원을 위한 클래스
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

        template<typename T>
        int getPosition(const std::vector<T>& vector, const T& target)
        // target의 vector 상의 인덱스를 반환함. 만약 target이 vector에 없는 경우 runtime error 발생.
        {
            auto it = std::find(vector.begin(), vector.end(), target);
            if (it == vector.end()) throw std::runtime_error("target is not in vector");
            return it - vector.begin();
        }

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

    class ChemBase
    /*
    화합물을 구성하는 기본 클래스

    사용 예시
    1) ChemBase("H2O");
    2) ChemBase("Water", "H2O");
    3) ChemBase("Water", "H2O", 18.0);
    4) ChemBase("Water", 18.0);
    */
    {
        private:

            // 다음 생성자 호출시 부여받는 _chemNum.
            static int nextChemNum;

            // ChemBase 객체들의 포인터를 저장함.
            static std::vector<ChemBase*> ChemList;

            inline void _setChemNum()
            // _chemNum을 할당함.
            {
                _ChemNum = nextChemNum++;
                ChemList.push_back(this);
            };

            int _ChemNum;
            std::string _Name;
            std::string _Abb;
            float _Mw;

        public:

            static auto getChemList() {return ChemList;}
            static auto getChemPtr(const int& i) {return ChemList[i];}

            // 생성자 정의부

            ChemBase() {_setChemNum();};
            
            ChemBase(const std::string& Abb):
                _Name(Abb), _Abb(Abb), _Mw(functions::calMw(_Abb)) {_setChemNum();};
            ChemBase(const std::string& Name, const std::string& Abb):
                _Name(Name), _Abb(Abb), _Mw(functions::calMw(_Abb)) {_setChemNum();};
            ChemBase(const std::string& Name, const float& Mw):
                _Name(Name), _Abb(""), _Mw(Mw) {_setChemNum();};
            ChemBase(const std::string& Name, const std::string& Abb, const float& Mw):
                _Name(Name), _Abb(Abb), _Mw(Mw) {_setChemNum();};

            // getter 정의부
            
            auto getName() {return _Name;}
            auto getAbb() {return _Abb;}
            auto getMw() {return _Mw;}
            auto getChemNum() {return _ChemNum;}
    };
    int ChemBase::nextChemNum = 0;

    class StreamBase
    /*
    화학공정흐름도에서 물질의 흐름(flow stream)을 표현하는 클래스.
    */
    {
        private:

            typedef std::vector<ChemBase*> ChemPtrVec;
            typedef std::vector<float> FloatVec;
            typedef std::vector<bool> BoolVec;
            typedef std::vector<int> IntVec;

            // _StreamNum은 0부터 시작함.
            int _StreamNum;

            // 다음 생성자 호출시 부여받는 _StreamNum.
            static int nextStreamNum;

            // StreamBase 객체들의 포인터를 저장함.
            static std::vector<StreamBase*> StreamList;

            ChemPtrVec _ChemList;

            // 해당 화학종의 몰수가 알려져 있으면 true, 아니면 false 값을 부여함.
            BoolVec _ChemIsKnown;

            // _ChemIsKnown[i] = false이면 _ChemMol[i] = 0.
            FloatVec _ChemMol;


            inline void _setStreamNum()
            // _streamNum을 설정함.
            {
                _StreamNum = StreamBase::nextStreamNum++;
                StreamList.push_back(this);
            }

            inline void _setInnerList(const ChemPtrVec& ChemList,
                const BoolVec& ChemIsKnown, const FloatVec& ChemMol)
            // _ChemList, _ChemMol, _ChemIsKnown을 설정함.
            {
                _ChemList = ChemList;
                _ChemIsKnown = ChemIsKnown;
                _ChemMol = ChemMol;
            }

            inline bool _addChem(ChemBase* Chem, const bool& IsKnown, const float& Mol)
            // StreamBase 객체에 화학종을 추가함. 성공시 true를 반환함.
            {
                if (inChemList(Chem)) return false;

                _ChemList.push_back(Chem);
                _ChemIsKnown.push_back(IsKnown);
                _ChemMol.push_back(Mol);

                return true;
            }

            inline bool _setChemMol(ChemBase* Chem, const float& Mol)
            // 화학종 Chem의 몰 유량을 Mol로 설정함. 성공시 true를 반환함.
            {
                if (!inChemList(Chem)) return false;
                int idx(functions::getPosition(_ChemList, Chem));

                _ChemIsKnown[idx] = true;
                _ChemMol[idx] = Mol;

                return true;
            }

            inline bool _setChemUnknown(ChemBase* Chem)
            // 화학종 Chem의 몰 유량을 미지수로 변경함. 성공시 true를 반환함.
            {
                if (!inChemList(Chem)) return false;
                int idx(functions::getPosition(_ChemList, Chem));

                _ChemIsKnown[idx] = false;
                _ChemMol[idx] = 0;

                return true;
            }

            inline bool _delChem(ChemBase* Chem)
            // StreamBase 객체에 i번째 인덱스의 화학종을 제거함. 성공시 true를 반환함.
            {
                if (!inChemList(Chem)) return false;
                int idx(functions::getPosition(_ChemList, Chem));

                _ChemList.erase(_ChemList.begin()+idx);
                _ChemIsKnown.erase(_ChemIsKnown.begin()+idx);
                _ChemMol.erase(_ChemMol.begin()+idx);

                return true;
            }

        public:

            static auto getStreamList() {return StreamList;}
            static auto getStreamPtr(const int& i) {return StreamList[i];}
            
            // 생성자 정의부

            StreamBase()
            {
                _setStreamNum();
            }            
            StreamBase(const ChemPtrVec& ChemList)
            {
                BoolVec ChemIsKnown(ChemList.size());
                FloatVec ChemMol(ChemList.size());
                for (auto i = 0; i < ChemList.size(); ++i)
                {
                    ChemIsKnown[i] = false;
                    ChemMol[i] = 0;
                }

                _setInnerList(ChemList, ChemIsKnown, ChemMol);
                _setStreamNum();
            }
            StreamBase(const IntVec& ChemNumList)
            // ChemBase::_chemNum을 활용한 생성자. 
            {
                ChemPtrVec ChemList(ChemNumList.size());
                BoolVec ChemIsKnown(ChemNumList.size());
                FloatVec ChemMol(ChemNumList.size());
                for (auto i = 0; i < ChemNumList.size(); ++i)
                {
                    ChemIsKnown[i] = false;
                    ChemMol[i] = 0;
                    ChemList[i] = ChemBase::getChemPtr(i);
                }

                _setInnerList(ChemList, ChemIsKnown, ChemMol);
                _setStreamNum();
            }
            StreamBase(const ChemPtrVec& ChemList, const FloatVec& ChemMol)
            // 모든 화학종의 몰 유량을 알고 있는 경우 사용함.
            {
                assert(ChemList.size() == ChemMol.size());

                BoolVec ChemIsKnown(ChemList.size());
                for (auto i = 0; i < ChemList.size(); ++i)
                {
                    ChemIsKnown[i] = true;
                }

                _setInnerList(ChemList, ChemIsKnown, ChemMol);
                _setStreamNum();
            }
            StreamBase(const IntVec& ChemNumList, const FloatVec& ChemMol)
            // ChemBase::_chemNum을 활용한 생성자. 모든 화학종의 몰 유량을 알고 있는 경우 사용함.
            {
                assert(ChemNumList.size() == ChemMol.size());

                BoolVec ChemIsKnown(ChemNumList.size());
                ChemPtrVec ChemList(ChemNumList.size());
                for (auto i = 0; i < ChemNumList.size(); ++i)
                {
                    ChemIsKnown[i] = true;
                    ChemList[i] = ChemBase::getChemPtr(i);
                }

                _setInnerList(ChemList, ChemIsKnown, ChemMol);
                _setStreamNum();
            }
            StreamBase(const ChemPtrVec& ChemList, const BoolVec ChemIsKnown, const FloatVec& ChemMol)
            {
                assert(ChemList.size() == ChemIsKnown.size() && ChemList.size() == ChemMol.size());

                _setInnerList(ChemList, ChemIsKnown, ChemMol);
            }
            StreamBase(const ChemPtrVec& ChemList, const std::map<ChemBase*, float>& ChemMolMap)
            // 일부 화학종의 몰 유량을 알 수 없는 경우 사용함.
            {
                BoolVec ChemIsKnown(ChemList.size());
                FloatVec ChemMol(ChemList.size());

                std::map<ChemBase*, float>::const_iterator it;

                for (auto i = 0; i < ChemList.size(); ++i)
                {
                    it = ChemMolMap.find(ChemList[i]);
                    if (it == ChemMolMap.end())
                    {
                        ChemIsKnown[i] = false;
                        ChemMol[i] = 0;
                    }
                    else
                    {
                        ChemIsKnown[i] = true;
                        ChemMol[i] = it->second;
                    }
                }

                _setInnerList(ChemList, ChemIsKnown, ChemMol);
                _setStreamNum();
            }
            StreamBase(const StreamBase& StreamBaseObj)
            // 대입 연산자를 위한 생성자.
            {
                _ChemList = StreamBaseObj._ChemList;
                _ChemIsKnown = StreamBaseObj._ChemIsKnown;
                _ChemMol = StreamBaseObj._ChemMol;

                _setStreamNum();
            }

            // getter/setter 정의부

            int getStreamNum() {return _StreamNum;}
            ChemPtrVec getChemList() {return _ChemList;}
            FloatVec getChemMol() {return _ChemMol;}
            BoolVec getChemIsKnown() {return _ChemIsKnown;}
            
            bool addChem(ChemBase* Chem)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                return _addChem(Chem, false, 0);
            }
            bool addChem(const int& ChemNum)
            // ChemBase::_chemNum을 활용해 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                return _addChem(ChemBase::getChemPtr(ChemNum), false, 0);
            }
            bool addChem(ChemBase* Chem, const float& Mol)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                return _addChem(Chem, true, Mol);
            }
            bool addChem(const int& ChemNum, const float& Mol)
            // ChemBase::_chemNum을 활용해 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                return _addChem(ChemBase::getChemPtr(ChemNum), true, Mol);
            }
            bool addChem(const ChemPtrVec& ChemList)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                bool res = true;
                for (auto Chem : ChemList)
                {
                    if (!_addChem(Chem, false, 0)) res = false;
                }

                return res;
            }
            bool addChem(const IntVec& ChemNumList)
            // ChemBase::_chemNum을 활용해 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                bool res = true;
                for (auto i : ChemNumList)
                {
                    if (!_addChem(ChemBase::getChemPtr(i), false, 0)) res = false;
                }

                return res;
            }
            bool addChem(const ChemPtrVec& ChemList, const FloatVec& ChemMol)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                assert(ChemList.size() == ChemMol.size());

                bool res = true;
                for (auto i = 0; i < ChemList.size(); ++i)
                {
                    if (!_addChem(ChemList[i], true, ChemMol[i])) res = false;
                }

                return res;
            }
            bool addChem(const IntVec& ChemNumList, const FloatVec& ChemMol)
            // ChemBase::_chemNum을 활용해 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                assert(ChemNumList.size() == ChemMol.size());

                bool res = true;
                for (auto i = 0; i < ChemNumList.size(); ++i)
                {
                    if (!_addChem(ChemBase::getChemPtr(ChemNumList[i]), true, ChemMol[i])) res = false;
                }

                return res;
            }
            bool addChem(const ChemPtrVec& ChemList, const std::map<ChemBase*, float>& ChemMolMap)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                bool res = true;
                bool buff;
                std::map<ChemBase*, float>::const_iterator it;
                for (auto Chem : ChemList)
                {
                    buff = true;
                    it = ChemMolMap.find(Chem);
                    if (it == ChemMolMap.end())
                    {
                        buff = _addChem(Chem, false, 0);
                    }
                    else
                    {
                        buff = _addChem(Chem, false, it->second);
                    }

                    if (!buff) res = false;
                }

                return res;
            }
            bool addChem(const ChemPtrVec& ChemList, const BoolVec& ChemIsKnown, const FloatVec& ChemMol)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {   
                assert(ChemList.size() == ChemIsKnown.size() && ChemIsKnown.size() == ChemMol.size());
                bool res = true;
                bool buff;
                for (auto i = 0; i < ChemList.size(); ++i)
                {
                    buff = true;
                    if (ChemIsKnown[i]) buff = _addChem(ChemList[i], true, ChemMol[i]);
                    else buff = _addChem(ChemList[i], false, 0);

                    if (!buff) res = false;
                }

                return res;
            }

            bool setChemMol(ChemBase* Chem, const float& ChemMol)
            // 화학종의 몰 유량을 설정함. 성공하면 true를 반환함.
            {
                return _setChemMol(Chem, ChemMol);
            }
            bool setChemMol(const int& ChemNum, const float& ChemMol)
            // ChemBase::_chemNum을 활용해 화학종의 몰 유량을 설정함. 성공하면 true를 반환함.
            {
                return _setChemMol(ChemBase::getChemPtr(ChemNum), ChemMol);
            }
            bool setChemMol(const ChemPtrVec& ChemList, const FloatVec& ChemMol)
            // 화학종의 몰 유량을 설정함. 성공하면 true를 반환함.
            {
                assert(ChemList.size() == ChemMol.size());

                bool res = true;
                for (auto i = 0; i < ChemList.size(); ++i)
                {
                    if (!setChemMol(ChemList[i], ChemMol[i])) res = false;
                }

                return res;
            }
            bool setChemMol(const IntVec& ChemNumList, const FloatVec& ChemMol)
            // ChemBase::_chemNum을 활용해 화학종의 몰 유량을 설정함. 성공하면 true를 반환함.
            {
                assert(ChemNumList.size() == ChemMol.size());

                bool res = true;
                for (auto i = 0; i < ChemNumList.size(); ++i)
                {
                    if (!setChemMol(ChemNumList[i], ChemMol[i])) res = false;
                }

                return res;
            }

            bool setChemUnknown(ChemBase* Chem)
            // 화학종의 몰 유량을 미지수로 변경함. 성공하면 true를 반환함.
            {
                return _setChemUnknown(Chem);
            }
            bool setChemUnknown(const int& ChemNum)
            // ChemBase::_chemNum을 활용해 화학종의 몰 유량을 미지수로 변경함. 성공하면 true를 반환함.
            {
                return _setChemUnknown(ChemBase::getChemPtr(ChemNum));
            }
            bool setChemUnknown(const ChemPtrVec& ChemList)
            // 화학종의 몰 유량을 미지수로 변경함. 성공하면 true를 반환함.
            {
                bool res = true;
                for (auto Chem : ChemList)
                {
                    if (!_setChemUnknown(Chem)) res = false;
                }

                return res;
            }
            bool setChemUnknown(const IntVec& ChemNumList)
            // ChemBase::_chemNum을 활용해 화학종의 몰 유량을 미지수로 변경함. 성공하면 true를 반환함.
            {
                bool res = true;
                for (auto ChemNum : ChemNumList)
                {
                    if (!_setChemUnknown(ChemBase::getChemPtr(ChemNum))) res = false;
                }

                return res;
            }

            bool delChem(ChemBase* Chem)
            // 스트림으로부터 화학종을 제거함. 성공하면 true를 반환함.
            {
                return _delChem(Chem);
            }
            bool delChem(const int& ChemNum)
            // ChemBase::_chemNum을 활용해 스트림으로부터 화학종을 제거함. 성공하면 true를 반환함.
            {
                return _delChem(ChemBase::getChemPtr(ChemNum));
            }
            bool delChem(const ChemPtrVec& ChemList)
            // 스트림으로부터 화학종을 제거함. 성공하면 true를 반환함.
            {
                bool res = true;
                for (auto Chem : ChemList)
                {
                    if (!_delChem(Chem)) res = false;
                }

                return res;
            }
            bool delChem(const IntVec& ChemNumList)
            // ChemBase::_chemNum을 활용해 스트림으로부터 화학종을 제거함. 성공하면 true를 반환함.
            {
                bool res = true;
                for (auto ChemNum : ChemNumList)
                {
                    if (!_delChem(ChemBase::getChemPtr(ChemNum))) res = false;
                }

                return res;
            }

            // 인스턴스 정의부

            bool inChemList(ChemBase* Chem)
            // _ChemList에 Chem이 있으면 true, 아니면 false를 반환함.
            {
                auto it = std::find(_ChemList.begin(), _ChemList.end(), Chem);
                if (it == _ChemList.end()) return false;
                else return true;
            }

            
    };
    int StreamBase::nextStreamNum = 0;

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
            static int nextProcObjNum;

            // ProcObjBase 객체들의 포인터를 저장함.
            static std::vector<ProcObjBase*> ProcObjList;
            
            StreamPtrVec _inStream;
            StreamPtrVec _outStream;

            inline void _setProcObjNum()
            // _ProcObjNum을 할당함.
            {
                _ProcObjNum = nextProcObjNum++;
                ProcObjList.push_back(this);
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

            std::string __Name;
            Eigen::MatrixXf __MainMat;
            FloatVec __ScalarVec;

        public:

            static auto getProcObjList() {return ProcObjList;}
            static auto getProcObjPtr(const int& i) {return ProcObjList[i];}

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
    int ProcObjBase::nextProcObjNum = 0;

    class RxtorBase : public ProcObjBase
    /*
    화학 반응기를 나타내는 클래스. ProcObjBase로부터 상속됨.
    ProcObjBase에서 화학 반응식 및 몰수지 부분을 추가함.
    이 클래스는 Eigen 3 라이브러리의 Solver 기능을 사용함.
    INTEL(R) Math Kernel Library, OpenBLAS 등이 있는 경우 CMake를 이용해 속도 향상이 가능함.
    */
    {   
        private:

            inline void _setRxtorNum()
            {
                _rxtorNum = nextRxtorNum++;
                std::pair<int, RxtorBase*> toRxtorMap;
                toRxtorMap.first = _rxtorNum;
                toRxtorMap.second = this;
                RxtorBase::RxtorMap.insert(toRxtorMap);
            }

        protected:

            int _rxtorNum;

            std::vector<ChemBase*> _chemList;
            std::map<ChemBase*, int> _chemMap;
            Eigen::MatrixXf _rxnRatio;

        public:

            // 다음 생성자 호출시 부여받는  _rxtorNum.
            static int nextRxtorNum;

            // Rxtorase 객체들의 _rxtorNum과 포인터를 저장함.
            static std::map<int, RxtorBase*> RxtorMap;

            // 생성자 선언부

            RxtorBase():
                ProcObjBase() {_setRxtorNum();};
            
            RxtorBase(const std::string& name):
                ProcObjBase(name) {_setRxtorNum();};
            
            RxtorBase(const std::string& name, StreamBase* inStream, StreamBase* outStream):
                ProcObjBase(name, {inStream}, {outStream}) {_setRxtorNum();};
            
            RxtorBase(const std::string& name, const int& inStreamNum, const int& outStreamNum):
                ProcObjBase(name, {inStreamNum}, {outStreamNum}) {_setRxtorNum();};
            
            RxtorBase(const std::string& name, StreamBase* inStream, StreamBase* outStream,
                const std::vector<float>& scalarVec):
                ProcObjBase(name, {inStream}, {outStream}, scalarVec) {_setRxtorNum();};
            
            RxtorBase(const std::string& name, const int& inStreamNum, const int& outStreamNum,
                const std::vector<float>& scalarVec):
                ProcObjBase(name, {inStreamNum}, {outStreamNum}, scalarVec) {_setRxtorNum();};
            
            //인스턴스 선언부

            bool checkStreamNum()
            // 입력 스트림과 출력 스트림이 모두 1개뿐인지 확인함.
            {
                if (_inStream.size() == 1 && _outStream.size() == 1) return true;
                else return false;
            }

            bool makeChemList()
            // 입력 스트림과 출력 스트림의 물질들을 모두 저장한 _chemList를 생성함. 성공한 경우 true, 실패한 경우 false 반환.
            {
                if (!checkStreamNum()) false;

                _chemList = std::vector<ChemBase*>();
                
                auto inStream = _inStream[0];
                auto outStream = _outStream[0];

                std::map<ChemBase*, bool>::const_iterator it_bool;
                std::map<ChemBase*, float>::const_iterator it_float;

                // to be continued.
            }

            
    };
}

#endif