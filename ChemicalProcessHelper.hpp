/*
ChemicalProcessHelper.hpp
복잡한 화학 공정에서의 계산을 편리하게 하기 위함.

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

// 이 코드는 Eigen 라이브러리를 필수적으로 요구함.
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
        protected:

            std::string _name;
            std::string _abb;
            float _mw;

        public:

            // 생성자 정의부

            ChemBase() = default;
            ChemBase(const std::string& abb):
                _name(abb), _abb(abb), _mw(functions::calMw(_abb)) {};
            ChemBase(const std::string& name, const std::string& abb):
                _name(name), _abb(abb), _mw(functions::calMw(_abb)) {};
            ChemBase(const std::string& name, const float& mw):
                _name(name), _abb(""), _mw(mw) {};
            ChemBase(const std::string& name, const std::string& abb, const float& mw):
                _name(name), _abb(abb), _mw(mw) {};
            
            // getter 정의부
            
            auto getName() {return _name;}
            auto getAbb() {return _abb;}
            auto getMw() {return _mw;}

    };

    class StreamBase
    /*
    화학공정흐름도에서 물질의 흐름(flow stream)을 표현하는 클래스.
    */
    {
        protected:

            // _streamNum은 0부터 시작함.
            int _streamNum;

            ProcObjBase* _startObj = nullptr;
            ProcObjBase* _endObj = nullptr;

            std::map<ChemBase*, float> _chemMol;

            // 해당 화학종의 몰수가 알려져 있으면 true, 아니면 false 값을 부여함.
            std::map<ChemBase*, bool> _chemIsKnown;

        public:

            // 다음 생성자 호출시 부여받는 _streamNum.
            static int nextStreamNum;

            // StreamBase 객체들의 _streamNum과 포인터를 저장함.
            static std::map<int, StreamBase*> StreamMap;

            // 생성자 정의부

            StreamBase()
            {
                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);
            }
            StreamBase(const std::vector<ChemBase*>& chemList)
            {
                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);

                std::pair<ChemBase*, bool> toChemIsKnown;
                for (const auto& chem : chemList)
                {
                    toChemIsKnown.first = chem;
                    toChemIsKnown.second = false;
                    _chemIsKnown.insert(toChemIsKnown);
                }
            }
            StreamBase(const std::vector<ChemBase*>& chemList, const std::map<ChemBase*, float>& chemMol)
            // 일부 화학종의 몰유량을 알 수가 없는 경우 사용함.
            {
                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);

                std::pair<ChemBase*, bool> toChemIsKnown;
                std::map<ChemBase*, float>::const_iterator it_chemMol;
                for (const auto& chem : chemList)
                {
                    toChemIsKnown.first = chem;
                    
                    it_chemMol = chemMol.find(chem);
                    if (it_chemMol != chemMol.end())
                    {
                        _chemMol.insert(*it_chemMol);
                        toChemIsKnown.second = true;
                    }
                    else
                    {
                        toChemIsKnown.second = false;
                    }
                    _chemIsKnown.insert(toChemIsKnown);
                }
            }
            StreamBase(const std::vector<ChemBase*>& chemList, const std::map<ChemBase*, float>& chemMol,
                ProcObjBase* startObj, ProcObjBase* endObj)
            // 일부 화학종의 몰유량을 알 수가 없는 경우 사용함.
            {
                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);

                std::pair<ChemBase*, bool> toChemIsKnown;
                std::map<ChemBase*, float>::const_iterator it_chemMol;
                for (const auto& chem : chemList)
                {
                    toChemIsKnown.first = chem;
                    
                    it_chemMol = chemMol.find(chem);
                    if (it_chemMol != chemMol.end())
                    {
                        _chemMol.insert(*it_chemMol);
                        toChemIsKnown.second = true;
                    }
                    else
                    {
                        toChemIsKnown.second = false;
                    }
                    _chemIsKnown.insert(toChemIsKnown);
                }

                _startObj = startObj;
                _endObj = endObj;
            }
            StreamBase(const std::vector<ChemBase*>& chemList, const std::vector<float>& chemMol)
            // 모든 화학종의 몰 유량을 알 수 있는 경우 사용함.
            {
                assert(chemList.size() == chemMol.size());

                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);

                std::pair<ChemBase*, bool> toChemIsKnown;
                std::pair<ChemBase*, float> toChemMol;
                for (auto i = 0; i < chemList.size(); ++i)
                {
                    toChemIsKnown.first = toChemMol.first = chemList[i];
                    toChemIsKnown.second = true;
                    toChemMol.second = chemMol[i];
                    _chemIsKnown.insert(toChemIsKnown);
                    _chemMol.insert(toChemMol);
                }
            }
            StreamBase(const std::vector<ChemBase*>& chemList, const std::vector<float>& chemMol,
                ProcObjBase* startObj, ProcObjBase* endObj)
            // 모든 화학종의 몰 유량을 알 수 있는 경우 사용함.
            {
                assert(chemList.size() == chemMol.size());

                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);

                std::pair<ChemBase*, bool> toChemIsKnown;
                std::pair<ChemBase*, float> toChemMol;
                for (auto i = 0; i < chemList.size(); ++i)
                {
                    toChemIsKnown.first = toChemMol.first = chemList[i];
                    toChemIsKnown.second = true;
                    toChemMol.second = chemMol[i];
                    _chemIsKnown.insert(toChemIsKnown);
                    _chemMol.insert(toChemMol);
                }

                _startObj = startObj;
                _endObj = endObj;
            }
            StreamBase(const std::map<ChemBase*, bool>& chemIsKnown, const std::map<ChemBase*, float>& chemMol)
            // 일부 화학종의 몰유량을 알 수가 없는 경우 사용함.
            {
                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);

                _chemIsKnown = chemIsKnown;
                _chemMol = chemMol;
            }
            StreamBase(const std::map<ChemBase*, bool>& chemIsKnown, const std::map<ChemBase*, float>& chemMol,
                ProcObjBase* startObj, ProcObjBase* endObj)
            // 일부 화학종의 몰유량을 알 수가 없는 경우 사용함.
            {
                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);

                _chemIsKnown = chemIsKnown;
                _chemMol = chemMol;

                _startObj = startObj;
                _endObj = endObj;
            }
            StreamBase(const StreamBase& StreamBaseObj)
            {
                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);

                _startObj = StreamBaseObj._startObj;
                _endObj = StreamBaseObj._endObj;
                _chemMol = StreamBaseObj._chemMol;
                _chemIsKnown = StreamBaseObj._chemIsKnown;
            }

            // getter/setter 정의부

            int getStreamNum() {return _streamNum;}
            ProcObjBase* getStartObj() {return _startObj;}
            ProcObjBase* getEndObj() {return _endObj;}
            std::map<ChemBase*, float> getChemMol() {return _chemMol;}
            std::map<ChemBase*, bool> getChemIsKnown() {return _chemIsKnown;}
            void setStartObj(ProcObjBase* startObj) {_startObj = startObj;}
            void setEndObj(ProcObjBase* endObj) {_endObj = endObj;}
            
            bool addChem(ChemBase* chem)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                std::map<ChemBase*, bool>::const_iterator it_chemIsKnown;
                it_chemIsKnown = _chemIsKnown.find(chem);
                if (it_chemIsKnown != _chemIsKnown.end()) return false;
                
                std::pair<ChemBase*, bool> toChemIsKnown;
                toChemIsKnown.first = chem;
                toChemIsKnown.second = false;
                _chemIsKnown.insert(toChemIsKnown);
                
                return true;
            }
            bool addChem(ChemBase* chem, const float& mol)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                std::map<ChemBase*, bool>::const_iterator it_chemIsKnown;
                it_chemIsKnown = _chemIsKnown.find(chem);
                if (it_chemIsKnown != _chemIsKnown.end()) return false;

                std::pair<ChemBase*, bool> toChemIsKnown;
                std::pair<ChemBase*, float> toChemMol;
                toChemIsKnown.first = toChemMol.first = chem;
                toChemIsKnown.second = false;
                toChemMol.second = mol;
                _chemIsKnown.insert(toChemIsKnown);
                _chemMol.insert(toChemMol);

                return true;
            }
            bool addChem(const std::vector<ChemBase*>& chemList)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                std::map<ChemBase*, bool>::const_iterator it_chemIsKnown;
                for (const auto& chem : chemList)
                {
                    it_chemIsKnown = _chemIsKnown.find(chem);
                    if (it_chemIsKnown == _chemIsKnown.end()) return false;
                }

                std::pair<ChemBase*, bool> toChemIsKnown;
                for (const auto& chem : chemList)
                {
                    toChemIsKnown.first = chem;
                    toChemIsKnown.second = false;
                    _chemIsKnown.insert(toChemIsKnown);
                }

                return true;
            }
            bool addChem(const std::vector<ChemBase*>& chemList, const std::map<ChemBase*, float>& chemMol)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                std::map<ChemBase*, bool>::const_iterator it_chemIsKnown;
                std::map<ChemBase*, float>::const_iterator it_chemMol;
                std::pair<ChemBase*, bool> toChemIsKnown;
                for (const auto& chem : chemList)
                {
                    it_chemIsKnown = _chemIsKnown.find(chem);
                    if (it_chemIsKnown == _chemIsKnown.end()) return false;

                    toChemIsKnown.first = chem;
                    it_chemMol = chemMol.find(chem);
                    if (it_chemMol == chemMol.end())
                    {
                        toChemIsKnown.second = false;
                        _chemIsKnown.insert(toChemIsKnown);
                    }
                    else
                    {
                        toChemIsKnown.second = true;
                        _chemIsKnown.insert(toChemIsKnown);
                        _chemMol.insert(*it_chemMol);
                    }
                }

                return true;
            }
            bool addChem(const std::map<ChemBase*, bool>& chemIsKnown, const std::map<ChemBase*, float>& chemMol)
            // 스트림에 화학종을 추가함. 성공하면 true를 반환함.
            {
                std::map<ChemBase*, bool>::const_iterator it_chemIsKnown;
                std::map<ChemBase*, float>::const_iterator it_chemMol;
                std::pair<ChemBase*, bool> toChemIsKnown;
                for (const auto& pair_chemIsKnown : chemIsKnown)
                {
                    it_chemIsKnown = _chemIsKnown.find(pair_chemIsKnown.first);
                    if (it_chemIsKnown == _chemIsKnown.end()) return false;

                    it_chemMol = chemMol.find(pair_chemIsKnown.first);
                    if (it_chemMol == chemMol.end())
                    {
                        toChemIsKnown.first = pair_chemIsKnown.first;
                        toChemIsKnown.second = false;
                        _chemIsKnown.insert(toChemIsKnown);
                    }
                    else
                    {
                        toChemIsKnown.first = pair_chemIsKnown.first;
                        toChemIsKnown.second = true;
                        _chemIsKnown.insert(toChemIsKnown);
                        _chemMol.insert(*it_chemMol);
                    }
                    
                    return true;
                }
            }

            bool setMol(ChemBase* chem, const float& mol)
            // 화학종의 몰 유량을 설정함.
            {
                std::map<ChemBase*, bool>::const_iterator it_chemIsKnown;
                it_chemIsKnown = _chemIsKnown.find(chem);
                if (it_chemIsKnown == _chemIsKnown.end()) return false;
                else if (it_chemIsKnown->second) _chemMol[chem] = mol;
                else
                {
                    _chemIsKnown[chem] = true;
                    _chemMol[chem] = mol;
                }
                
                return true;
            }
            
            // 인스턴스 정의부

            bool isFullyConnected()
            // _startObj, _endObj가 모두 제대로 정의되었으면 true, 아니면 false를 반환함.
            {
                if (_startObj == nullptr || _endObj == nullptr)
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }

            
    };
    int StreamBase::nextStreamNum = 0;

    template<typename Scalar>
    class ProcObjBase
    /*
    MixerBase, RxtorBase, SpliterBase의 상위 클래스.
    */
    {
        public:

            typedef typename Eigen::Matrix<Eigen::Dynamic, Eigen::Dynamic, Scalar> DynMat;

        protected:

            std::vector<StreamBase*> _inStream;
            std::vector<StreamBase*> _outStream;
            DynMat _mainMat;
            std::vector<float> _ratioVec;
    };

}

#endif