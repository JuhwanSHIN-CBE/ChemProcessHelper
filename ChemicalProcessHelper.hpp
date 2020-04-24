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

            inline void _setChemNum()
            // _chemNum을 할당함.
            {
                _chemNum = nextChemNum++;
                std::pair<int, ChemBase*> toChemMap;
                toChemMap.first = _chemNum;
                toChemMap.second = this;
                ChemBase::ChemMap.insert(toChemMap);
            };


        protected:

            int _chemNum;
            std::string _name;
            std::string _abb;
            float _mw;

        public:

            // 다음 생성자 호출시 부여받는 _chemNum.
            static int nextChemNum;

            // ChemBase 객체들의 _chemNum과 포인터를 저장함.
            static std::map<int, ChemBase*> ChemMap;

            // 생성자 정의부

            ChemBase() {_setChemNum();};
            
            ChemBase(const std::string& abb):
                _name(abb), _abb(abb), _mw(functions::calMw(_abb)) {_setChemNum();};
            ChemBase(const std::string& name, const std::string& abb):
                _name(name), _abb(abb), _mw(functions::calMw(_abb)) {_setChemNum();};
            ChemBase(const std::string& name, const float& mw):
                _name(name), _abb(""), _mw(mw) {_setChemNum();};
            ChemBase(const std::string& name, const std::string& abb, const float& mw):
                _name(name), _abb(abb), _mw(mw) {_setChemNum();};

            // getter 정의부
            
            auto getName() {return _name;}
            auto getAbb() {return _abb;}
            auto getMw() {return _mw;}
    };
    int ChemBase::nextChemNum = 0;

    class StreamBase
    /*
    화학공정흐름도에서 물질의 흐름(flow stream)을 표현하는 클래스.
    */
    {
        private:

            inline void _setStreamNum()
            // _streamNum을 설정함.
            {
                _streamNum = StreamBase::nextStreamNum++;
                std::pair<int, StreamBase*> toStreamMap;
                toStreamMap.first = _streamNum;
                toStreamMap.second = this;
                StreamBase::StreamMap.insert(toStreamMap);
            }

        protected:

            // _streamNum은 0부터 시작함.
            int _streamNum;

            std::map<ChemBase*, float> _chemMol;

            // 해당 화학종의 몰수가 알려져 있으면 true, 아니면 false 값을 부여함.
            std::map<ChemBase*, bool> _chemIsKnown;

        public:

            // 다음 생성자 호출시 부여받는 _streamNum.
            static int nextStreamNum;

            // StreamBase 객체들의 _streamNum과 포인터를 저장함.
            static std::map<int, StreamBase*> StreamMap;

            // 생성자 정의부

            StreamBase() {_setStreamNum();};
            
            StreamBase(const std::vector<ChemBase*>& chemList)
            {
                std::pair<ChemBase*, bool> toChemIsKnown;
                for (const auto& chem : chemList)
                {
                    toChemIsKnown.first = chem;
                    toChemIsKnown.second = false;
                    _chemIsKnown.insert(toChemIsKnown);
                }

                _setStreamNum();
            }
            StreamBase(const std::vector<ChemBase*>& chemList, const std::map<ChemBase*, float>& chemMol)
            // 일부 화학종의 몰유량을 알 수가 없는 경우 사용함.
            {
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

                _setStreamNum();
            }
            StreamBase(const std::vector<int>& chemNumList, const std::map<ChemBase*, float>& chemMol)
            // ChemBase::_chemNum을 활용한 생성자. 일부 화학종의 몰유량을 알 수가 없는 경우 사용함.
            {
                std::vector<ChemBase*> chemList;
                std::map<int, ChemBase*>::const_iterator it;
                for (auto i : chemNumList)
                {
                    it = ChemBase::ChemMap.find(i);
                    if (it == ChemBase::ChemMap.end()) throw std::runtime_error("wrong ChemBase::_chemNum");
                    chemList.push_back(it->second);
                }

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

                _setStreamNum();
            }
            StreamBase(const std::vector<ChemBase*>& chemList, const std::vector<float>& chemMol)
            // 모든 화학종의 몰 유량을 알 수 있는 경우 사용함.
            {
                assert(chemList.size() == chemMol.size());
                
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

                _setStreamNum();
            }
            StreamBase(const std::vector<int>& chemNumList, const std::vector<float>& chemMol)
            // ChemBase::_chemNum을 활용한 생성자. 모든 화학종의 몰 유량을 알 수 있는 경우 사용함.
            {
                assert(chemNumList.size() == chemMol.size());

                std::vector<ChemBase*> chemList;
                std::map<int, ChemBase*>::const_iterator it;
                for (auto i : chemNumList)
                {
                    it = ChemBase::ChemMap.find(i);
                    if (it == ChemBase::ChemMap.end()) throw std::runtime_error("wrong ChemBase::_chemNum");
                    chemList.push_back(it->second);
                }

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

                _setStreamNum();
            }
            StreamBase(const std::map<ChemBase*, bool>& chemIsKnown, const std::map<ChemBase*, float>& chemMol)
            // 일부 화학종의 몰유량을 알 수가 없는 경우 사용함.
            {
                _chemIsKnown = chemIsKnown;
                _chemMol = chemMol;

                _setStreamNum();
            }
            StreamBase(const StreamBase& StreamBaseObj)
            {
                _chemMol = StreamBaseObj._chemMol;
                _chemIsKnown = StreamBaseObj._chemIsKnown;

                _setStreamNum();
            }

            // getter/setter 정의부

            int getStreamNum() {return _streamNum;}
            std::map<ChemBase*, float> getChemMol() {return _chemMol;}
            std::map<ChemBase*, bool> getChemIsKnown() {return _chemIsKnown;}
            
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

            
    };
    int StreamBase::nextStreamNum = 0;

    class ProcObjBase
    /*
    MixerBase, RxtorBase, SpliterBase의 상위 클래스.
    */
    {
        private:

            inline void _setObjNum()
            // _objNum을 할당함.
            {
                _objNum = nextObjNum++;
                std::pair<int, ProcObjBase*> toObjMap;
                toObjMap.first = _objNum;
                toObjMap.second = this;
                ProcObjBase::ObjMap.insert(toObjMap);
            }

        protected:

            int _objNum;
            std::string _name;
            std::vector<StreamBase*> _inStream;
            std::vector<StreamBase*> _outStream;
            Eigen::MatrixXf _mainMat;
            std::vector<float> _scalarVec;

        public:

            // 다음 생성자 호출시 부여받는 _objNum.
            static int nextObjNum;

            // StreamBase 객체들의 _objNum과 포인터를 저장함.
            static std::map<int, ProcObjBase*> ObjMap;

            // 생성자 선언부

            ProcObjBase() {_setObjNum();};

            ProcObjBase(const std::string& name):
                _name(name) {_setObjNum();};

            ProcObjBase(const std::string& name, const Eigen::MatrixXf& mainMat):
                _name(name), _mainMat(mainMat) {_setObjNum();};

            ProcObjBase(const std::string& name, const std::vector<StreamBase*>& inStream,
                const std::vector<StreamBase*>& outStream):
                _name(name), _inStream(inStream), _outStream(outStream) {_setObjNum();};

            ProcObjBase(const std::string& name, const std::vector<int>& inStreamNum,
                const std::vector<int>& outStreamNum):
                _name(name)
            // StreamNum을 활용한 생성자 호출.
            {
                std::map<int, StreamBase*>::const_iterator it;
                for (auto i : inStreamNum)
                {
                    it = StreamBase::StreamMap.find(i);
                    if (it == StreamBase::StreamMap.end()) throw std::runtime_error("wrong Stream::_streamNum");
                    _inStream.push_back(it->second);
                }
                for (auto o : outStreamNum)
                {
                    it = StreamBase::StreamMap.find(o);
                    if (it == StreamBase::StreamMap.end()) throw std::runtime_error("wrong Stream::_streamNum");
                    _outStream.push_back(it->second);
                }

                _setObjNum();
            };

            ProcObjBase(const std::string& name, const std::vector<StreamBase*>& inStream,
                const std::vector<StreamBase*>& outStream, const std::vector<float>& scalarVec):
                _name(name), _inStream(inStream), _outStream(outStream), _scalarVec(scalarVec) {_setObjNum();};

            ProcObjBase(const std::string& name, const std::vector<int>& inStreamNum,
                const std::vector<int>& outStreamNum, const std::vector<float>& scalarVec):
                _name(name), _scalarVec(scalarVec)
            {
                std::map<int, StreamBase*>::const_iterator it;
                for (auto i : inStreamNum)
                {
                    it = StreamBase::StreamMap.find(i);
                    if (it == StreamBase::StreamMap.end()) throw std::runtime_error("wrong Stream::_streamNum");
                    _inStream.push_back(it->second);
                }
                for (auto o : outStreamNum)
                {
                    it = StreamBase::StreamMap.find(o);
                    if (it == StreamBase::StreamMap.end()) throw std::runtime_error("wrong Stream::_streamNum");
                    _outStream.push_back(it->second);
                }

                _setObjNum();
            }

            // setter/getter 정의부

            std::string getName() {return _name;}
            std::vector<StreamBase*> getInStream() {return _inStream;}
            std::vector<StreamBase*> getOutStream() {return _outStream;}
            Eigen::MatrixXf getMainMat() {return _mainMat;}
            std::vector<float> getScalarVec() {return _scalarVec;}

            void setName(const std::string& name) {_name = name;}
            void setMainMat(const Eigen::MatrixXf& mainMat) {_mainMat = mainMat;}
            void setScalarVec(const std::vector<float>& ScalarVec) {_scalarVec = ScalarVec;}
            
            bool addStream(StreamBase* stream, bool direction)
            // direction = true이면 in, false이면 out. 추가에 성공하면 true, 아닌 경우 false 반환.
            {
                std::vector<StreamBase*>::const_iterator it;
                if (direction)
                {
                    it = std::find(_inStream.begin(), _inStream.end(), stream);
                    if (it != _inStream.end()) return false;
                    _inStream.push_back(stream);
                }
                else
                {
                    it = std::find(_outStream.begin(), _outStream.end(), stream);
                    if (it != _inStream.end()) return false;
                    _outStream.push_back(stream);
                }

                return true;
            }
            bool addStream(int streamNum, bool direction)
            // StreamBase::_streamNum을 사용함. direction = true이면 in, false이면 out. 추가에 성공하면 true, 아닌 경우 false 반환.
            {
                std::map<int, StreamBase*>::const_iterator c_it = StreamBase::StreamMap.find(streamNum);
                if (c_it == StreamBase::StreamMap.end()) throw std::runtime_error("wrong streamNum");
                auto stream = c_it->second;

                std::vector<StreamBase*>::const_iterator it;
                if (direction)
                {
                    it = std::find(_inStream.begin(), _inStream.end(), stream);
                    if (it != _inStream.end()) return false;
                    _inStream.push_back(stream);
                }
                else
                {
                    it = std::find(_outStream.begin(), _outStream.end(), stream);
                    if (it != _inStream.end()) return false;
                    _outStream.push_back(stream);
                }

                return true;                
            }

            bool delStream(StreamBase* stream, bool direction)
            // direction = true이면 in, false이면 out. 제거에 성공하면 ture, 아닌 경우 false 반환.
            {   
                std::vector<StreamBase*>::iterator it;
                if (direction)
                {
                    it = std::find(_inStream.begin(), _inStream.end(), stream);
                    if (it == _inStream.end()) return false;
                    _inStream.erase(it);
                }
                else
                {
                    it = std::find(_outStream.begin(), _outStream.end(), stream);
                    if (it == _outStream.end()) return false;
                    _outStream.erase(it);
                }
                
                return true;
            }
            bool delStream(int streamNum, bool direction)
            // StreamBase::_streamNum을 사용함. direction = true이면 in, false이면 out. 제거에 성공하면 ture, 아닌 경우 false 반환.
            {
                std::map<int, StreamBase*>::const_iterator c_it = StreamBase::StreamMap.find(streamNum);
                if (c_it == StreamBase::StreamMap.end()) throw std::runtime_error("wrong streamNum");
                auto stream = c_it->second;

                std::vector<StreamBase*>::iterator it;
                if (direction)
                {
                    it = std::find(_inStream.begin(), _inStream.end(), stream);
                    if (it == _inStream.end()) return false;
                    _inStream.erase(it);
                }
                else
                {
                    it = std::find(_outStream.begin(), _outStream.end(), stream);
                    if (it == _outStream.end()) return false;
                    _outStream.erase(it);
                }
                
                return true;
            }
    };
    int ProcObjBase::nextObjNum = 0;

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