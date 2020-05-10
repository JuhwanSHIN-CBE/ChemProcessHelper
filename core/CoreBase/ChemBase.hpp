/*
core/CoreBase/ChemBase.hpp
--------------------------
화합물의 기초가 되는 ChemBase 클래스를 정의함.
*/
#ifndef _CHEMPROCHELPER_CHEMBASE
#define _CHEMPROCHELPER_CHEMBASE

namespace chemprochelper
{
    /*
    화합물을 지정하는 기본 클래스
    ---------------------------
    ChemBase는 다음과 같은 멤버 변수를 가짐.
    private:
        _Name : 화학종의 이름을 저장함.
        _Abb : 화학종의 축약형을 저장함.
    */
    class ChemBase
    {
        private:
            
            // 화합물의 명칭을 저장함.
            std::string _Name;

            // 화합물의 축약형을 저장함.
            std::string _Abb;

            // 화합물의 축약형과 그에 해당하는 ChemBase 객체의 포인터를 저장함.
            static std::unordered_map<std::string, ChemBase*> _AbbMap;

        public:

            // 생성자 정의부

            // 디폴트 생성자
            ChemBase() = default;

            // 축약형만 입력된 경우 _Name과 _Abb에 축약형을 저장함.
            ChemBase(const std::string& Abb):
                _Name(Abb), _Abb(Abb)
            {
                auto it = _AbbMap.find(Abb);
                if (it != _AbbMap.end()) throw std::runtime_error("Abb "+Abb+" is already in ChemBase::_ChemList");
                _AbbMap[Abb] = this;
            }

            ChemBase(const std::string& Name, const std::string& Abb):
                    _Name(Name), _Abb(Abb)
            {
                auto it = _AbbMap.find(Abb);
                if (it != _AbbMap.end()) throw std::runtime_error("Abb "+Abb+" is already in ChemBase::_ChemList");
                _AbbMap[Abb] = this;
            }

            ~ChemBase()
            {
                auto it = _AbbMap.find(_Abb);
                if (it->second == this)
                {
                    _AbbMap.erase(it);
                }
            }

            // 대입 연산자 정의부

            // other이 우측값인 경우 대입 직후 other이 소멸하므로, ChemBase::_AbbMap에 등록된 포인터 주소를 변경함.
            ChemBase& operator=(ChemBase&& other)
            {
                _Name = other._Name;
                _Abb = other._Abb;

                _AbbMap[_Abb] = this;

                return *this;
            }

            // other이 좌측값인 경우 대입 이후에도 other이 잔존하므로, ChemBase::_AbbMap을 변경하지 않음.
            ChemBase& operator=(ChemBase& other)
            {
                _Name = other._Name;
                _Abb = other._Name;

                return *this;
            }

            // getter 정의부

            auto getName() {return _Name;}
            auto getAbb() {return _Abb;}

            // 인스턴스 정의부

            // 축약형에 맞는 포인터를 반환함.
            static auto getChemPtr(const std::string& Abb)
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

            // iostream 지원을 위한 함수
            std::string getString() const
            {
                std::string retVal = "Name : " + _Name + ", Abb : " + _Abb;
                return retVal;
            }
    };
    std::unordered_map<std::string, ChemBase*> ChemBase::_AbbMap;

    std::ostream& operator<<(std::ostream& os, const ChemBase& obj)
    {
        os << obj.getString();
        return os;
    }
} // namespace chemprochelper

#endif