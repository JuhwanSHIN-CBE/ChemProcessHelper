/*
core/RxnFamily/SpeedRxnBase.hpp
-------------------------------
화학 반응 속도식을 포함하는 화학 반응식 클래스인 SpeedRxnBase를 정의함.
*/
#ifndef _CHEMPROCHELPER_SPEEDRXNBASE
#define _CHEMPROCHELPER_SPEEDRXNBASE

namespace chemprochelper
{
    /*
    속도식을 포함한 화학 반응식을 지정하는 기본 클래스.
    -----------------------------------------------
    SpeedRxnBase는 다음과 같은 멤버변수를 가짐.
    private:
        _SpeedFunc : 화학 반응 속도식(std::functional)을 저장함.
    */
    class SpeedRxnBase : public RxnBase
    {
        public:
            
            // 화학 반응 속도식을 저장함.
            
    }
} // namespace chemprochelper


#endif