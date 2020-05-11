/*
core/CSTR.hpp
-----------------
CSTR에 대한 클래스를 정의함.
*/
#ifndef _CHEMPROCHELPER_CSTR
#define _CHEMPROCHELPER_CSTR

namespace chemprochelper
{
    /*
    CSTR을 지정하는 기본 클래스.
    */
    class CSTR : public RxtorBase
    {
        private:

            // _ChemMol

            // __ChemIdx로부터 _ChemMold의 값을 0으로 초기화함.
            void _resetChemMol()
            {
                __ChemMol = std::vector<float>(__ChemIdx.size(), 0);
            }

        public:

            // 생성자 정의부

            // 임시 객체를 위한 생성자
            CSTR():
                RxtorBase() {}
            
            // Comment를 사용하지 않는 경우
            CSTR(StreamBase* inStreamPtr, StreamBase* outStreamPtr, RxnBase* RxnPtr):
                RxtorBase(inStreamPtr, outStreamPtr, RxnPtr)
            {
                _resetChemMol();
            }

            // Comment를 사용하는 경우
            CSTR(StreamBase* inStreamPtr, StreamBase* outStreamPtr, RxnBase* RxnPtr,
                std::string& Comment):
                RxtorBase(inStreamPtr, outStreamPtr, RxnPtr, Comment)
            {
                _resetChemMol();
            }

            
    };
} // namespace chemprochelper


#endif