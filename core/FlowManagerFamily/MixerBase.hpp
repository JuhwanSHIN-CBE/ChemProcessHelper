/*
core/FlowManagerFamily/MixerBase.hpp
------------------------------------
혼합기의 기초가 되는 MixerBase 클래스를 정의함.
*/
#ifndef _CHEMPROCHELPER_MIXERBASE
#define _CHEMPROCHELPER_MIXERBASE

namespace chemprochelper
{
    /*
    분리기를 지정하는 기본 클래스
    */
    class MixerBase : public ProcObjBase
    {
        private:

            // 출력 스트림에 포함된 모든 화학종이 입력 스트림의 모든 화학종과 동일한지 확인함.
            bool _checkStreamValid()
            {
                std::set<ChemBase*> inChemSet;
                std::set<ChemBase*> outChemSet;

                for (auto inStreamPtr : getInStreamIdx())
                {
                    for (auto inChemIdx : inStreamPtr->getChemIdx())
                    {
                        if (!functions::inSet(inChemSet, inChemIdx)) inChemSet.insert(inChemIdx);
                    }
                }

                for (auto outStreamPtr : getOutStreamIdx())
                {
                    for (auto outChemIdx : outStreamPtr->getChemIdx())
                    {
                        if (!functions::inSet(outChemSet, outChemIdx)) outChemSet.insert(outChemIdx);
                    }
                }

                return (inChemSet == outChemSet);
            }
        
        public:

            // 생성자 정의부

            // 임시 객체를 위한 생성자.
            MixerBase():
                ProcObjBase() {}
            
            // StreamBase* 포인터를 이용함. 입/출력 스트림이 정의된 경우
            MixerBase(const std::vector<StreamBase*>& inStreamPtr, StreamBase* outStreamPtr):
                ProcObjBase(inStreamPtr, std::vector<StreamBase*>(1, outStreamPtr)) {}

            // StreamBase* 포인터를 이용함. 입/출력 스트림이 정의되고, 코멘트도 입력하는 경우.
            MixerBase(const std::vector<StreamBase*>& inStreamPtr, StreamBase* outStreamPtr,
                const std::string& Comment):
                ProcObjBase(inStreamPtr, std::vector<StreamBase*>(1, outStreamPtr), Comment) {}
    };
} // namespace chemprochelper

#endif