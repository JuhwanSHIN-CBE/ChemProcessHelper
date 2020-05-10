namespace chemprochelper
{
    /*
    MixerBase, RxtorBase, SpliterBase의 상위 클래스.
    */
    class ProcObjBase
    {
        private:

            // 반응기 등을 구성하는 입력 스트림의 포인터를 저장함.
            std::vector<StreamBase*> _inStreamIdx;

            // 반응기 등을 구성하는 출력 스트림의 포인터를 저장함.
            std::vector<StreamBase*> _outStreamIdx;

            // 반응기 등에 대한 간단한 메모를 할 수 있음.
            std::string _Comment;

        protected:

            // 반응기 등을 구성하는 화합물들의 포인터를 저장함.
            std::vector<ChemBase*> __ChemIdx;
            
            // 반응기 등에서 중요한 정보들(전화율 등)을 저장함. 자녀 클래스마다 저장하는 값이 다름.
            std::vector<float> __ScalarVec;

            // 반응기 등에서 중요한 정보들(행렬 등)을 저장함. 자녀 클래스마다 저장하는 값이 다름.
            Eigen::MatrixXf __MainMat;

        public:

            // 생성자 정의부

            // 디폴트 생성자
            ProcObjBase() = default;

            // 입/출력 스트림이 정의된 경우
            ProcObjBase(const std::vector<StreamBase*>& inStreamIdx,
                const std::vector<StreamBase*>& outStreamIdx):
                _inStreamIdx(inStreamIdx), _outStreamIdx(outStreamIdx) {}

            // 입/출력 스트림이 정의되고 코멘트 또한 남기는 경우
            ProcObjBase(const std::vector<StreamBase*>& inStreamIdx,
                const std::vector<StreamBase*>& outStreamIdx, const std::string& Comment):
                _inStreamIdx(inStreamIdx), _outStreamIdx(outStreamIdx), _Comment(Comment) {};
            

            // getter 정의부

            auto getInStreamIdx() {return _inStreamIdx;}
            auto getOutStreamIdx() {return _outStreamIdx;}
            auto getComment() {return _Comment;}
            auto getChemIdx() {return __ChemIdx;}
            auto getScalarVec() {return __ScalarVec;}
            auto getMainMat() {return __MainMat;}
    };
} // namespace chemprochelper
