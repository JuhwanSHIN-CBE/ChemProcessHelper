/*
core/CoreBase/StreamBase.hpp
----------------------------
공정흐름도에서 물질 스트림을 표현하는 StreamBase 클래스를 정의함.
*/
#ifndef _CHEMPROCHELPER_STREAMBASE
#define _CHEMPROCHELPER_STREAMBASE

namespace chemprochelper
{
    class StreamBase
    /*
    화학공정흐름도에서 물질의 흐름(flow stream)을 표현하는 클래스.
    ----------------------------------------------------------
    StreamBase는 다음과 같은 멤버 변수를 가짐.
    private:
        _ChemIdx : 해당 객체에 연관된 ChemBase 객체들의 포인터를 저장함.
        _ChemMask : 해당 화학종의 몰 유량을 알고 있는지의 여부를 저장함.
        _ChemMol : 해당 화학종의 몰 유량을 저장함.
    */
    {
        private:

            // 흐름을 구성하는 화학종의 포인터를 저장함.
            std::vector<ChemBase*> _ChemIdx;

            // 화학종의 몰 유량이 알려져 있으면 true, 아닌 경우 false를 부여함.
            std::vector<bool> _ChemMask;

            // 화학종의 몰 유량을 저장함. 몰 유량을 알 수 없는 경우 0을 저장함.
            std::vector<float> _ChemMol;

            // StreamBase 객체에 화학종을 추가함. 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool _updateChem(ChemBase* ChemIdx, const bool& ChemMask, const float& ChemMol)
            {
                auto it = std::find(_ChemIdx.begin(), _ChemIdx.end(), ChemIdx);
                if (it == _ChemIdx.end())
                {
                    auto idx = it - _ChemIdx.begin();
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
            bool _updateChem(const std::vector<ChemBase*>& ChemIdx,
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
            bool _delChem(ChemBase* ChemIdx)
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
            bool _setChemUnkown(ChemBase* ChemIdx)
            {
                auto it = std::find(_ChemIdx.begin(), _ChemIdx.end(), ChemIdx);
                if (it == _ChemIdx.end()) return false;

                auto idx = it - _ChemIdx.begin();
                _ChemMask[idx] = false;
                _ChemMol[idx] = 0;

                return true;
            }

        public:

            // 생성자 정의부

            // 디폴트 생성자
            StreamBase() = default;

            // 모든 물질의 몰 유량을 모르는 경우
            StreamBase(const std::vector<ChemBase*>& ChemIdx)            
            {
                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    _ChemIdx.push_back(ChemIdx[i]);
                    _ChemMask.push_back(false);
                    _ChemMol.push_back(0);
                }
            }

            // 모든 물질의 몰 유량을 아는 경우
            StreamBase(const std::vector<ChemBase*>& ChemIdx, const std::vector<float>& ChemMol)            
            {
                assert(ChemIdx.size() == ChemMol.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    _ChemIdx.push_back(ChemIdx[i]);
                    _ChemMask.push_back(true);
                    _ChemMol.push_back(ChemMol[i]);
                }
            }

            // 모든 물질의 몰 유량을 아는 경우
            StreamBase(const std::unordered_map<ChemBase*, float>& ChemMol)
            {
                _ChemIdx.resize(ChemMol.size());
                _ChemMask.resize(ChemMol.size());
                _ChemMol.resize(ChemMol.size());

                int idx = 0;
                for (const auto& it : ChemMol)
                {
                    _ChemIdx[idx] = it.first;
                    _ChemMask[idx] = true;
                    _ChemMol[idx] = it.second;

                    ++idx;
                }
            }

            // 일부 물질의 몰 유량만을 아는 경우.
            StreamBase(const std::vector<ChemBase*>& ChemIdx, const std::vector<bool>& ChemMask,
                const std::vector<float>& ChemMol)
            {
                assert(ChemIdx.size() == ChemMask.size() && ChemMask.size() == ChemMol.size());

                _ChemIdx = ChemIdx;
                _ChemMask = ChemMask;
                _ChemMol = ChemMol;
            }

            // 일부 물질의 몰 유량만을 아는 경우.
            StreamBase(const std::vector<ChemBase*>& ChemIdx,
                const std::unordered_map<ChemBase*, float>& ChemMol)
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
            }

            // getter 정의부

            auto getChemIdx() {return _ChemIdx;}
            auto getChemMask() {return _ChemMask;}
            auto getChemMol() {return _ChemMol;}

            // 인스턴스 정의부

            // 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(ChemBase* ChemIdx)            
            {
                return _updateChem(ChemIdx, false, 0);
            }

            // 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<ChemBase*>& ChemIdx)            
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

            // 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(ChemBase* ChemIdx, const float& ChemMol)
            {
                return _updateChem(ChemIdx, true, ChemMol);
            }

            // 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<ChemBase*>& ChemIdx, const std::vector<float>& ChemMol)
            {
                assert(ChemIdx.size() == ChemMol.size());

                std::vector<bool> ChemMask(ChemIdx.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    ChemMask[i] = true;
                }

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<ChemBase*>& ChemIdx, const std::vector<bool>& ChemMask,
                const std::vector<float>& ChemMol)
            {
                assert(ChemIdx.size() == ChemMask.size() && ChemMask.size() == ChemMol.size());

                return _updateChem(ChemIdx, ChemMask, ChemMol);
            }

            // 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::vector<ChemBase*>& ChemIdx, std::unordered_map<ChemBase*, float>& ChemMolMap)
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

            // 기존 화학종의 값을 덮어쓴 경우 false를 반환함.
            bool updateChem(const std::unordered_map<ChemBase*, float>& ChemMolMap)
            {
                std::vector<ChemBase*> ChemIdx(ChemMolMap.size());
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

            // StreamBase 객체에서 화학종을 제거함. 성공한 경우 true를 반환함.
            bool delChem(ChemBase* ChemIdx)
            {
                return _delChem(ChemIdx);
            }

            // StreamBase 객체에서 화학종을 제거함. 성공한 경우 true를 반환함.
            bool delChem(const std::vector<ChemBase*>& ChemIdx)
            {
                bool res = true;

                for (auto idx : ChemIdx)
                {
                    if (!_delChem(idx)) res = false;
                }

                return res;
            }

            // StreamBase 객체의 모든 화학종을 미지수로 변경함.
            void setAllUnknown()
            {
                std::for_each(_ChemMask.begin(), _ChemMask.end(), [](auto& k)->void{k = false;});
                std::for_each(_ChemMol.begin(), _ChemMol.end(), [](auto& k)->void{k = 0;});
            }

            // StreamBase 객체에서 특정 화학종을 미지수로 변경함. 성공한 경우 true를 반환함.
            bool setChemUnknown(ChemBase* ChemIdx)
            {
                return _setChemUnkown(ChemIdx);
            }

            // StreamBase 객체에서 특정 화학종을 미지수로 변경함. 성공한 경우 true를 반환함.
            bool setChemUnknown(const std::vector<ChemBase*>& ChemIdx)
            {
                bool res = true;

                for (auto idx : ChemIdx)
                {
                    if (!_setChemUnkown(idx)) res = false;
                }

                return res;
            }

            // StreamBase 객체에 특정 화학종이 스트림에 포함된 경우 true를 반환한다.
            bool inChemList(ChemBase* ChemIdx)
            {
                return functions::inVector(_ChemIdx, ChemIdx);
            }

            // StreamBase 객체에 특정 화학종이 스트림에 포함된 경우 true를 반환함.
            std::vector<bool> inChemList(const std::vector<ChemBase*>& ChemIdx)
            {
                std::vector<bool> ChemMask(ChemIdx.size());

                for (auto i = 0; i < ChemIdx.size(); ++i)
                {
                    ChemMask[i] = functions::inVector(_ChemIdx, ChemIdx[i]);
                }

                return ChemMask;
            }
    };
} // namespace chemprochelper

#endif