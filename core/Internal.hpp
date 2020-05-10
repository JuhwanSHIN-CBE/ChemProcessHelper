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
        const std::regex pat_bra("([A-Z][a-z]?|\\((?:[^()]*(?:\\(.*\\))?[^()]*)+\\))(\\d*)");
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

        // target이 set에 있으면 true, 아니면 false를 반환한다.
        template<typename T>
        bool inSet(const std::set<T>& set, const T& target)
        {
            auto it = set.find(target);
            if (it == set.end()) return false;
            else return true;
        }

        // target이 unordered_map의 key에 있으면 true, 아니면 false를 반환한다.
        template<typename K, typename V>
        bool inMap(const std::unordered_map<K, V>& map, const K& target)
        {
            auto it = map.find(target);
            if (it == map.end()) return false;
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

        // 화학식으로부터 원자 조성을 추출함.
        std::unordered_map<std::string, int> _getElemComp(const std::string& chem)
        {
            std::unordered_map<std::string, int> elemIdx, buff;

            int effi;
            std::string elem;

            for (auto m : _RegexIter(chem, const_variables::pat_bra))
            {
                if (m[2] == "") effi = 1;
                else effi = std::stoi(m[2]);

                elem = m[1];
                if (elem[0] == '(')
                {
                    buff = _getElemComp(elem);

                    for (auto pair : buff)
                    {
                        pair.second *= effi;

                        if (inMap(elemIdx, pair.first))
                        {
                            elemIdx[pair.first] += pair.second;
                        }
                        else
                        {
                            elemIdx.insert(pair);
                        }
                    }
                }
                else
                {
                    if (inMap(elemIdx, elem)) elemIdx[elem] += effi;
                    else elemIdx[elem] = effi;
                }
            }

            return elemIdx;
        }

        // 균형 잡힌 화학식을 반환함.
        std::string _balRxnEqn(const std::string& eqn)
        {
            auto strIdx = eqn.find("=");
            if (strIdx == -1) throw std::runtime_error("Invalid chemical Reaction has entered.");

            std::string reac = eqn.substr(0, strIdx);
            std::string prod = eqn.substr(strIdx + 1);

            std::vector<std::string> reacVec, prodVec;

            std::vector<std::unordered_map<std::string, int>> reacIdx, prodIdx;

            std::unordered_map<std::string, int> chemIdx, buff;
            
            // 반응물 부분
            for (const auto& m : _RegexIter(reac, const_variables::pat_big))
            {
                reacVec.push_back(m[2]);
                buff = _getElemComp(m[2]);

                printf("Reactant : %s\n", m[2].str().c_str());
                
                for (const auto& pair : buff)
                {
                    if (!inMap(chemIdx, pair.first)) chemIdx[pair.first] = chemIdx.size();
                }

                reacIdx.push_back(buff);
            }

            // 생성물 부분
            for (const auto& m : _RegexIter(prod, const_variables::pat_big))
            {
                prodVec.push_back(m[2]);
                buff = _getElemComp(m[2]);

                printf("Product : %s\n", m[2].str().c_str());

                for (const auto& pair : buff)
                {
                    if (!inMap(chemIdx, pair.first)) chemIdx[pair.first] = chemIdx.size();
                }

                prodIdx.push_back(buff);
            }

            std::printf("The size of reacIdx is %d\n", reacIdx.size());
            std::printf("The size of prodIdx is %d\n", prodIdx.size());

            Eigen::MatrixXf mat(chemIdx.size(), reacVec.size() + prodVec.size());
            mat.setZero();
            Eigen::VectorXf ans(chemIdx.size());
            ans.setZero();

            // 반응물 부분
            for (auto j = 0; j < reacIdx.size(); ++j)
            {
                for (const auto& pair : reacIdx[j])
                {
                    mat(chemIdx[pair.first], j) = pair.second;
                }
            }

            // 생성물 부분
            for (auto j = 0; j < prodIdx.size(); ++j)
            {
                for (const auto& pair : prodIdx[j])
                {
                    mat(chemIdx[pair.first], j + reacVec.size()) = -1 * pair.second;
                }
            }

            ans = -1 * mat.col(0);

            mat.block(0,0,mat.rows(), mat.cols()-1) = mat.block(0,1,mat.rows(),mat.cols()-1);
            mat.conservativeResize(mat.rows(), mat.cols()-1);

            Eigen::VectorXf res = mat.colPivHouseholderQr().solve(ans);

            res.conservativeResize(res.size() + 1);
            for (auto i = res.size() - 1; i > 0; --i) res[i] = res[i-1];
            res[0] = 1;

            std::stringstream ss;
            ss.precision(2);    // 소수점 아래 둘째 자리까지만 남김.

            reac.clear();
            prod.clear();

            for (auto i = 0; i < reacVec.size(); ++i)
            {
                if (abs(res[i] - 1) < 0.01)
                {
                    reac.append(reacVec[i]);
                }
                else
                {
                    ss << res[i];
                    reac.append(ss.str());
                    reac.append(reacVec[i]);
                }
                reac.append(" + ");
            }
            ss.str("");

            for (auto i = 0; i < prodVec.size(); ++i)
            {
                if (abs(res[i+reacVec.size()] - 1) < 0.01)
                {
                    prod.append(prodVec[i]);
                }
                else
                {
                    ss << res[i+reacVec.size()];
                    prod.append(ss.str());
                    prod.append(prodVec[i]);
                }
                prod.append(" + ");
            }
            ss.str("");

            reac.erase(reac.end()-3, reac.end());
            prod.erase(prod.end()-3, prod.end());
            
            reac.append(" = ");
            reac.append(prod);

            return reac;
        }
    }
}