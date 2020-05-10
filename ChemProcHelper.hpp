/*
ChemProcHelper.hpp
----------------------------------------
복잡한 화학 공정에서의 계산을 물질 흐름을 중심으로 빠르고 편리하게 계산함.

본 헤더 파일은 thread-safe 하지 않음.

주요 최상위 클래스 : ChemBase, RxnBase, ProcObjBase
----------------------------------------
ProcObjBase
    <= MixerBase, RxtorBase, SpliterBase
RxnBase
    <= SpeedRxnBase, EnergyRxnBase, StateRxnBase
ChemBase
*/

#ifndef _CHEMPROCHELPER_
#define _CHEMPROCHELPER_

/*
행렬연산이 필요한 경우,
#include _INCLUDE_CHEMPROCHELPER_SOLVER
를 꼭 써줄것!
(컴파일 타임이 길어져서 별도 옵션으로 분리함.)
*/

// 표준 라이브러리
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <string>
#include <regex>
#include <sstream>
#include <algorithm>
#include <set>

/*
이 라이브러리는 Eigen 3 라이브러리를 필수로 요구함.
https://gitlab.com/libeigen/eigen 참조.
INTEL(R) MKL 등이 있는 경우 CMake를 이용할 것.
*/
#include <Eigen/Dense>

/*
이 라이브러리는 boost 라이브러리를 필수로 요구함.
https://www.boost.org/ 참조.
*/
#include <boost/numeric/odeint.hpp>

// 내부 헤더 파일 연결부
#include "core/Internal.hpp"
#include "core/CoreBase.hpp"
#include "core/RxtorFamily.hpp"
#include "core/FlowManagerFamily.hpp"

#endif