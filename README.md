# EMP-agmpc 
![arm](https://github.com/emp-toolkit/emp-agmpc/workflows/arm/badge.svg)
![x86](https://github.com/emp-toolkit/emp-agmpc/workflows/x86/badge.svg)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/emp-toolkit/emp-agmpc.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/emp-toolkit/emp-agmpc/alerts/)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/emp-toolkit/emp-agmpc.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/emp-toolkit/emp-agmpc/context:cpp)

<img src="https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/art/logo-full.jpg" width=300px/>

## Global-Scale Secure Multiparty Computation

More details of the protocol can be found in the [paper](https://eprint.iacr.org/2017/189).

<img src="https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/art/logo-full.jpg" width=300px/>


# Installation
1. `wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py`
2. `python[3] install.py --deps --tool --ot --agmpc`
    1. By default it will build for Release. `-DCMAKE_BUILD_TYPE=[Release|Debug]` option is also available.
    2. No sudo? Change [`CMAKE_INSTALL_PREFIX`](https://cmake.org/cmake/help/v2.8.8/cmake.html#variable%3aCMAKE_INSTALL_PREFIX).


### Question
Please send email to wangxiao1254@gmail.com

