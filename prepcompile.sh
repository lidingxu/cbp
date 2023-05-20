# only some specific gcc compilor version (9.4) can speed up the for loop in oracleSepValue
PREFIX=/home/lxu/software/gcc94
export C_INCLUDE_PATH=$PREFIX/include
export CPLUS_INCLUDE_PATH=$PREFIX/include
export LD_LIBRARY_PATH=$PREFIX/lib:$PREFIX/lib64:$LD_LIBRARY_PATH
export PATH=$PREFIX/bin:$PATH
export CXX=$PREFIX/bin/g++-9.4
