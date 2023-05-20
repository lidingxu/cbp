# mkdir build
# cd build
# cmake ..
# make
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CMAKE_MAKE_PROGRAM"${CMAKE_MAKE_PROGRAM} -j ${N}")      
endif()
cmake_minimum_required(VERSION 3.3)
set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3 -march=native -flto -pthread")
project(Binpacking)


include_directories(${SCIP_DIR}/src  ${SCIP_DIR}/build ${CPLEX_DIR}/concert/include ${CPLEX_DIR}/cplex/include) 
link_directories(${SCIP_DIR}/build/lib ${CPLEX_DIR}/cplex/lib/x86-64_linux/static_pic ${CPLEX_DIR}/concert/lib/x86-64_linux/static_pic)



add_executable(cbp
  srccpp/branch_ryanfoster.cpp
  srccpp/cppmain.cpp
  srccpp/cons_samediff.cpp
  srccpp/pricer_cknap.cpp
  srccpp/probdata_cbp.cpp
  srccpp/reader_cbp.cpp
  srccpp/reader_bp.cpp
  srccpp/conic_knapsack.cpp
  srccpp/conic_knapsack_heur.cpp
  srccpp/conic_knapsack_pwlbc.cpp
  srccpp/estimator.cpp
  srccpp/statistics_table.cpp 
  srccpp/rmp_heur.cpp
)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()


target_link_libraries(cbp -lm -lpthread -ldl -lilocplex -lconcert -lcplex -ldl -lscip ${LIBM})

