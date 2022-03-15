include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
  set(CMAKE_MAKE_PROGRAM"${CMAKE_MAKE_PROGRAM} -j ${N}")      
endif()
cmake_minimum_required(VERSION 3.3)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_FLAGS "-g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g ")
project(cbp)


find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})
include_directories("/opt/ibm/ILOG/CPLEX_Studio201/cplex/include" "/opt/ibm/ILOG/CPLEX_Studio201/concert/include")
link_directories("/opt/ibm/ILOG/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic" "/opt/ibm/ILOG/CPLEX_Studio201/concert/lib/x86-64_linux/static_pic")


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
   srccpp/conic_knapsack_approx.cpp
   srccpp/estimator.cpp
   srccpp/statistics_table.cpp
   srccpp/rmp_heur.cpp
)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()


target_link_libraries(cbp -lilocplex -lconcert -lcplex -lm -lpthread -ldl ${SCIP_LIBRARIES} ${LIBM})

