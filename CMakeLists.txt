cmake_minimum_required(VERSION 3.14)
project(DRMSimulationModel
        VERSION 1.0
        LANGUAGES CXX)

set(USE_POLYSCOPE false)
#Download external dependencies
include(FetchContent)
if(${USE_POLYSCOPE})
        add_compile_definitions(POLYSCOPE_DEFINED)
        ##polyscope
        FetchContent_Declare(polyscope
                                 GIT_REPOSITORY     https://github.com/nmwsharp/polyscope.git
                                 GIT_TAG	    master
          )
        FetchContent_MakeAvailable(polyscope)
endif()

##vcglib devel branch
FetchContent_Declare(vcglib
                 GIT_REPOSITORY     https://gitea-s2i2s.isti.cnr.it/manolas/vcglib.git
                 GIT_TAG	    devel
  )
FetchContent_MakeAvailable(vcglib)

##matplot++ lib
FetchContent_Declare(matplot
                 GIT_REPOSITORY     https://github.com/alandefreitas/matplotplusplus
                 GIT_TAG		master
  )
FetchContent_MakeAvailable(matplot)

##threed-beam-fea
FetchContent_Declare(threed-beam-fea
                 GIT_REPOSITORY      https://gitea-s2i2s.isti.cnr.it/manolas/threed-beam-fea.git
                 GIT_TAG		master
  )
FetchContent_MakeAvailable(threed-beam-fea)

set(TBB_BUILD_SHARED          OFF)
set(TBB_BUILD_STATIC          ON)
set(TBB_BUILD_TBBMALLOC       OFF)
set(TBB_BUILD_TBBMALLOC_PROXY OFF)
set(TBB_BUILD_TESTS           OFF)
set(TBB_NO_DATE               OFF)
set(TBB_BUILD_PYTHON          OFF)
set(TBB_SET_SOVERSION         OFF)
set(TBB_BUILD_TESTS           OFF)
FetchContent_Declare(tbb
                         GIT_REPOSITORY      https://github.com/wjakob/tbb.git
                         GIT_TAG             master
  )
FetchContent_MakeAvailable(tbb)
FetchContent_Declare(eigen
                         GIT_REPOSITORY      https://gitlab.com/libeigen/eigen.git
                         GIT_TAG             master
  )
FetchContent_MakeAvailable(eigen)

##Eigen3
find_package(Eigen3 3.4 REQUIRED)
if(MSVC)
        add_compile_definitions(_HAS_STD_BYTE=0)
endif(MSVC)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

if(${PROJECT_IS_TOP_LEVEL})
        set(SHOULD_BUILD_DEMO true)
else()
        set(SHOULD_BUILD_DEMO false)
endif()

if(${SHOULD_BUILD_DEMO})
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/demo)
endif()
