cmake_minimum_required(VERSION 3.14)
project(DRMSimulationModel_lib)
#set(LIB_NAME DRMSimulationModel_lib)

file(GLOB LIB_SOURCES *.cpp)
file(GLOB LIB_HEADERS ../include/*.hpp)
add_library(${PROJECT_NAME} STATIC ${LIB_SOURCES} ${LIB_HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC ../include/)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
else()
        add_compile_definitions(ENABLE_PARALLEL_DRM)
endif()

##nlhoman json
target_include_directories(${PROJECT_NAME} PUBLIC ../deps)
##polyscope
if(${USE_POLYSCOPE})
        target_link_libraries(${PROJECT_NAME} PUBLIC polyscope)
        target_include_directories(${PROJECT_NAME} PUBLIC ${polyscope_SOURCE_DIR}/include)
        target_sources(${PROJECT_NAME} PUBLIC ${polyscope_SOURCE_DIR}/deps/imgui/imgui/misc/cpp/imgui_stdlib.h ${polyscope_SOURCE_DIR}/deps/imgui/imgui/misc/cpp/imgui_stdlib.cpp)
endif()
##vcglib
target_include_directories(${PROJECT_NAME} PUBLIC ${vcglib_SOURCE_DIR})
target_sources(${PROJECT_NAME} PUBLIC ${vcglib_SOURCE_DIR}/wrap/ply/plylib.cpp)
##matplot++
target_include_directories(${PROJECT_NAME} PRIVATE ${matplot_SOURCE_DIR}/source)
target_link_libraries(${PROJECT_NAME} PRIVATE matplot)
##threedbeamfea
target_link_libraries(${PROJECT_NAME} PUBLIC ThreedBeamFEA)
target_include_directories(${PROJECT_NAME} PUBLIC ${threed-beam-fea_SOURCE_DIR})
##tbb
target_link_libraries(${PROJECT_NAME} PRIVATE tbb_static)
##eigen3
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
