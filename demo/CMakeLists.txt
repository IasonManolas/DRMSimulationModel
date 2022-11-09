cmake_minimum_required(VERSION 3.0)
#set(CMAKE_CXX_STANDARD 20)
set(PROJECT_NAME ${PROJECT_NAME}_demo)

file(GLOB DEMO_FILES ${CMAKE_CURRENT_LIST_DIR}/main.cpp)
add_executable(${PROJECT_NAME} ${DEMO_FILES})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

target_link_libraries(${PROJECT_NAME} PRIVATE DRMSimulationModel_lib)
target_include_directories(${PROJECT_NAME} PRIVATE DRMSimulationModel_lib)
target_include_directories(${PROJECT_NAME} PRIVATE ../deps)
