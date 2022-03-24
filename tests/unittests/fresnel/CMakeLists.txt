project(TestFresnel)
cmake_minimum_required(VERSION 3.0)
set (CMAKE_C_STANDARD 99)
set (ADDA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/)
set (SOURCE_FILES TestFresnel.c ${ADDA_SOURCE_DIR}/cmplx.c)
find_library(CHECK_LIBRARY REQUIRED
        NAMES Check check)
set(CHECK_LIBRARIES "${CHECK_LIBRARY}")
include_directories(${CHECK_INCLUDE_DIRS})
MARK_AS_ADVANCED( CHECK_INCLUDE_DIR CHECK_LIBRARIES )

add_executable(test_fresnel ${SOURCE_FILES})
target_link_libraries(test_fresnel ${CHECK_LIBRARIES})