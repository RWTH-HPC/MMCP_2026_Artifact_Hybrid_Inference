cmake_minimum_required(VERSION 2.8.1)

include_directories(${INCLUDE_DIR_ABS})
include_directories(SYSTEM ${INCLUDE_DIRS})
include_directories(${SRC_DIR_ABS})
include_directories(${CMAKE_SOURCE_DIR}/${DOCTEST_DIR})

# Select all .cpp files located in the test dir to be built
file(GLOB tst_lib CONFIGURE_DEPENDS
  ${TST_DIR_ABS}/*.cpp
  )

add_executable(${TST_TARGET} $<TARGET_OBJECTS:${OBJ_LIB_TARGET}> ${tst_lib})
set_target_properties(${TST_TARGET} PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(${TST_TARGET} PROPERTIES OUTPUT_NAME ${TST_TARGET})
target_link_libraries(${TST_TARGET} ${LIBRARIES}) #doctest)
