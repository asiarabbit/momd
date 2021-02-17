######################################################
# SUNNY project, Anyang Normal University, IMP-CAS
#
# CMakeLists.txt for user-composed programs
# Author: SUN Yazhou, asia.rabbit@163.com
# Created: 2020/08/30
# Last modified: 2020/10/04, SUN Yazhou
######################################################

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS}) # ROOT header path
include_directories(${PROJECT_SOURCE_DIR}/fci/inc) # SUNNY header path
include_directories(${PROJECT_SOURCE_DIR}/momd/inc)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIB_LIST libmomd libsunny ${ROOT_LIBRARIES})

add_executable(te test.cxx)
target_link_libraries(te ${LIB_LIST})

# many-body basis generation
add_executable(mbsd mbsd.cxx)
target_link_libraries(mbsd ${LIB_LIST})

# hamiltonian matrix establishment
add_executable(ham hamiltonian.cxx)
target_link_libraries(ham ${LIB_LIST})
