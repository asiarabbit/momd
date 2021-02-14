######################################################
# SUNNY project, Anyang Normal University, IMP-CAS
#
# CMakeLists.txt for compilation of library sunny
# Author SUN Yazhou
# Created: 2020/02/02
# Last modified: 2020/10/04
######################################################

# Using ROOT libraries
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

include_directories(inc)
include_directories(../momd/inc)
aux_source_directory(src LIB_SRCS)

add_library(libsunny SHARED ${LIB_SRCS})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set_target_properties(libsunny PROPERTIES OUTPUT_NAME "sunny")

set(LIB_LIST libmomd)
target_link_libraries(libsunny ${LIB_LIST})
