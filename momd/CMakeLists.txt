##################################################################################################
# Data Analysis Code Project for the External Target Facility, CSR-HIRFL, @IMP
#
# CMakeLists.txt for compilation of object files and dynamically-linked libraries.
#
# Created Time: 2020/08/30
# Author: SUN Yazhou, asia.rabbit@163.com
# Last Modified: 2017/08/30, SUN Yazhou
##################################################################################################

include_directories(inc)
include_directories(../fci/inc)
aux_source_directory(src LIB_SRCS)

add_library(libmomd SHARED ${LIB_SRCS})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set_target_properties(libmomd PROPERTIES OUTPUT_NAME "momd")
