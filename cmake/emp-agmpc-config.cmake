find_package(emp-tool)
find_package(emp-ot)
find_path(EMP-AGMPC_INCLUDE_DIR emp-agmpc/emp-agmpc.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(EMP-AGMPC DEFAULT_MSG EMP-AGMPC_INCLUDE_DIR)

if(EMP-AGMPC_FOUND)
        set(EMP-AGMPC_INCLUDE_DIRS ${EMP-AGMPC_INCLUDE_DIR} ${EMP-OT_INCLUDE_DIRS})
        set(EMP-AGMPC_LIBRARIES ${EMP-OT_LIBRARIES})
endif()