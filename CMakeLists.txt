#----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.9)
project(LEPAM)

#----------------------------------------------------------------------------

find_package(ROOT 6.30 CONFIG REQUIRED)

#----------------------------------------------------------------------------

find_package(OpenMP)

#----------------------------------------------------------------------------

include_directories(${PROJECT_SOURCE_DIR}/include
		    ${PROJECT_SOURCE_DIR}/shared/include
        ${ROOT_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
		  ${PROJECT_SOURCE_DIR}/shared/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
		  ${PROJECT_SOURCE_DIR}/shared/include/*.hh)

#----------------------------------------------------------------------------

separate_arguments(ROOT_EXE_LINKER_FLAGS)

#----------------------------------------------------------------------------

set(CMAKE_CXX_COMPILER g++)

#----------------------------------------------------------------------------

add_compile_options(-std=c++17 -Wall -O2)

#----------------------------------------------------------------------------

add_executable(LEPAM Mainfile.cc ${sources} ${headers})
target_link_libraries(LEPAM PUBLIC ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS} ROOT::Geom OpenMP::OpenMP_CXX)

#----------------------------------------------------------------------------

set(LEPAM_SCRIPTS run.param)

foreach(_script ${LEPAM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------

install(TARGETS LEPAM DESTINATION ${PROJECT_SOURCE_DIR})

#----------------------------------------------------------------------------

add_custom_target("uninstall" COMMENT "Uninstall installed files")
add_custom_command(
    TARGET "uninstall"
    POST_BUILD
    COMMENT "Uninstall files with install_manifest.txt"
    COMMAND xargs rm -vf < install_manifest.txt || echo Nothing in
            install_manifest.txt to be uninstalled!
)

#----------------------------------------------------------------------------