cmake_minimum_required(VERSION 3.13)

# Use vcpkg if available
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

project(eos CXX)

option(EOS_BUILD_TEST "Build unit tests" ON)

add_library(eos INTERFACE)
target_compile_features(eos
  INTERFACE
    cxx_std_14
  )
target_compile_options(eos
  INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
  )
target_include_directories(eos
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  )
target_compile_definitions(eos
  INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX _USE_MATH_DEFINES>
  )

if(EOS_BUILD_TEST)
  find_package(GSL REQUIRED)

  # Disable googletest options
  option(BUILD_GMOCK "Build gmock" OFF)
  option(INSTALL_GTEST "Install gtest" OFF)

  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  if(WIN32)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  endif()

  add_subdirectory(third-party/googletest)

  enable_testing()
  add_subdirectory(test)
endif()