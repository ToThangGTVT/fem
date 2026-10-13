cmake_minimum_required(VERSION 3.20)
project(FEM2D VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(FEM_BUILD_GUI "Xây dựng giao diện SDL3 + Dear ImGui" ON)
option(FEM_BUILD_TESTS "Xây dựng bộ kiểm tra lõi tính toán" ON)

# ---------------------------------------------------------------- lõi tính toán (không phụ thuộc giao diện)
add_library(femcore STATIC
  src/core/Model.cpp
  src/core/Skyline.cpp
  src/core/Elements.cpp
  src/core/Analysis.cpp
  src/core/Library.cpp
  src/core/Mesher.cpp
  src/core/ModelIO.cpp
  src/core/Examples.cpp
  src/core/PostProcess.cpp
  src/core/Eigen.cpp
  src/core/Kinematics.cpp
  src/core/SectionGeom.cpp
  src/core/ForceMethod.cpp
  src/core/Dynamics.cpp
  src/core/Design.cpp
  src/core/Dxf.cpp
  src/core/Report.cpp
  src/core/PlaneElements.cpp
)
target_include_directories(femcore PUBLIC src)
if(MSVC)
  target_compile_options(femcore PRIVATE /W4 /utf-8)
  target_compile_definitions(femcore PUBLIC _USE_MATH_DEFINES)
else()
  target_compile_options(femcore PRIVATE -Wall -Wextra)
endif()

# ---------------------------------------------------------------- kiểm tra
if(FEM_BUILD_TESTS)
  enable_testing()
  add_executable(fem_tests tests/test_core.cpp)
  target_link_libraries(fem_tests PRIVATE femcore)
  add_test(NAME core COMMAND fem_tests)
endif()

# ---------------------------------------------------------------- giao diện
if(FEM_BUILD_GUI)
  include(FetchContent)

  find_package(SDL3 CONFIG QUIET)
  if(NOT SDL3_FOUND)
    message(STATUS "Không tìm thấy SDL3 trong hệ thống — tải mã nguồn SDL3")
    FetchContent_Declare(SDL3
      GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
      GIT_TAG release-3.2.8
      GIT_SHALLOW TRUE)
    set(SDL_SHARED OFF CACHE BOOL "" FORCE)
    set(SDL_STATIC ON CACHE BOOL "" FORCE)
    set(SDL_TEST_LIBRARY OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(SDL3)
  endif()

  set(FEM_IMGUI_DIR "" CACHE PATH "Thư mục mã nguồn Dear ImGui có sẵn (để trống: tự tải)")
  if(FEM_IMGUI_DIR)
    set(imgui_SOURCE_DIR ${FEM_IMGUI_DIR})
  else()
    FetchContent_Declare(imgui
      GIT_REPOSITORY https://github.com/ocornut/imgui.git
      GIT_TAG v1.92.9b-docking
      GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(imgui)
  endif()

  add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_sdl3.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_sdlrenderer3.cpp
    ${imgui_SOURCE_DIR}/misc/cpp/imgui_stdlib.cpp
  )
  target_include_directories(imgui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends ${imgui_SOURCE_DIR}/misc/cpp)
  target_link_libraries(imgui PUBLIC SDL3::SDL3)

  add_executable(fem2d
    src/app/main.cpp
    src/app/App.cpp
    src/app/Panels.cpp
    src/app/View.cpp
    src/app/Commands.cpp
    src/app/Ribbon.cpp
    src/app/SectionUI.cpp
    src/app/Textbook.cpp
    src/app/DynamicsUI.cpp
    src/app/DesignUI.cpp
  )
  target_link_libraries(fem2d PRIVATE femcore imgui)
  if(MSVC)
    target_compile_options(fem2d PRIVATE /utf-8)
  else()
    target_compile_options(fem2d PRIVATE -Wall -Wextra)
  endif()
  if(WIN32)
    set_target_properties(fem2d PROPERTIES WIN32_EXECUTABLE ON)
  endif()
endif()
