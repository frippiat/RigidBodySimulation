# ----------------------------------------------------------------------------
# CMakeLists.txt
#
#  Created on: 13 Dec 2020
#      Author: Kiwon Um
#        Mail: kiwon.um@telecom-paris.fr
#
# Description: CMake configuration
#
# Copyright 2020-2024 Kiwon Um
#
# The copyright to the computer program(s) herein is the property of Kiwon Um,
# Telecom Paris, France. The program(s) may be used and/or copied only with
# the written permission of Kiwon Um or in accordance with the terms and
# conditions stipulated in the agreement/contract under which the program(s)
# have been supplied.
# ----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.10)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(tpRigid)

add_executable(
  ${PROJECT_NAME}
  src/main.cpp
  # src/Error.cpp # You can include Error.cpp if your system supports OpenGL 4.3 or later; don't forget to replace glad.
  src/Mesh.cpp
  src/ShaderProgram.cpp)

target_sources(${PROJECT_NAME} PRIVATE dep/glad/src/gl.c)
target_include_directories(${PROJECT_NAME} PRIVATE dep/glad/include/)

add_subdirectory(dep/glfw)
target_link_libraries(${PROJECT_NAME} PRIVATE glfw)

add_subdirectory(dep/glm)
target_link_libraries(${PROJECT_NAME} PRIVATE glm::glm-header-only)
target_include_directories(${PROJECT_NAME} PRIVATE dep/glm/../)

target_link_libraries(${PROJECT_NAME} PRIVATE ${CMAKE_DL_LIBS})

add_custom_command(TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_CURRENT_SOURCE_DIR})
