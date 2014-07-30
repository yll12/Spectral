# Install script for directory: D:/Spectral/Spectral/Spectral/lapack-3.5.0/lapacke

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/LAPACK")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES
    "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lib/liblapacke.dll.a"
    "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lib/liblapacke.lib"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/bin/liblapacke.dll")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/liblapacke.dll" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/liblapacke.dll")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "C:/MinGW/bin/strip.exe" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/liblapacke.dll")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "D:/Spectral/Spectral/Spectral/lapack-3.5.0/lapacke/include/lapacke.h"
    "D:/Spectral/Spectral/Spectral/lapack-3.5.0/lapacke/include/lapacke_config.h"
    "D:/Spectral/Spectral/Spectral/lapack-3.5.0/lapacke/include/lapacke_utils.h"
    "D:/Spectral/Spectral/Spectral/lapack-3.5.0/lapacke/include/lapacke_mangling.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/LAPACK/lib/pkgconfig/lapacke.pc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "C:/Program Files (x86)/LAPACK/lib/pkgconfig" TYPE FILE FILES "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lapacke/lapacke.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lapacke/include/cmake_install.cmake")
  include("D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lapacke/src/cmake_install.cmake")
  include("D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lapacke/utils/cmake_install.cmake")
  include("D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lapacke/example/cmake_install.cmake")

endif()

