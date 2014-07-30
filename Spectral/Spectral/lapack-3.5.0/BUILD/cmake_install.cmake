# Install script for directory: D:/Spectral/Spectral/Spectral/lapack-3.5.0

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
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "C:/Program Files (x86)/LAPACK/lib/pkgconfig/lapack.pc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "C:/Program Files (x86)/LAPACK/lib/pkgconfig" TYPE FILE FILES "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lapack.pc")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/lapack-3.4.2" TYPE FILE FILES
    "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/CMakeFiles/lapack-config.cmake"
    "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lapack-config-version.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/lapack-3.4.2/lapack-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/lapack-3.4.2/lapack-targets.cmake"
         "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/CMakeFiles/Export/lib/cmake/lapack-3.4.2/lapack-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/lapack-3.4.2/lapack-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/lapack-3.4.2/lapack-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/lapack-3.4.2" TYPE FILE FILES "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/CMakeFiles/Export/lib/cmake/lapack-3.4.2/lapack-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/lapack-3.4.2" TYPE FILE FILES "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/CMakeFiles/Export/lib/cmake/lapack-3.4.2/lapack-targets-noconfig.cmake")
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/BLAS/cmake_install.cmake")
  include("D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/SRC/cmake_install.cmake")
  include("D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/Testing/cmake_install.cmake")
  include("D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/lapacke/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

file(WRITE "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/${CMAKE_INSTALL_MANIFEST}" "")
foreach(file ${CMAKE_INSTALL_MANIFEST_FILES})
  file(APPEND "D:/Spectral/Spectral/Spectral/lapack-3.5.0/BUILD/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
endforeach()
