# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files (x86)\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files (x86)\CMake\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\Spectral\Spectral\Spectral\lapack-3.5.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\Spectral\Spectral\Spectral\lapack-3.5.0\BUILD

# Include any dependencies generated for this target.
include BLAS/TESTING/CMakeFiles/xblat3d.dir/depend.make

# Include the progress variables for this target.
include BLAS/TESTING/CMakeFiles/xblat3d.dir/progress.make

# Include the compile flags for this target's objects.
include BLAS/TESTING/CMakeFiles/xblat3d.dir/flags.make

BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj: BLAS/TESTING/CMakeFiles/xblat3d.dir/flags.make
BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj: ../BLAS/TESTING/dblat3.f
	$(CMAKE_COMMAND) -E cmake_progress_report D:\Spectral\Spectral\Spectral\lapack-3.5.0\BUILD\CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj"
	cd /d D:\Spectral\Spectral\Spectral\lapack-3.5.0\BUILD\BLAS\TESTING && C:\MinGW\bin\gfortran.exe  $(Fortran_DEFINES) $(Fortran_FLAGS) -c D:\Spectral\Spectral\Spectral\lapack-3.5.0\BLAS\TESTING\dblat3.f -o CMakeFiles\xblat3d.dir\dblat3.f.obj

BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj.requires:
.PHONY : BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj.requires

BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj.provides: BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj.requires
	$(MAKE) -f BLAS\TESTING\CMakeFiles\xblat3d.dir\build.make BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj.provides.build
.PHONY : BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj.provides

BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj.provides.build: BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj

# Object files for target xblat3d
xblat3d_OBJECTS = \
"CMakeFiles/xblat3d.dir/dblat3.f.obj"

# External object files for target xblat3d
xblat3d_EXTERNAL_OBJECTS =

bin/xblat3d.exe: BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj
bin/xblat3d.exe: BLAS/TESTING/CMakeFiles/xblat3d.dir/build.make
bin/xblat3d.exe: lib/libblas.dll.a
bin/xblat3d.exe: BLAS/TESTING/CMakeFiles/xblat3d.dir/objects1.rsp
bin/xblat3d.exe: BLAS/TESTING/CMakeFiles/xblat3d.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking Fortran executable ..\..\bin\xblat3d.exe"
	cd /d D:\Spectral\Spectral\Spectral\lapack-3.5.0\BUILD\BLAS\TESTING && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\xblat3d.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
BLAS/TESTING/CMakeFiles/xblat3d.dir/build: bin/xblat3d.exe
.PHONY : BLAS/TESTING/CMakeFiles/xblat3d.dir/build

BLAS/TESTING/CMakeFiles/xblat3d.dir/requires: BLAS/TESTING/CMakeFiles/xblat3d.dir/dblat3.f.obj.requires
.PHONY : BLAS/TESTING/CMakeFiles/xblat3d.dir/requires

BLAS/TESTING/CMakeFiles/xblat3d.dir/clean:
	cd /d D:\Spectral\Spectral\Spectral\lapack-3.5.0\BUILD\BLAS\TESTING && $(CMAKE_COMMAND) -P CMakeFiles\xblat3d.dir\cmake_clean.cmake
.PHONY : BLAS/TESTING/CMakeFiles/xblat3d.dir/clean

BLAS/TESTING/CMakeFiles/xblat3d.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\Spectral\Spectral\Spectral\lapack-3.5.0 D:\Spectral\Spectral\Spectral\lapack-3.5.0\BLAS\TESTING D:\Spectral\Spectral\Spectral\lapack-3.5.0\BUILD D:\Spectral\Spectral\Spectral\lapack-3.5.0\BUILD\BLAS\TESTING D:\Spectral\Spectral\Spectral\lapack-3.5.0\BUILD\BLAS\TESTING\CMakeFiles\xblat3d.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : BLAS/TESTING/CMakeFiles/xblat3d.dir/depend

