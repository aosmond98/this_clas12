# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.26.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.26.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/osmondal/GitHub/this_clas12

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/osmondal/GitHub/this_clas12/build

# Include any dependencies generated for this target.
include CMakeFiles/clas12_analysis.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/clas12_analysis.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/clas12_analysis.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/clas12_analysis.dir/flags.make

CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o: CMakeFiles/clas12_analysis.dir/flags.make
CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o: /Users/osmondal/GitHub/this_clas12/src/exe/clas12_analysis.cpp
CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o: CMakeFiles/clas12_analysis.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/osmondal/GitHub/this_clas12/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o -MF CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o.d -o CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o -c /Users/osmondal/GitHub/this_clas12/src/exe/clas12_analysis.cpp

CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/osmondal/GitHub/this_clas12/src/exe/clas12_analysis.cpp > CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.i

CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/osmondal/GitHub/this_clas12/src/exe/clas12_analysis.cpp -o CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.s

# Object files for target clas12_analysis
clas12_analysis_OBJECTS = \
"CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o"

# External object files for target clas12_analysis
clas12_analysis_EXTERNAL_OBJECTS =

clas12_analysis: CMakeFiles/clas12_analysis.dir/src/exe/clas12_analysis.cpp.o
clas12_analysis: CMakeFiles/clas12_analysis.dir/build.make
clas12_analysis: libclas12lib.a
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libCore.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libImt.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libRIO.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libNet.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libHist.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libGraf.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libGraf3d.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libGpad.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libROOTDataFrame.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libTree.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libTreePlayer.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libRint.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libPostscript.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libMatrix.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libPhysics.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libMathCore.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libThread.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libMultiProc.so
clas12_analysis: /usr/local/Cellar/root/6.26.06_2/lib/root/libROOTVecOps.so
clas12_analysis: CMakeFiles/clas12_analysis.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/osmondal/GitHub/this_clas12/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable clas12_analysis"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/clas12_analysis.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/clas12_analysis.dir/build: clas12_analysis
.PHONY : CMakeFiles/clas12_analysis.dir/build

CMakeFiles/clas12_analysis.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/clas12_analysis.dir/cmake_clean.cmake
.PHONY : CMakeFiles/clas12_analysis.dir/clean

CMakeFiles/clas12_analysis.dir/depend:
	cd /Users/osmondal/GitHub/this_clas12/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/osmondal/GitHub/this_clas12 /Users/osmondal/GitHub/this_clas12 /Users/osmondal/GitHub/this_clas12/build /Users/osmondal/GitHub/this_clas12/build /Users/osmondal/GitHub/this_clas12/build/CMakeFiles/clas12_analysis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/clas12_analysis.dir/depend

