# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build

# Include any dependencies generated for this target.
include src/CMakeFiles/wd.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/wd.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/wd.dir/flags.make

src/CMakeFiles/wd.dir/main.c.o: src/CMakeFiles/wd.dir/flags.make
src/CMakeFiles/wd.dir/main.c.o: ../src/main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/wd.dir/main.c.o"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/wd.dir/main.c.o   -c /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/main.c

src/CMakeFiles/wd.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wd.dir/main.c.i"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/main.c > CMakeFiles/wd.dir/main.c.i

src/CMakeFiles/wd.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wd.dir/main.c.s"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/main.c -o CMakeFiles/wd.dir/main.c.s

src/CMakeFiles/wd.dir/weil.c.o: src/CMakeFiles/wd.dir/flags.make
src/CMakeFiles/wd.dir/weil.c.o: ../src/weil.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/CMakeFiles/wd.dir/weil.c.o"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/wd.dir/weil.c.o   -c /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/weil.c

src/CMakeFiles/wd.dir/weil.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wd.dir/weil.c.i"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/weil.c > CMakeFiles/wd.dir/weil.c.i

src/CMakeFiles/wd.dir/weil.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wd.dir/weil.c.s"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/weil.c -o CMakeFiles/wd.dir/weil.c.s

src/CMakeFiles/wd.dir/weil_out.c.o: src/CMakeFiles/wd.dir/flags.make
src/CMakeFiles/wd.dir/weil_out.c.o: ../src/weil_out.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/CMakeFiles/wd.dir/weil_out.c.o"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/wd.dir/weil_out.c.o   -c /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/weil_out.c

src/CMakeFiles/wd.dir/weil_out.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wd.dir/weil_out.c.i"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/weil_out.c > CMakeFiles/wd.dir/weil_out.c.i

src/CMakeFiles/wd.dir/weil_out.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wd.dir/weil_out.c.s"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/weil_out.c -o CMakeFiles/wd.dir/weil_out.c.s

src/CMakeFiles/wd.dir/semaev.c.o: src/CMakeFiles/wd.dir/flags.make
src/CMakeFiles/wd.dir/semaev.c.o: ../src/semaev.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/CMakeFiles/wd.dir/semaev.c.o"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/wd.dir/semaev.c.o   -c /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev.c

src/CMakeFiles/wd.dir/semaev.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wd.dir/semaev.c.i"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev.c > CMakeFiles/wd.dir/semaev.c.i

src/CMakeFiles/wd.dir/semaev.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wd.dir/semaev.c.s"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev.c -o CMakeFiles/wd.dir/semaev.c.s

src/CMakeFiles/wd.dir/semaev_out.c.o: src/CMakeFiles/wd.dir/flags.make
src/CMakeFiles/wd.dir/semaev_out.c.o: ../src/semaev_out.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/CMakeFiles/wd.dir/semaev_out.c.o"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/wd.dir/semaev_out.c.o   -c /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev_out.c

src/CMakeFiles/wd.dir/semaev_out.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wd.dir/semaev_out.c.i"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev_out.c > CMakeFiles/wd.dir/semaev_out.c.i

src/CMakeFiles/wd.dir/semaev_out.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wd.dir/semaev_out.c.s"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev_out.c -o CMakeFiles/wd.dir/semaev_out.c.s

src/CMakeFiles/wd.dir/semaev_masks.c.o: src/CMakeFiles/wd.dir/flags.make
src/CMakeFiles/wd.dir/semaev_masks.c.o: ../src/semaev_masks.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object src/CMakeFiles/wd.dir/semaev_masks.c.o"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/wd.dir/semaev_masks.c.o   -c /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev_masks.c

src/CMakeFiles/wd.dir/semaev_masks.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wd.dir/semaev_masks.c.i"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev_masks.c > CMakeFiles/wd.dir/semaev_masks.c.i

src/CMakeFiles/wd.dir/semaev_masks.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wd.dir/semaev_masks.c.s"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/semaev_masks.c -o CMakeFiles/wd.dir/semaev_masks.c.s

src/CMakeFiles/wd.dir/params.c.o: src/CMakeFiles/wd.dir/flags.make
src/CMakeFiles/wd.dir/params.c.o: ../src/params.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object src/CMakeFiles/wd.dir/params.c.o"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/wd.dir/params.c.o   -c /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/params.c

src/CMakeFiles/wd.dir/params.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wd.dir/params.c.i"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/params.c > CMakeFiles/wd.dir/params.c.i

src/CMakeFiles/wd.dir/params.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wd.dir/params.c.s"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/params.c -o CMakeFiles/wd.dir/params.c.s

src/CMakeFiles/wd.dir/vect_bin.c.o: src/CMakeFiles/wd.dir/flags.make
src/CMakeFiles/wd.dir/vect_bin.c.o: ../src/vect_bin.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object src/CMakeFiles/wd.dir/vect_bin.c.o"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/wd.dir/vect_bin.c.o   -c /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/vect_bin.c

src/CMakeFiles/wd.dir/vect_bin.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/wd.dir/vect_bin.c.i"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/vect_bin.c > CMakeFiles/wd.dir/vect_bin.c.i

src/CMakeFiles/wd.dir/vect_bin.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/wd.dir/vect_bin.c.s"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src/vect_bin.c -o CMakeFiles/wd.dir/vect_bin.c.s

# Object files for target wd
wd_OBJECTS = \
"CMakeFiles/wd.dir/main.c.o" \
"CMakeFiles/wd.dir/weil.c.o" \
"CMakeFiles/wd.dir/weil_out.c.o" \
"CMakeFiles/wd.dir/semaev.c.o" \
"CMakeFiles/wd.dir/semaev_out.c.o" \
"CMakeFiles/wd.dir/semaev_masks.c.o" \
"CMakeFiles/wd.dir/params.c.o" \
"CMakeFiles/wd.dir/vect_bin.c.o"

# External object files for target wd
wd_EXTERNAL_OBJECTS =

../wd: src/CMakeFiles/wd.dir/main.c.o
../wd: src/CMakeFiles/wd.dir/weil.c.o
../wd: src/CMakeFiles/wd.dir/weil_out.c.o
../wd: src/CMakeFiles/wd.dir/semaev.c.o
../wd: src/CMakeFiles/wd.dir/semaev_out.c.o
../wd: src/CMakeFiles/wd.dir/semaev_masks.c.o
../wd: src/CMakeFiles/wd.dir/params.c.o
../wd: src/CMakeFiles/wd.dir/vect_bin.c.o
../wd: src/CMakeFiles/wd.dir/build.make
../wd: src/CMakeFiles/wd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking C executable ../../wd"
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/wd.dir/build: ../wd

.PHONY : src/CMakeFiles/wd.dir/build

src/CMakeFiles/wd.dir/clean:
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src && $(CMAKE_COMMAND) -P CMakeFiles/wd.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/wd.dir/clean

src/CMakeFiles/wd.dir/depend:
	cd /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/src /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src /Users/monika/casspair/Index_calculus/create_benchmarks_S4_github/build/src/CMakeFiles/wd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/wd.dir/depend
