# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/build

# Include any dependencies generated for this target.
include CMakeFiles/earthcoords.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/earthcoords.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/earthcoords.dir/flags.make

CMakeFiles/earthcoords.dir/src/EarthCoords.o: CMakeFiles/earthcoords.dir/flags.make
CMakeFiles/earthcoords.dir/src/EarthCoords.o: ../src/EarthCoords.cpp
CMakeFiles/earthcoords.dir/src/EarthCoords.o: ../manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/roslang/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/roscpp/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/stacks/bullet/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/geometry_msgs/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/sensor_msgs/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/rosconsole/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/stacks/geometry/angles/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/rospy/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/rostest/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/roswtf/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/message_filters/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/stacks/geometry/tf/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /home/loki/SwarmSimX/trunk/stacks/ssx_deps/armadillo_wrapper/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/share/std_msgs/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/manifest.xml
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/stacks/geometry/tf/msg_gen/generated
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /opt/ros/fuerte/stacks/geometry/tf/srv_gen/generated
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/msg_gen/generated
CMakeFiles/earthcoords.dir/src/EarthCoords.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/srv_gen/generated
	$(CMAKE_COMMAND) -E cmake_progress_report /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/earthcoords.dir/src/EarthCoords.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -o CMakeFiles/earthcoords.dir/src/EarthCoords.o -c /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/src/EarthCoords.cpp

CMakeFiles/earthcoords.dir/src/EarthCoords.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/earthcoords.dir/src/EarthCoords.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -E /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/src/EarthCoords.cpp > CMakeFiles/earthcoords.dir/src/EarthCoords.i

CMakeFiles/earthcoords.dir/src/EarthCoords.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/earthcoords.dir/src/EarthCoords.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -S /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/src/EarthCoords.cpp -o CMakeFiles/earthcoords.dir/src/EarthCoords.s

CMakeFiles/earthcoords.dir/src/EarthCoords.o.requires:
.PHONY : CMakeFiles/earthcoords.dir/src/EarthCoords.o.requires

CMakeFiles/earthcoords.dir/src/EarthCoords.o.provides: CMakeFiles/earthcoords.dir/src/EarthCoords.o.requires
	$(MAKE) -f CMakeFiles/earthcoords.dir/build.make CMakeFiles/earthcoords.dir/src/EarthCoords.o.provides.build
.PHONY : CMakeFiles/earthcoords.dir/src/EarthCoords.o.provides

CMakeFiles/earthcoords.dir/src/EarthCoords.o.provides.build: CMakeFiles/earthcoords.dir/src/EarthCoords.o

# Object files for target earthcoords
earthcoords_OBJECTS = \
"CMakeFiles/earthcoords.dir/src/EarthCoords.o"

# External object files for target earthcoords
earthcoords_EXTERNAL_OBJECTS =

../lib/libearthcoords.so: CMakeFiles/earthcoords.dir/src/EarthCoords.o
../lib/libearthcoords.so: CMakeFiles/earthcoords.dir/build.make
../lib/libearthcoords.so: CMakeFiles/earthcoords.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../lib/libearthcoords.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/earthcoords.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/earthcoords.dir/build: ../lib/libearthcoords.so
.PHONY : CMakeFiles/earthcoords.dir/build

CMakeFiles/earthcoords.dir/requires: CMakeFiles/earthcoords.dir/src/EarthCoords.o.requires
.PHONY : CMakeFiles/earthcoords.dir/requires

CMakeFiles/earthcoords.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/earthcoords.dir/cmake_clean.cmake
.PHONY : CMakeFiles/earthcoords.dir/clean

CMakeFiles/earthcoords.dir/depend:
	cd /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/build /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/build /home/loki/Documents/ros_workspace/mikrokopter/mk_behaviors/build/CMakeFiles/earthcoords.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/earthcoords.dir/depend

