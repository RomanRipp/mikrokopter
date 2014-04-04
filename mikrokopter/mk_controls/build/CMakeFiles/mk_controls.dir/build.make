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
CMAKE_SOURCE_DIR = /home/loki/Documents/ros_workspace/mikrokopter/mk_controls

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/build

# Include any dependencies generated for this target.
include CMakeFiles/mk_controls.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mk_controls.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mk_controls.dir/flags.make

CMakeFiles/mk_controls.dir/src/pidcontroller.o: CMakeFiles/mk_controls.dir/flags.make
CMakeFiles/mk_controls.dir/src/pidcontroller.o: ../src/pidcontroller.cpp
CMakeFiles/mk_controls.dir/src/pidcontroller.o: ../manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/roslang/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/roscpp/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/std_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/geometry_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /home/loki/SwarmSimX/trunk/stacks/ssx_deps/armadillo_wrapper/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/stacks/bullet/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/sensor_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/rosconsole/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/stacks/geometry/angles/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/rospy/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/rostest/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/roswtf/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/share/message_filters/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/stacks/geometry/tf/manifest.xml
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/msg_gen/generated
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/srv_gen/generated
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/stacks/geometry/tf/msg_gen/generated
CMakeFiles/mk_controls.dir/src/pidcontroller.o: /opt/ros/fuerte/stacks/geometry/tf/srv_gen/generated
	$(CMAKE_COMMAND) -E cmake_progress_report /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mk_controls.dir/src/pidcontroller.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -o CMakeFiles/mk_controls.dir/src/pidcontroller.o -c /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/pidcontroller.cpp

CMakeFiles/mk_controls.dir/src/pidcontroller.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mk_controls.dir/src/pidcontroller.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -E /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/pidcontroller.cpp > CMakeFiles/mk_controls.dir/src/pidcontroller.i

CMakeFiles/mk_controls.dir/src/pidcontroller.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mk_controls.dir/src/pidcontroller.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -S /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/pidcontroller.cpp -o CMakeFiles/mk_controls.dir/src/pidcontroller.s

CMakeFiles/mk_controls.dir/src/pidcontroller.o.requires:
.PHONY : CMakeFiles/mk_controls.dir/src/pidcontroller.o.requires

CMakeFiles/mk_controls.dir/src/pidcontroller.o.provides: CMakeFiles/mk_controls.dir/src/pidcontroller.o.requires
	$(MAKE) -f CMakeFiles/mk_controls.dir/build.make CMakeFiles/mk_controls.dir/src/pidcontroller.o.provides.build
.PHONY : CMakeFiles/mk_controls.dir/src/pidcontroller.o.provides

CMakeFiles/mk_controls.dir/src/pidcontroller.o.provides.build: CMakeFiles/mk_controls.dir/src/pidcontroller.o

CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: CMakeFiles/mk_controls.dir/flags.make
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: ../src/bkstepcontroller.cpp
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: ../manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/roslang/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/roscpp/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/std_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/geometry_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /home/loki/SwarmSimX/trunk/stacks/ssx_deps/armadillo_wrapper/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/stacks/bullet/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/sensor_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/rosconsole/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/stacks/geometry/angles/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/rospy/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/rostest/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/roswtf/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/share/message_filters/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/stacks/geometry/tf/manifest.xml
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/msg_gen/generated
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/srv_gen/generated
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/stacks/geometry/tf/msg_gen/generated
CMakeFiles/mk_controls.dir/src/bkstepcontroller.o: /opt/ros/fuerte/stacks/geometry/tf/srv_gen/generated
	$(CMAKE_COMMAND) -E cmake_progress_report /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mk_controls.dir/src/bkstepcontroller.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -o CMakeFiles/mk_controls.dir/src/bkstepcontroller.o -c /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/bkstepcontroller.cpp

CMakeFiles/mk_controls.dir/src/bkstepcontroller.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mk_controls.dir/src/bkstepcontroller.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -E /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/bkstepcontroller.cpp > CMakeFiles/mk_controls.dir/src/bkstepcontroller.i

CMakeFiles/mk_controls.dir/src/bkstepcontroller.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mk_controls.dir/src/bkstepcontroller.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -S /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/bkstepcontroller.cpp -o CMakeFiles/mk_controls.dir/src/bkstepcontroller.s

CMakeFiles/mk_controls.dir/src/bkstepcontroller.o.requires:
.PHONY : CMakeFiles/mk_controls.dir/src/bkstepcontroller.o.requires

CMakeFiles/mk_controls.dir/src/bkstepcontroller.o.provides: CMakeFiles/mk_controls.dir/src/bkstepcontroller.o.requires
	$(MAKE) -f CMakeFiles/mk_controls.dir/build.make CMakeFiles/mk_controls.dir/src/bkstepcontroller.o.provides.build
.PHONY : CMakeFiles/mk_controls.dir/src/bkstepcontroller.o.provides

CMakeFiles/mk_controls.dir/src/bkstepcontroller.o.provides.build: CMakeFiles/mk_controls.dir/src/bkstepcontroller.o

CMakeFiles/mk_controls.dir/src/drivercontroller.o: CMakeFiles/mk_controls.dir/flags.make
CMakeFiles/mk_controls.dir/src/drivercontroller.o: ../src/drivercontroller.cpp
CMakeFiles/mk_controls.dir/src/drivercontroller.o: ../manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/roslang/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/roscpp/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/std_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/geometry_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /home/loki/SwarmSimX/trunk/stacks/ssx_deps/armadillo_wrapper/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/stacks/bullet/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/sensor_msgs/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/rosconsole/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/stacks/geometry/angles/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/rospy/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/rostest/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/roswtf/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/share/message_filters/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/stacks/geometry/tf/manifest.xml
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/msg_gen/generated
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/srv_gen/generated
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/stacks/geometry/tf/msg_gen/generated
CMakeFiles/mk_controls.dir/src/drivercontroller.o: /opt/ros/fuerte/stacks/geometry/tf/srv_gen/generated
	$(CMAKE_COMMAND) -E cmake_progress_report /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mk_controls.dir/src/drivercontroller.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -o CMakeFiles/mk_controls.dir/src/drivercontroller.o -c /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/drivercontroller.cpp

CMakeFiles/mk_controls.dir/src/drivercontroller.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mk_controls.dir/src/drivercontroller.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -E /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/drivercontroller.cpp > CMakeFiles/mk_controls.dir/src/drivercontroller.i

CMakeFiles/mk_controls.dir/src/drivercontroller.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mk_controls.dir/src/drivercontroller.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -DBT_USE_DOUBLE_PRECISION -DBT_EULER_DEFAULT_ZYX -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -S /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/src/drivercontroller.cpp -o CMakeFiles/mk_controls.dir/src/drivercontroller.s

CMakeFiles/mk_controls.dir/src/drivercontroller.o.requires:
.PHONY : CMakeFiles/mk_controls.dir/src/drivercontroller.o.requires

CMakeFiles/mk_controls.dir/src/drivercontroller.o.provides: CMakeFiles/mk_controls.dir/src/drivercontroller.o.requires
	$(MAKE) -f CMakeFiles/mk_controls.dir/build.make CMakeFiles/mk_controls.dir/src/drivercontroller.o.provides.build
.PHONY : CMakeFiles/mk_controls.dir/src/drivercontroller.o.provides

CMakeFiles/mk_controls.dir/src/drivercontroller.o.provides.build: CMakeFiles/mk_controls.dir/src/drivercontroller.o

# Object files for target mk_controls
mk_controls_OBJECTS = \
"CMakeFiles/mk_controls.dir/src/pidcontroller.o" \
"CMakeFiles/mk_controls.dir/src/bkstepcontroller.o" \
"CMakeFiles/mk_controls.dir/src/drivercontroller.o"

# External object files for target mk_controls
mk_controls_EXTERNAL_OBJECTS =

../lib/libmk_controls.so: CMakeFiles/mk_controls.dir/src/pidcontroller.o
../lib/libmk_controls.so: CMakeFiles/mk_controls.dir/src/bkstepcontroller.o
../lib/libmk_controls.so: CMakeFiles/mk_controls.dir/src/drivercontroller.o
../lib/libmk_controls.so: CMakeFiles/mk_controls.dir/build.make
../lib/libmk_controls.so: CMakeFiles/mk_controls.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../lib/libmk_controls.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mk_controls.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mk_controls.dir/build: ../lib/libmk_controls.so
.PHONY : CMakeFiles/mk_controls.dir/build

CMakeFiles/mk_controls.dir/requires: CMakeFiles/mk_controls.dir/src/pidcontroller.o.requires
CMakeFiles/mk_controls.dir/requires: CMakeFiles/mk_controls.dir/src/bkstepcontroller.o.requires
CMakeFiles/mk_controls.dir/requires: CMakeFiles/mk_controls.dir/src/drivercontroller.o.requires
.PHONY : CMakeFiles/mk_controls.dir/requires

CMakeFiles/mk_controls.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mk_controls.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mk_controls.dir/clean

CMakeFiles/mk_controls.dir/depend:
	cd /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/loki/Documents/ros_workspace/mikrokopter/mk_controls /home/loki/Documents/ros_workspace/mikrokopter/mk_controls /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/build /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/build /home/loki/Documents/ros_workspace/mikrokopter/mk_controls/build/CMakeFiles/mk_controls.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mk_controls.dir/depend

