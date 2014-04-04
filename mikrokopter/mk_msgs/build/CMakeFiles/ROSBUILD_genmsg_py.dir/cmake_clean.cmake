FILE(REMOVE_RECURSE
  "../msg_gen"
  "../srv_gen"
  "../src/mk_msgs/msg"
  "../src/mk_msgs/srv"
  "../msg_gen"
  "../srv_gen"
  "CMakeFiles/ROSBUILD_genmsg_py"
  "../src/mk_msgs/msg/__init__.py"
  "../src/mk_msgs/msg/_sensorData.py"
  "../src/mk_msgs/msg/_motorCommands.py"
  "../src/mk_msgs/msg/_stateData.py"
  "../src/mk_msgs/msg/_ARMarker.py"
  "../src/mk_msgs/msg/_ARMarkers.py"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ROSBUILD_genmsg_py.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
