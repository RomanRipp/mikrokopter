FILE(REMOVE_RECURSE
  "../msg_gen"
  "../srv_gen"
  "../src/mk_msgs/msg"
  "../src/mk_msgs/srv"
  "../msg_gen"
  "../srv_gen"
  "CMakeFiles/ROSBUILD_genmsg_cpp"
  "../msg_gen/cpp/include/mk_msgs/sensorData.h"
  "../msg_gen/cpp/include/mk_msgs/motorCommands.h"
  "../msg_gen/cpp/include/mk_msgs/stateData.h"
  "../msg_gen/cpp/include/mk_msgs/ARMarker.h"
  "../msg_gen/cpp/include/mk_msgs/ARMarkers.h"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ROSBUILD_genmsg_cpp.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
