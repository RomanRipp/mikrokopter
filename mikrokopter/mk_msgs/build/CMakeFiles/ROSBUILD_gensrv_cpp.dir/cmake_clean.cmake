FILE(REMOVE_RECURSE
  "../msg_gen"
  "../srv_gen"
  "../src/mk_msgs/msg"
  "../src/mk_msgs/srv"
  "../msg_gen"
  "../srv_gen"
  "CMakeFiles/ROSBUILD_gensrv_cpp"
  "../srv_gen/cpp/include/mk_msgs/SetSerialPoti.h"
  "../srv_gen/cpp/include/mk_msgs/SetGPSWaypoint.h"
  "../srv_gen/cpp/include/mk_msgs/GetRendezvousPoint.h"
  "../srv_gen/cpp/include/mk_msgs/listwp.h"
  "../srv_gen/cpp/include/mk_msgs/SetLanding.h"
  "../srv_gen/cpp/include/mk_msgs/SetTakeoff.h"
  "../srv_gen/cpp/include/mk_msgs/SetNavigation.h"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ROSBUILD_gensrv_cpp.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
