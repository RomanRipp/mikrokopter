FILE(REMOVE_RECURSE
  "../msg_gen"
  "../srv_gen"
  "../src/mk_msgs/msg"
  "../src/mk_msgs/srv"
  "../msg_gen"
  "../srv_gen"
  "CMakeFiles/ROSBUILD_gensrv_py"
  "../src/mk_msgs/srv/__init__.py"
  "../src/mk_msgs/srv/_SetSerialPoti.py"
  "../src/mk_msgs/srv/_SetGPSWaypoint.py"
  "../src/mk_msgs/srv/_GetRendezvousPoint.py"
  "../src/mk_msgs/srv/_listwp.py"
  "../src/mk_msgs/srv/_SetLanding.py"
  "../src/mk_msgs/srv/_SetTakeoff.py"
  "../src/mk_msgs/srv/_SetNavigation.py"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ROSBUILD_gensrv_py.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
