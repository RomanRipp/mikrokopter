FILE(REMOVE_RECURSE
  "../msg_gen"
  "../srv_gen"
  "../src/mk_msgs/msg"
  "../src/mk_msgs/srv"
  "../msg_gen"
  "../srv_gen"
  "CMakeFiles/ROSBUILD_gensrv_lisp"
  "../srv_gen/lisp/SetSerialPoti.lisp"
  "../srv_gen/lisp/_package.lisp"
  "../srv_gen/lisp/_package_SetSerialPoti.lisp"
  "../srv_gen/lisp/SetGPSWaypoint.lisp"
  "../srv_gen/lisp/_package.lisp"
  "../srv_gen/lisp/_package_SetGPSWaypoint.lisp"
  "../srv_gen/lisp/GetRendezvousPoint.lisp"
  "../srv_gen/lisp/_package.lisp"
  "../srv_gen/lisp/_package_GetRendezvousPoint.lisp"
  "../srv_gen/lisp/listwp.lisp"
  "../srv_gen/lisp/_package.lisp"
  "../srv_gen/lisp/_package_listwp.lisp"
  "../srv_gen/lisp/SetLanding.lisp"
  "../srv_gen/lisp/_package.lisp"
  "../srv_gen/lisp/_package_SetLanding.lisp"
  "../srv_gen/lisp/SetTakeoff.lisp"
  "../srv_gen/lisp/_package.lisp"
  "../srv_gen/lisp/_package_SetTakeoff.lisp"
  "../srv_gen/lisp/SetNavigation.lisp"
  "../srv_gen/lisp/_package.lisp"
  "../srv_gen/lisp/_package_SetNavigation.lisp"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ROSBUILD_gensrv_lisp.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
