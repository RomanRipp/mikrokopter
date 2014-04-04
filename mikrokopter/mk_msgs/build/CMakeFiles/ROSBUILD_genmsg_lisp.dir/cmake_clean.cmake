FILE(REMOVE_RECURSE
  "../msg_gen"
  "../srv_gen"
  "../src/mk_msgs/msg"
  "../src/mk_msgs/srv"
  "../msg_gen"
  "../srv_gen"
  "CMakeFiles/ROSBUILD_genmsg_lisp"
  "../msg_gen/lisp/sensorData.lisp"
  "../msg_gen/lisp/_package.lisp"
  "../msg_gen/lisp/_package_sensorData.lisp"
  "../msg_gen/lisp/motorCommands.lisp"
  "../msg_gen/lisp/_package.lisp"
  "../msg_gen/lisp/_package_motorCommands.lisp"
  "../msg_gen/lisp/stateData.lisp"
  "../msg_gen/lisp/_package.lisp"
  "../msg_gen/lisp/_package_stateData.lisp"
  "../msg_gen/lisp/ARMarker.lisp"
  "../msg_gen/lisp/_package.lisp"
  "../msg_gen/lisp/_package_ARMarker.lisp"
  "../msg_gen/lisp/ARMarkers.lisp"
  "../msg_gen/lisp/_package.lisp"
  "../msg_gen/lisp/_package_ARMarkers.lisp"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ROSBUILD_genmsg_lisp.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
