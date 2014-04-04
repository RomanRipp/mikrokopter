FILE(REMOVE_RECURSE
  "../msg_gen"
  "../msg_gen"
  "../src/mk_vision/msg"
  "CMakeFiles/ROSBUILD_genmsg_lisp"
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
