FILE(REMOVE_RECURSE
  "../msg_gen"
  "../msg_gen"
  "../src/mk_vision/msg"
  "CMakeFiles/ROSBUILD_genmsg_py"
  "../src/mk_vision/msg/__init__.py"
  "../src/mk_vision/msg/_ARMarker.py"
  "../src/mk_vision/msg/_ARMarkers.py"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ROSBUILD_genmsg_py.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
