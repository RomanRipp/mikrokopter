FILE(REMOVE_RECURSE
  "../msg_gen"
  "../srv_gen"
  "../src/mk_msgs/msg"
  "../src/mk_msgs/srv"
  "../msg_gen"
  "../srv_gen"
  "CMakeFiles/clean-test-results"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/clean-test-results.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
