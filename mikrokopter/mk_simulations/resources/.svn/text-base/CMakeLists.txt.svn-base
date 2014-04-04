############ Copy media folder to SSX architecture ############
# uncomment this iff media is non-empty
#add_custom_target(media ALL COMMAND mkdir -p ${SSX_HOME}/share/SwarmSimX/media COMMAND cp -R ${PROJECT_SOURCE_DIR}/resources/media/* ${SSX_HOME}/share/SwarmSimX/media/ COMMAND cd ${SSX_HOME}/share/SwarmSimX/media/ && find -name '*.svn' | xargs rm -rf COMMENT "Copying media files to destination: ${SSX_HOME}/share/SwarmSimX/media/")
###############################################################

############ Copy artifacts to SSX architecture ############
# uncomment this iff media is non-empty
add_custom_target(artifacts ALL COMMAND mkdir -p ${SSX_HOME}/share/SwarmSimX/Artifacts COMMAND cp -R ${PROJECT_SOURCE_DIR}/resources/artifacts/* ${SSX_HOME}/share/SwarmSimX/Artifacts/ COMMAND cd ${SSX_HOME}/share/SwarmSimX/Artifacts/ && find -name '*.svn' | xargs rm -rf COMMENT "Copying Artifacts to destination: ${SSX_HOME}/share/SwarmSimX/Artifacts/")
###############################################################

############ Copy configuration to SSX architecture ############
# uncomment this iff media is non-empty
add_custom_target(config ALL COMMAND mkdir -p ${SSX_HOME}/share/SwarmSimX/Config COMMAND cp -R ${PROJECT_SOURCE_DIR}/resources/config/* ${SSX_HOME}/share/SwarmSimX/Config/ COMMAND cd ${SSX_HOME}/share/SwarmSimX/Config/ && find -name '*.svn' | xargs rm -rf COMMENT "Copying Configuration files to destination: ${SSX_HOME}/share/SwarmSimX/Config/")
###############################################################

############ Copy scenes to SSX architecture ############
# uncomment this iff media is non-empty
add_custom_target(scenes ALL COMMAND mkdir -p ${SSX_HOME}/share/SwarmSimX/Scenes COMMAND cp -R ${PROJECT_SOURCE_DIR}/resources/scenes/* ${SSX_HOME}/share/SwarmSimX/Scenes/ COMMAND cd ${SSX_HOME}/share/SwarmSimX/Scenes/ && find -name '*.svn' | xargs rm -rf COMMENT "Copying Scenes to destination: ${SSX_HOME}/share/SwarmSimX/Config/")
###############################################################
