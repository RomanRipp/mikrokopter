
(cl:in-package :asdf)

(defsystem "mk_msgs-msg"
  :depends-on (:roslisp-msg-protocol :roslisp-utils :geometry_msgs-msg
               :std_msgs-msg
)
  :components ((:file "_package")
    (:file "sensorData" :depends-on ("_package_sensorData"))
    (:file "_package_sensorData" :depends-on ("_package"))
    (:file "motorCommands" :depends-on ("_package_motorCommands"))
    (:file "_package_motorCommands" :depends-on ("_package"))
    (:file "stateData" :depends-on ("_package_stateData"))
    (:file "_package_stateData" :depends-on ("_package"))
    (:file "ARMarker" :depends-on ("_package_ARMarker"))
    (:file "_package_ARMarker" :depends-on ("_package"))
    (:file "ARMarkers" :depends-on ("_package_ARMarkers"))
    (:file "_package_ARMarkers" :depends-on ("_package"))
  ))