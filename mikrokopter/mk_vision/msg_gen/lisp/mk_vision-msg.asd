
(cl:in-package :asdf)

(defsystem "mk_vision-msg"
  :depends-on (:roslisp-msg-protocol :roslisp-utils :geometry_msgs-msg
               :std_msgs-msg
)
  :components ((:file "_package")
    (:file "ARMarker" :depends-on ("_package_ARMarker"))
    (:file "_package_ARMarker" :depends-on ("_package"))
    (:file "ARMarkers" :depends-on ("_package_ARMarkers"))
    (:file "_package_ARMarkers" :depends-on ("_package"))
  ))