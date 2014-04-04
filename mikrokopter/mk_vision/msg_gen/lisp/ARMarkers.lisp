; Auto-generated. Do not edit!


(cl:in-package mk_vision-msg)


;//! \htmlinclude ARMarkers.msg.html

(cl:defclass <ARMarkers> (roslisp-msg-protocol:ros-message)
  ((header
    :reader header
    :initarg :header
    :type std_msgs-msg:Header
    :initform (cl:make-instance 'std_msgs-msg:Header))
   (markers
    :reader markers
    :initarg :markers
    :type (cl:vector mk_vision-msg:ARMarker)
   :initform (cl:make-array 0 :element-type 'mk_vision-msg:ARMarker :initial-element (cl:make-instance 'mk_vision-msg:ARMarker))))
)

(cl:defclass ARMarkers (<ARMarkers>)
  ())

(cl:defmethod cl:initialize-instance :after ((m <ARMarkers>) cl:&rest args)
  (cl:declare (cl:ignorable args))
  (cl:unless (cl:typep m 'ARMarkers)
    (roslisp-msg-protocol:msg-deprecation-warning "using old message class name mk_vision-msg:<ARMarkers> is deprecated: use mk_vision-msg:ARMarkers instead.")))

(cl:ensure-generic-function 'header-val :lambda-list '(m))
(cl:defmethod header-val ((m <ARMarkers>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_vision-msg:header-val is deprecated.  Use mk_vision-msg:header instead.")
  (header m))

(cl:ensure-generic-function 'markers-val :lambda-list '(m))
(cl:defmethod markers-val ((m <ARMarkers>))
  (roslisp-msg-protocol:msg-deprecation-warning "Using old-style slot reader mk_vision-msg:markers-val is deprecated.  Use mk_vision-msg:markers instead.")
  (markers m))
(cl:defmethod roslisp-msg-protocol:serialize ((msg <ARMarkers>) ostream)
  "Serializes a message object of type '<ARMarkers>"
  (roslisp-msg-protocol:serialize (cl:slot-value msg 'header) ostream)
  (cl:let ((__ros_arr_len (cl:length (cl:slot-value msg 'markers))))
    (cl:write-byte (cl:ldb (cl:byte 8 0) __ros_arr_len) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 8) __ros_arr_len) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 16) __ros_arr_len) ostream)
    (cl:write-byte (cl:ldb (cl:byte 8 24) __ros_arr_len) ostream))
  (cl:map cl:nil #'(cl:lambda (ele) (roslisp-msg-protocol:serialize ele ostream))
   (cl:slot-value msg 'markers))
)
(cl:defmethod roslisp-msg-protocol:deserialize ((msg <ARMarkers>) istream)
  "Deserializes a message object of type '<ARMarkers>"
  (roslisp-msg-protocol:deserialize (cl:slot-value msg 'header) istream)
  (cl:let ((__ros_arr_len 0))
    (cl:setf (cl:ldb (cl:byte 8 0) __ros_arr_len) (cl:read-byte istream))
    (cl:setf (cl:ldb (cl:byte 8 8) __ros_arr_len) (cl:read-byte istream))
    (cl:setf (cl:ldb (cl:byte 8 16) __ros_arr_len) (cl:read-byte istream))
    (cl:setf (cl:ldb (cl:byte 8 24) __ros_arr_len) (cl:read-byte istream))
  (cl:setf (cl:slot-value msg 'markers) (cl:make-array __ros_arr_len))
  (cl:let ((vals (cl:slot-value msg 'markers)))
    (cl:dotimes (i __ros_arr_len)
    (cl:setf (cl:aref vals i) (cl:make-instance 'mk_vision-msg:ARMarker))
  (roslisp-msg-protocol:deserialize (cl:aref vals i) istream))))
  msg
)
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql '<ARMarkers>)))
  "Returns string type for a message object of type '<ARMarkers>"
  "mk_vision/ARMarkers")
(cl:defmethod roslisp-msg-protocol:ros-datatype ((msg (cl:eql 'ARMarkers)))
  "Returns string type for a message object of type 'ARMarkers"
  "mk_vision/ARMarkers")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql '<ARMarkers>)))
  "Returns md5sum for a message object of type '<ARMarkers>"
  "b35e1e178a9cd7039dbb63cf2764131a")
(cl:defmethod roslisp-msg-protocol:md5sum ((type (cl:eql 'ARMarkers)))
  "Returns md5sum for a message object of type 'ARMarkers"
  "b35e1e178a9cd7039dbb63cf2764131a")
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql '<ARMarkers>)))
  "Returns full string definition for message of type '<ARMarkers>"
  (cl:format cl:nil "Header header~%ARMarker[] markers~%~%================================================================================~%MSG: std_msgs/Header~%# Standard metadata for higher-level stamped data types.~%# This is generally used to communicate timestamped data ~%# in a particular coordinate frame.~%# ~%# sequence ID: consecutively increasing ID ~%uint32 seq~%#Two-integer timestamp that is expressed as:~%# * stamp.secs: seconds (stamp_secs) since epoch~%# * stamp.nsecs: nanoseconds since stamp_secs~%# time-handling sugar is provided by the client library~%time stamp~%#Frame this data is associated with~%# 0: no frame~%# 1: global frame~%string frame_id~%~%================================================================================~%MSG: mk_vision/ARMarker~%Header header~%uint32 id~%geometry_msgs/PoseWithCovariance pose~%uint32 confidence~%~%================================================================================~%MSG: geometry_msgs/PoseWithCovariance~%# This represents a pose in free space with uncertainty.~%~%Pose pose~%~%# Row-major representation of the 6x6 covariance matrix~%# The orientation parameters use a fixed-axis representation.~%# In order, the parameters are:~%# (x, y, z, rotation about X axis, rotation about Y axis, rotation about Z axis)~%float64[36] covariance~%~%================================================================================~%MSG: geometry_msgs/Pose~%# A representation of pose in free space, composed of postion and orientation. ~%Point position~%Quaternion orientation~%~%================================================================================~%MSG: geometry_msgs/Point~%# This contains the position of a point in free space~%float64 x~%float64 y~%float64 z~%~%================================================================================~%MSG: geometry_msgs/Quaternion~%# This represents an orientation in free space in quaternion form.~%~%float64 x~%float64 y~%float64 z~%float64 w~%~%~%"))
(cl:defmethod roslisp-msg-protocol:message-definition ((type (cl:eql 'ARMarkers)))
  "Returns full string definition for message of type 'ARMarkers"
  (cl:format cl:nil "Header header~%ARMarker[] markers~%~%================================================================================~%MSG: std_msgs/Header~%# Standard metadata for higher-level stamped data types.~%# This is generally used to communicate timestamped data ~%# in a particular coordinate frame.~%# ~%# sequence ID: consecutively increasing ID ~%uint32 seq~%#Two-integer timestamp that is expressed as:~%# * stamp.secs: seconds (stamp_secs) since epoch~%# * stamp.nsecs: nanoseconds since stamp_secs~%# time-handling sugar is provided by the client library~%time stamp~%#Frame this data is associated with~%# 0: no frame~%# 1: global frame~%string frame_id~%~%================================================================================~%MSG: mk_vision/ARMarker~%Header header~%uint32 id~%geometry_msgs/PoseWithCovariance pose~%uint32 confidence~%~%================================================================================~%MSG: geometry_msgs/PoseWithCovariance~%# This represents a pose in free space with uncertainty.~%~%Pose pose~%~%# Row-major representation of the 6x6 covariance matrix~%# The orientation parameters use a fixed-axis representation.~%# In order, the parameters are:~%# (x, y, z, rotation about X axis, rotation about Y axis, rotation about Z axis)~%float64[36] covariance~%~%================================================================================~%MSG: geometry_msgs/Pose~%# A representation of pose in free space, composed of postion and orientation. ~%Point position~%Quaternion orientation~%~%================================================================================~%MSG: geometry_msgs/Point~%# This contains the position of a point in free space~%float64 x~%float64 y~%float64 z~%~%================================================================================~%MSG: geometry_msgs/Quaternion~%# This represents an orientation in free space in quaternion form.~%~%float64 x~%float64 y~%float64 z~%float64 w~%~%~%"))
(cl:defmethod roslisp-msg-protocol:serialization-length ((msg <ARMarkers>))
  (cl:+ 0
     (roslisp-msg-protocol:serialization-length (cl:slot-value msg 'header))
     4 (cl:reduce #'cl:+ (cl:slot-value msg 'markers) :key #'(cl:lambda (ele) (cl:declare (cl:ignorable ele)) (cl:+ (roslisp-msg-protocol:serialization-length ele))))
))
(cl:defmethod roslisp-msg-protocol:ros-message-to-list ((msg <ARMarkers>))
  "Converts a ROS message object to a list"
  (cl:list 'ARMarkers
    (cl:cons ':header (header msg))
    (cl:cons ':markers (markers msg))
))
