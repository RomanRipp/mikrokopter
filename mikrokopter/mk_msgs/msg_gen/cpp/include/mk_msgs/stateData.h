/* Auto-generated by genmsg_cpp for file /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/msg/stateData.msg */
#ifndef MK_MSGS_MESSAGE_STATEDATA_H
#define MK_MSGS_MESSAGE_STATEDATA_H
#include <string>
#include <vector>
#include <map>
#include <ostream>
#include "ros/serialization.h"
#include "ros/builtin_message_traits.h"
#include "ros/message_operations.h"
#include "ros/time.h"

#include "ros/macros.h"

#include "ros/assert.h"

#include "std_msgs/Header.h"

namespace mk_msgs
{
template <class ContainerAllocator>
struct stateData_ {
  typedef stateData_<ContainerAllocator> Type;

  stateData_()
  : header()
  , x(0.0)
  , y(0.0)
  , z(0.0)
  , roll(0.0)
  , pitch(0.0)
  , yaw(0.0)
  {
  }

  stateData_(const ContainerAllocator& _alloc)
  : header(_alloc)
  , x(0.0)
  , y(0.0)
  , z(0.0)
  , roll(0.0)
  , pitch(0.0)
  , yaw(0.0)
  {
  }

  typedef  ::std_msgs::Header_<ContainerAllocator>  _header_type;
   ::std_msgs::Header_<ContainerAllocator>  header;

  typedef float _x_type;
  float x;

  typedef float _y_type;
  float y;

  typedef float _z_type;
  float z;

  typedef float _roll_type;
  float roll;

  typedef float _pitch_type;
  float pitch;

  typedef float _yaw_type;
  float yaw;


  typedef boost::shared_ptr< ::mk_msgs::stateData_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::mk_msgs::stateData_<ContainerAllocator>  const> ConstPtr;
  boost::shared_ptr<std::map<std::string, std::string> > __connection_header;
}; // struct stateData
typedef  ::mk_msgs::stateData_<std::allocator<void> > stateData;

typedef boost::shared_ptr< ::mk_msgs::stateData> stateDataPtr;
typedef boost::shared_ptr< ::mk_msgs::stateData const> stateDataConstPtr;


template<typename ContainerAllocator>
std::ostream& operator<<(std::ostream& s, const  ::mk_msgs::stateData_<ContainerAllocator> & v)
{
  ros::message_operations::Printer< ::mk_msgs::stateData_<ContainerAllocator> >::stream(s, "", v);
  return s;}

} // namespace mk_msgs

namespace ros
{
namespace message_traits
{
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::stateData_<ContainerAllocator> > : public TrueType {};
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::stateData_<ContainerAllocator>  const> : public TrueType {};
template<class ContainerAllocator>
struct MD5Sum< ::mk_msgs::stateData_<ContainerAllocator> > {
  static const char* value() 
  {
    return "047123678a922f9085d012f623f64b55";
  }

  static const char* value(const  ::mk_msgs::stateData_<ContainerAllocator> &) { return value(); } 
  static const uint64_t static_value1 = 0x047123678a922f90ULL;
  static const uint64_t static_value2 = 0x85d012f623f64b55ULL;
};

template<class ContainerAllocator>
struct DataType< ::mk_msgs::stateData_<ContainerAllocator> > {
  static const char* value() 
  {
    return "mk_msgs/stateData";
  }

  static const char* value(const  ::mk_msgs::stateData_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator>
struct Definition< ::mk_msgs::stateData_<ContainerAllocator> > {
  static const char* value() 
  {
    return "Header header\n\
\n\
float32 x\n\
float32 y\n\
float32 z\n\
\n\
float32 roll\n\
float32 pitch\n\
float32 yaw\n\
\n\
\n\
================================================================================\n\
MSG: std_msgs/Header\n\
# Standard metadata for higher-level stamped data types.\n\
# This is generally used to communicate timestamped data \n\
# in a particular coordinate frame.\n\
# \n\
# sequence ID: consecutively increasing ID \n\
uint32 seq\n\
#Two-integer timestamp that is expressed as:\n\
# * stamp.secs: seconds (stamp_secs) since epoch\n\
# * stamp.nsecs: nanoseconds since stamp_secs\n\
# time-handling sugar is provided by the client library\n\
time stamp\n\
#Frame this data is associated with\n\
# 0: no frame\n\
# 1: global frame\n\
string frame_id\n\
\n\
";
  }

  static const char* value(const  ::mk_msgs::stateData_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator> struct HasHeader< ::mk_msgs::stateData_<ContainerAllocator> > : public TrueType {};
template<class ContainerAllocator> struct HasHeader< const ::mk_msgs::stateData_<ContainerAllocator> > : public TrueType {};
} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

template<class ContainerAllocator> struct Serializer< ::mk_msgs::stateData_<ContainerAllocator> >
{
  template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
  {
    stream.next(m.header);
    stream.next(m.x);
    stream.next(m.y);
    stream.next(m.z);
    stream.next(m.roll);
    stream.next(m.pitch);
    stream.next(m.yaw);
  }

  ROS_DECLARE_ALLINONE_SERIALIZER;
}; // struct stateData_
} // namespace serialization
} // namespace ros

namespace ros
{
namespace message_operations
{

template<class ContainerAllocator>
struct Printer< ::mk_msgs::stateData_<ContainerAllocator> >
{
  template<typename Stream> static void stream(Stream& s, const std::string& indent, const  ::mk_msgs::stateData_<ContainerAllocator> & v) 
  {
    s << indent << "header: ";
s << std::endl;
    Printer< ::std_msgs::Header_<ContainerAllocator> >::stream(s, indent + "  ", v.header);
    s << indent << "x: ";
    Printer<float>::stream(s, indent + "  ", v.x);
    s << indent << "y: ";
    Printer<float>::stream(s, indent + "  ", v.y);
    s << indent << "z: ";
    Printer<float>::stream(s, indent + "  ", v.z);
    s << indent << "roll: ";
    Printer<float>::stream(s, indent + "  ", v.roll);
    s << indent << "pitch: ";
    Printer<float>::stream(s, indent + "  ", v.pitch);
    s << indent << "yaw: ";
    Printer<float>::stream(s, indent + "  ", v.yaw);
  }
};


} // namespace message_operations
} // namespace ros

#endif // MK_MSGS_MESSAGE_STATEDATA_H

