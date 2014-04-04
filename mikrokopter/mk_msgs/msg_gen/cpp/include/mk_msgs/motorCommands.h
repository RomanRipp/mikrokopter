/* Auto-generated by genmsg_cpp for file /home/loki/Documents/ros_workspace/mikrokopter/mk_msgs/msg/motorCommands.msg */
#ifndef MK_MSGS_MESSAGE_MOTORCOMMANDS_H
#define MK_MSGS_MESSAGE_MOTORCOMMANDS_H
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
struct motorCommands_ {
  typedef motorCommands_<ContainerAllocator> Type;

  motorCommands_()
  : header()
  , pitch(0.0)
  , roll(0.0)
  , yaw(0.0)
  , throttle(0)
  {
  }

  motorCommands_(const ContainerAllocator& _alloc)
  : header(_alloc)
  , pitch(0.0)
  , roll(0.0)
  , yaw(0.0)
  , throttle(0)
  {
  }

  typedef  ::std_msgs::Header_<ContainerAllocator>  _header_type;
   ::std_msgs::Header_<ContainerAllocator>  header;

  typedef double _pitch_type;
  double pitch;

  typedef double _roll_type;
  double roll;

  typedef double _yaw_type;
  double yaw;

  typedef int64_t _throttle_type;
  int64_t throttle;


  typedef boost::shared_ptr< ::mk_msgs::motorCommands_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::mk_msgs::motorCommands_<ContainerAllocator>  const> ConstPtr;
  boost::shared_ptr<std::map<std::string, std::string> > __connection_header;
}; // struct motorCommands
typedef  ::mk_msgs::motorCommands_<std::allocator<void> > motorCommands;

typedef boost::shared_ptr< ::mk_msgs::motorCommands> motorCommandsPtr;
typedef boost::shared_ptr< ::mk_msgs::motorCommands const> motorCommandsConstPtr;


template<typename ContainerAllocator>
std::ostream& operator<<(std::ostream& s, const  ::mk_msgs::motorCommands_<ContainerAllocator> & v)
{
  ros::message_operations::Printer< ::mk_msgs::motorCommands_<ContainerAllocator> >::stream(s, "", v);
  return s;}

} // namespace mk_msgs

namespace ros
{
namespace message_traits
{
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::motorCommands_<ContainerAllocator> > : public TrueType {};
template<class ContainerAllocator> struct IsMessage< ::mk_msgs::motorCommands_<ContainerAllocator>  const> : public TrueType {};
template<class ContainerAllocator>
struct MD5Sum< ::mk_msgs::motorCommands_<ContainerAllocator> > {
  static const char* value() 
  {
    return "fab271cfc8e36c109ccbcbfe60210ac3";
  }

  static const char* value(const  ::mk_msgs::motorCommands_<ContainerAllocator> &) { return value(); } 
  static const uint64_t static_value1 = 0xfab271cfc8e36c10ULL;
  static const uint64_t static_value2 = 0x9ccbcbfe60210ac3ULL;
};

template<class ContainerAllocator>
struct DataType< ::mk_msgs::motorCommands_<ContainerAllocator> > {
  static const char* value() 
  {
    return "mk_msgs/motorCommands";
  }

  static const char* value(const  ::mk_msgs::motorCommands_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator>
struct Definition< ::mk_msgs::motorCommands_<ContainerAllocator> > {
  static const char* value() 
  {
    return "Header header\n\
\n\
# between -pi to +pi\n\
float64 pitch\n\
float64 roll\n\
float64 yaw\n\
\n\
# between 0 to 255\n\
int64 throttle\n\
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

  static const char* value(const  ::mk_msgs::motorCommands_<ContainerAllocator> &) { return value(); } 
};

template<class ContainerAllocator> struct HasHeader< ::mk_msgs::motorCommands_<ContainerAllocator> > : public TrueType {};
template<class ContainerAllocator> struct HasHeader< const ::mk_msgs::motorCommands_<ContainerAllocator> > : public TrueType {};
} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

template<class ContainerAllocator> struct Serializer< ::mk_msgs::motorCommands_<ContainerAllocator> >
{
  template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
  {
    stream.next(m.header);
    stream.next(m.pitch);
    stream.next(m.roll);
    stream.next(m.yaw);
    stream.next(m.throttle);
  }

  ROS_DECLARE_ALLINONE_SERIALIZER;
}; // struct motorCommands_
} // namespace serialization
} // namespace ros

namespace ros
{
namespace message_operations
{

template<class ContainerAllocator>
struct Printer< ::mk_msgs::motorCommands_<ContainerAllocator> >
{
  template<typename Stream> static void stream(Stream& s, const std::string& indent, const  ::mk_msgs::motorCommands_<ContainerAllocator> & v) 
  {
    s << indent << "header: ";
s << std::endl;
    Printer< ::std_msgs::Header_<ContainerAllocator> >::stream(s, indent + "  ", v.header);
    s << indent << "pitch: ";
    Printer<double>::stream(s, indent + "  ", v.pitch);
    s << indent << "roll: ";
    Printer<double>::stream(s, indent + "  ", v.roll);
    s << indent << "yaw: ";
    Printer<double>::stream(s, indent + "  ", v.yaw);
    s << indent << "throttle: ";
    Printer<int64_t>::stream(s, indent + "  ", v.throttle);
  }
};


} // namespace message_operations
} // namespace ros

#endif // MK_MSGS_MESSAGE_MOTORCOMMANDS_H
