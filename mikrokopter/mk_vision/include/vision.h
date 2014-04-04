#ifndef AR_POSE_AR_SINGLE_H
#define AR_POSE_AR_SINGLE_H

#include <string.h>
#include <stdarg.h>

#include <AR/config.h>
#include <AR/param.h>
#include <AR/ar.h>
#include <AR/video.h>

#include <ros/ros.h>
#include <ros/package.h>
#include <ros/console.h>
#include <geometry_msgs/TransformStamped.h>
#include <tf/transform_broadcaster.h>
#include <image_transport/image_transport.h>
#include <sensor_msgs/CameraInfo.h>
#include <visualization_msgs/Marker.h>
#include <resource_retriever/retriever.h>

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <cv_bridge/CvBridge.h>

#include <ar_pose/ARMarker.h>

const std::string cameraImageTopic_ = "/usb_cam/image_raw";
const std::string cameraInfoTopic_  = "/usb_cam/camera_info";

const double AR_TO_ROS = 0.001;

namespace ar_pose
{
  class ARSinglePublisher
  {
  public:
    ARSinglePublisher (ros::NodeHandle & n, double x, double y, double z, double rx, double ry, double rz);
    ~ARSinglePublisher (void);
  private:
    void arInit ();
    void getTransformationCallback (const sensor_msgs::ImageConstPtr &);
    void camInfoCallback (const sensor_msgs::CameraInfoConstPtr &);
	IplImage* processFrame(IplImage*);

    ros::NodeHandle n_;
    ros::Subscriber sub_;
    tf::TransformBroadcaster broadcaster_;
    ros::Publisher arMarkerPub_;

    btTransform transform_;

    ar_pose::ARMarker ar_pose_marker_;
    image_transport::ImageTransport it_;
    image_transport::Subscriber cam_sub_;
    sensor_msgs::CvBridge bridge_;
    sensor_msgs::CameraInfo cam_info_;

    // **** parameters

    //std::string cameraFrame_;
    std::string markerFrame_;
    bool publishTf_;
    bool publishVisualMarkers_;
    bool useHistory_;
    int threshold_;
    double markerWidth_;        // Size of the AR Marker in mm

    ARParam cam_param_;         // Camera Calibration Parameters
    int patt_id_;               // AR Marker Pattern
    char pattern_filename_[FILENAME_MAX];
    bool reverse_transform_;    // Reverse direction of transform marker -> cam

    double marker_center_[2];   // Physical Center of the Marker
    double marker_trans_[3][4]; // Marker Transform

    // **** for visualisation in rviz
    ros::Publisher rvizMarkerPub_;
    visualization_msgs::Marker rvizMarker_;
    
    int contF;
    bool getCamInfo_;
    CvSize sz_;
    IplImage *capture_;
  };                            // end class ARSinglePublisher
}                               // end namespace ar_pose

#endif
