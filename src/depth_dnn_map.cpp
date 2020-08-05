//Header files
#include "ros/ros.h"
#include <ros/console.h>
#include "stdlib.h"
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <geometry_msgs/Point.h>
#include <mypkg/BoundingBoxes.h>
#include <mypkg/BoundingBox.h>
#include <mypkg/list.h>
#include <mypkg/list_data.h>
#include <std_msgs/Float64.h>
#include <sensor_msgs/image_encodings.h>
#include <sensor_msgs/Image.h>
#include <image_transport/image_transport.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include "opencv2/core/version.hpp"
#include <math.h> 
#include <mypkg/observations.h>
#include <mypkg/observations_list.h>

//Namespaces
using namespace message_filters;

//Variables from message files and ROS
mypkg::observations_list observations;
mypkg::observations data1;
mypkg::list list;
mypkg::list_data data2;
geometry_msgs::Point debug;

//Simple variables
float lxmin,lxmax,lymin,lymax;
int cont1=0,cont2=0,center_x,center_y,size;

//Camera params
double T=0.12000, fx= 350.386, fy= 350.386, cx= 348.071, cy= 191.094;
int resx=672, resy=376; //resolution of ZED

//double T=0.12000, fx= 700.772, fy= 700.772, cx= 665.141, cy= 367.188;
//int resx=1280, resy=720; //resolution of ZED
//cv::Mat depth_img;

void callback1(const sensor_msgs::ImageConstPtr& msg)
{  
/*   depth_img= cv_bridge::toCvShare(msg,sensor_msgs::image_encodings::TYPE_32FC1)->image;
  cv::Size s = depth_img.size();
  int rows = s.height;
  int cols = s.width;
  ROS_INFO_STREAM("rows"<< rows);
  ROS_INFO_STREAM("cols"<< cols); */

  if (cont1==1){
    observations.data.clear(); // if new detections occur, make sure previous are deleted
    float* disparity = (float*)(&msg->data[0]);
    for (int j = 0; j < size; ++j){
      int centerIdx = list.data[j].center_x + msg->width * list.data[j].center_y;
      float z = T * fx / disparity[centerIdx];
      if (z!= INFINITY){
        float x = ((list.data[j].center_x- cx) * z) / (fx);
        float y = ((list.data[j].center_y- cy) * z) / (fy);
        float range = sqrt(x*x + y*y + z*z);
        float bearing = atan(x/z);
        data1.label= list.data[j].label;
        data1.x= x;
        data1.y= z;
        //ROS_INFO_STREAM(z);
         if (list.data[j].label=="target"){
          debug.x= x;
          debug.y= z;
        }
        //ROS_INFO_STREAM("target_x: "<< debug.y); 
        //ROS_INFO_STREAM("target_y: "<< debug.x);
        observations.data.push_back(data1);
        cont1=0;
        cont2=1;
        }
      }
    }
  }      

void callback2(const mypkg::BoundingBoxes::ConstPtr& objects)
{
  cont1=1;
  size=objects-> bounding_boxes.size();
  list.data.clear();// if new detections occur, make sure previous are deleted
  for (int j = 0; j < objects-> bounding_boxes.size(); ++j){
      lxmin= objects-> bounding_boxes[j].xmin;
      lxmax= objects-> bounding_boxes[j].xmax;
      lymin= objects-> bounding_boxes[j].ymin;
      lymax= objects-> bounding_boxes[j].ymax;
      data2.label= objects-> bounding_boxes[j].Class;
      data2.center_x= (((lxmax-lxmin)/2+lxmin)*513/resx); //numbers are resolution of depth image from redtail
      data2.center_y= (((lymax-lymin)/2+lymin)*257/resy);
      list.data.push_back(data2);
  }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "image_converter");
  //cv::namedWindow(OPENCV_WINDOW);
  ros::NodeHandle n;
  ros::Publisher pub = n.advertise<mypkg::observations_list>("/observations", 1000); 
  ros::Publisher pub2 = n.advertise<geometry_msgs::Point>("/debug", 1000); 

  image_transport::ImageTransport it_(n);
  image_transport::Subscriber sub1 = it_.subscribe("/stereo_dnn_ros/network/output", 10, callback1);
  ros::Subscriber sub2 = n.subscribe("/left/darknet_ros/bounding_boxes", 1000, callback2);

  ros::Rate loop_rate(10);
  while (ros::ok())
  {
    ros::spinOnce();
    if (cont2==1){
      pub.publish(observations);
      pub2.publish(debug);
      cont2=0;
    }
    loop_rate.sleep();
  }
  //cv::destroyWindow(OPENCV_WINDOW);
  return 0;
}


