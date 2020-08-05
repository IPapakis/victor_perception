//Header files
#include "ros/ros.h"
#include <ros/console.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <nav_msgs/Odometry.h>
#include <eigen3/Eigen/Dense>
#include "geometry_msgs/Vector3.h"
#include "geometry_msgs/Quaternion.h"
#include <mypkg/landmarks.h>
#include <mypkg/landmarks_list.h>
#include <mypkg/observations.h>
#include <mypkg/observations_list.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <nav_msgs/Odometry.h>
#include <visualization_msgs/MarkerArray.h>
#include <tf/transform_broadcaster.h>
#include "geometry_msgs/Point.h"
#include "std_msgs/Float64MultiArray.h"
#include <complex>
#include <string>


double robot_moved= 9;//by robot move
double include_dist= 1;


//Namespaces
using namespace message_filters;
using namespace Eigen;

//Variables from message files and ROS
tf::Quaternion quat;
geometry_msgs::Quaternion msg;
mypkg::landmarks lm_database_data; //all observed landmarks database
mypkg::landmarks_list lm_database_msg;
mypkg::landmarks lm_current_data;  //all currently observed landmarks
mypkg::landmarks_list lm_current_msg;
nav_msgs::Odometry husky; //publish
nav_msgs::Odometry odometry;
visualization_msgs::MarkerArray Markerarr;
visualization_msgs::MarkerArray robot_marker;
std_msgs::Float64MultiArray output;

//Simple variables
bool first_loop= true, first_correction=false, odom_first=false, target_include=false,target_seen=false, target_first_time=true;
double lx,ly,Dist,dx,dy,rl,th_l,th_r,q;
double robot_roll, robot_pitch, robot_yaw, robot_x, robot_y;
double rx, ry, ryaw, previous_robot_x, previous_robot_y, previous_robot_yaw,rxc,ryc,ryawc,ox,oy;
int IDlimit,i,total_landmarks,current_landmark;
double distance_param=10;
double previous_target_x=0, previous_target_y=0,target_x,target_y, target_measure, map_measure;
std::string label;
int target_ID=-1;
double target_ratio;
double previous_measure_slam=0;
double previous_measure_target=0;
double target_covariance;
double motion,final_cov,previous_cov;


//Matrix Variables
//fixed size
Matrix<double, 2, 3> M1;
Matrix<double, 2, 2> M3;
Matrix<double, 2, 1> z;
Matrix<double, 2, 2> I2;
Matrix<double, 3, 3> I3;
Matrix<double, 2, 1> delta; 
Matrix<double, 2, 1> meas; 
Matrix<double, 2, 1> pred; 
Matrix<double, 2, 5> H1;
Matrix<double, 2, 2> sensor_noise;
Matrix<double, 2, 2> sensor_noise_original;
Matrix<double, 2, 2> covariance_target;
Matrix<double, 2, 5> H_temp; 
Matrix<double, 5, 2> K_temp; 
Matrix<double, 5, 5> covariance_temp;

//target tracker
Matrix<double, 2, 1> target; 
Matrix<double, 2, 1> target_output; 
Matrix<double, 2, 2> H_target; 
Matrix<double, 2, 2> K_target; 
Matrix<double, 2, 3> H_target2; 
Matrix<double, 3, 2> K_target2; 
Matrix<double, 3, 3> temp_covariance; 
MatrixXd I5(5,5);



//dynamic size matrices
MatrixXd state(3,1);
MatrixXd covariance(3,3);
MatrixXd covariance_ll(2,2);
MatrixXd covariance_lr(2,3);
MatrixXd H2(2,5); 
MatrixXd M2(2,2); 
MatrixXd K(5,2); 
MatrixXd Fc(5,3);
MatrixXd eigenvalues_m(1,1);
MatrixXd eigenvectors_m(1,1);

//System parameters
Matrix<double, 3, 3> motion_covariance;
double sensor_close= 0.4;

void SLAM();
void noise();


void noise(){
    sensor_noise= sensor_noise_original;
    sensor_noise= sensor_noise_original*pow((-0.000216586927084062*pow(sqrt(q),3)+0.0440127941248540*pow(sqrt(q),2)-0.231574299239173*sqrt(q)+1.14503948437718),2);
    sensor_noise(1,1)=0.2*0.2;
}
//Process with SLAM
void SLAM(){
    state(0,0)=rx;
    state(1,0)=ry;
    state(2,0)=ryaw;

    //first deal with those landmarks of the current list that are old, they exist in database from previous step
    
    for (int j = 0; j < lm_current_msg.data.size(); ++j){     
        if (lm_current_msg.data[j].ID<=total_landmarks){  
            //if (j==0){
            //   map_measure=10000;
            //}
            first_correction=true;
            current_landmark= lm_current_msg.data[j].ID-1;
            dx= state((current_landmark*2+3),0)-state(0,0); 
            dy= state((current_landmark*2+4),0)-state(1,0);
            delta<< dx,
                    dy;
            q=delta.transpose()*delta;
            pred<< sqrt(q),
                    atan(dy/dx)-state(2,0);
            meas<< sqrt(pow((lm_current_msg.data[j].x-state(0,0)),2)+pow((lm_current_msg.data[j].y-state(1,0)),2)),
                   atan((lm_current_msg.data[j].y-state(1,0))/(lm_current_msg.data[j].x-state(0,0)))-state(2,0);
            Fc.setZero();
            Fc.topLeftCorner(3,3)=I3;
            Fc.block(3,current_landmark*2+3,2,2)=I2;
            H1<< -sqrt(q)*dx, -sqrt(q)*dy, 0, sqrt(q)*dx, sqrt(q)*dy,
                    dy,          -dx,     -q,  -dy,         dx;
            H1=H1/q;
            H2=H1*Fc;
            noise();
            //zero correlations
            ROS_INFO_STREAM("covariance \n"<<covariance);

            if (target_include==true){
              previous_cov= covariance.block(3+2*(target_ID-1),3+2*(target_ID-1),2,2).determinant();
            } 
            K=  covariance*H2.transpose()*(H2*covariance*H2.transpose()+sensor_noise).inverse();
            if (K==K){
                state=state+K*(meas-pred);
                M2= MatrixXd::Identity(covariance.rows(),covariance.cols());
                covariance= (M2-K*H2)*covariance;
            }else{
            }
            if (target_include==true){
              final_cov= previous_cov/(covariance.block(3+2*(target_ID-1),3+2*(target_ID-1),2,2).determinant());
            } 
        }
    }    


    //now deal with the new landmarks in the current list
    for (int j = 0; j < lm_current_msg.data.size(); ++j){  
        if (lm_current_msg.data[j].ID>total_landmarks){  
            state.conservativeResize(3+2*lm_current_msg.data[j].ID,1); 
            covariance.conservativeResize(3+2*lm_current_msg.data[j].ID,3+2*lm_current_msg.data[j].ID); 
            Fc.resize(5,3+2*lm_current_msg.data[j].ID); 
            H2.resize(2,3+2*lm_current_msg.data[j].ID); 
            K.resize(3+2*lm_current_msg.data[j].ID,2);
            covariance_lr.resize(2,(lm_current_msg.data[j].ID)*2+1);
            if (lm_current_msg.data[j].ID==target_ID){
                current_landmark= lm_current_msg.data[j].ID-1;
                state((current_landmark*2+3),0)= lm_current_msg.data[j].x;//target(0,0);
                state((current_landmark*2+4),0)= lm_current_msg.data[j].y;//target(1,0);
                double dx=state((current_landmark*2+3),0)-state(0,0);
                double dy=state((current_landmark*2+4),0)-state(1,0);
                delta<< dx,
                        dy;
                q=delta.transpose()*delta;
                rl= sqrt(pow(dx,2)+pow(dy,2));
                th_l= atan(dy/dx);
                th_r= state(2,0);
                M1 << 1.0, 0.0, -rl*sin(th_r+th_l),
                        0.0, 1.0, rl*cos(th_r+th_l); 

                M2<< cos(th_r+th_l), -rl*sin(th_r+th_l),
                    sin(th_r+th_l), rl*cos(th_r+th_l);

                noise();
                covariance_ll(0,0)= covariance_target(0,0);
                covariance_ll(1,0)= covariance_target(1,0);
                covariance_ll(0,1)= covariance_target(0,1);
                covariance_ll(1,1)= covariance_target(1,1);
                covariance.bottomRightCorner(2,2)= covariance_ll;
                covariance.block(3+2*(lm_current_msg.data[j].ID-1),0,2,3)= covariance_temp.block(3,0,2,3);
                covariance.block(0,3+2*(lm_current_msg.data[j].ID-1), 3, 2)= covariance_temp.block(3,0,2,3).transpose();
                covariance.block(3,3+(target_ID-1)*2,(lm_current_msg.data[j].ID)*2+1-3,2)= (M1*covariance.block(0,3,3,(lm_current_msg.data[j].ID)*2+1-3)).transpose();
                covariance.block(3+(target_ID-1)*2,3,2,(lm_current_msg.data[j].ID)*2+1-3)= (M1*covariance.block(0,3,3,(lm_current_msg.data[j].ID)*2+1-3));
                //ROS_INFO_STREAM("targetID"<<target_ID);
            }else{
                current_landmark= lm_current_msg.data[j].ID-1;
                state((current_landmark*2+3),0)= lm_current_msg.data[j].x;
                state((current_landmark*2+4),0)= lm_current_msg.data[j].y;
                double dx=lm_current_msg.data[j].x-state(0,0);
                double dy=lm_current_msg.data[j].y-state(1,0);
                delta<< dx,
                        dy;
                q=delta.transpose()*delta;
                rl= sqrt(pow(dx,2)+pow(dy,2));
                th_l= atan(dy/dx);
                th_r= state(2,0);
                M1 << 1.0, 0.0, -rl*sin(th_r+th_l),
                        0.0, 1.0, rl*cos(th_r+th_l); 

                M2<< cos(th_r+th_l), -rl*sin(th_r+th_l),
                    sin(th_r+th_l), rl*cos(th_r+th_l);

                noise();
                covariance_ll=M1*covariance.block(0,0,3,3)*M1.transpose()+M2*sensor_noise*M2.transpose();
                covariance.bottomRightCorner(2,2)= covariance_ll;
                covariance_lr=M1*covariance.block(0,0,3,(lm_current_msg.data[j].ID)*2+1);
                covariance.block(3+2*(lm_current_msg.data[j].ID-1),0,2,3+2*lm_current_msg.data[j].ID-2)=covariance_lr;
                covariance.block(0,3+2*(lm_current_msg.data[j].ID-1), 3+2*lm_current_msg.data[j].ID-2, 2)=covariance_lr.transpose();
            }
        }
    }
     //pass the matrix values to the variables 
    //robot corrected position after landmark observation. changes only when landmark is observed
    rxc= state(0,0);
    ryc= state(1,0);
    ryawc= state(2,0);
    //robot corrected position after landmark observation. changes when landmark is observed and when odometry is inputted.
    rx= state(0,0);
    ry= state(1,0);
    ryaw= state(2,0);
    //update the landmarks database with corrected positions
    for (int j=0; j< total_landmarks; ++j){   
        lm_database_msg.data[j].x=state(j*2+3,0);
        lm_database_msg.data[j].y=state(j*2+4,0); 
    }
}

void callback1(const mypkg::observations_list::ConstPtr& observations){ // observation input and association step

    if (odom_first==true){
        total_landmarks= lm_database_msg.data.size();
        previous_robot_x= robot_x; // robot states at the time of every landmark observation
        previous_robot_y= robot_y;
        previous_robot_yaw= robot_yaw;
        
        for (i = 0; i < observations->data.size(); ++i)
        {
            ox= observations->data[i].y;
            oy= -observations->data[i].x;
            label= observations->data[i].label;

            //target tracking
            if (label=="target" && target_include==false && (robot_x>robot_moved)){ //&& (robot_x>14)
                if (target_seen==false){
                    target(0,0)= (cos(ryaw)*ox-sin(ryaw)*oy) + rx;
                    target(1,0)= (sin(ryaw)*ox+cos(ryaw)*oy) + ry;
                    rl= sqrt(pow(ox,2)+pow(oy,2));
                    th_l= atan(oy/ox);
                    th_r= ryaw;
                    M1 << 1.0, 0.0, -rl*sin(th_r+th_l),
                            0.0, 1.0, rl*cos(th_r+th_l);
                    M2<< cos(th_r+th_l), -rl*sin(th_r+th_l),
                        sin(th_r+th_l), rl*cos(th_r+th_l);
                    dx= target(0,0)-rx;
                    dy= target(1,0)-ry;
                    delta<< dx,
                            dy;
                    q=delta.transpose()*delta;    
                    noise();
                    covariance_target= M1*covariance.block(0,0,3,3)*M1.transpose()+M2*sensor_noise*M2.transpose();
                }else{
                    dx= target(0,0)-rx;
                    dy= target(1,0)-ry;
                    delta<< dx,
                            dy;
                    q=delta.transpose()*delta;
                    pred<< sqrt(q),
                            atan(dy/dx)-ryaw;
                    meas<< sqrt(pow(oy,2)+pow(ox,2)),
                            atan(oy/ox);

                    noise();                       
                    covariance_temp.bottomRightCorner(2,2)= covariance_target;
                    H_temp<< -sqrt(q)*dx, -sqrt(q)*dy, 0, sqrt(q)*dx, sqrt(q)*dy,
                                dy,          -dx,     -q,  -dy,         dx;
                    H_temp= H_temp/q;

                    covariance_temp.block(0,0,3,3)= covariance.block(0,0,3,3); 
                    K_temp= (covariance_temp*H_temp.transpose()*(H_temp*covariance_temp*H_temp.transpose()+sensor_noise).inverse());
                    //K_target= (covariance_target)*H_target.transpose()*(H_target*(covariance_target)*H_target.transpose()+sensor_noise).inverse();

                    if (K_temp==K_temp){
                        target=target+K_temp.block(3,0,2,2)*(meas-pred);
                        I5= MatrixXd::Identity(5, 5);
                        previous_cov= covariance_temp.block(3,3,2,2).determinant();
                        
                        covariance_temp= (I5-(K_temp*H_temp))*covariance_temp; 
                        covariance_target= covariance_temp.block(3,3,2,2);
                        final_cov= previous_cov/(covariance_target.determinant());
                        //(I2-K_target*H_target)*covariance_target;
                    }
                    //ROS_INFO_STREAM("tracking \n"<<covariance_temp);
                }
                target_seen=true;
            }
    

            ///map management and slam
            if ( (sqrt(pow(ox,2)+pow(oy,2))<4.0 && sqrt(pow(ox,2)+pow(oy,2))!=0 && label=="map_landmark") || (label=="target" && target_include==true)){  //

                if (first_loop== true) //Only for the first landmark
                {
                    first_loop=false;
                    lm_database_data.x= (cos(ryaw)*ox-sin(ryaw)*oy) + rx;
                    lm_database_data.y= (sin(ryaw)*ox+cos(ryaw)*oy) + ry;
                    lm_database_data.ID=1;
                    lm_database_msg.data.push_back(lm_database_data);
                    lm_current_data.x= lm_database_data.x;
                    lm_current_data.y= lm_database_data.y;
                    lm_current_data.ID=1;
                    lm_current_msg.data.push_back(lm_current_data);
                }else{
                    IDlimit= lm_database_msg.data.size();
                    lm_current_data.x= (cos(ryaw)*ox-sin(ryaw)*oy) + rx;
                    lm_current_data.y= (sin(ryaw)*ox+cos(ryaw)*oy) + ry;
                    if (label=="target" && target_include==true){
                        if (target_first_time==true){
                            target_ID= IDlimit+1;
                            lm_database_data.x= lm_current_data.x;
                            lm_database_data.y= lm_current_data.y;
                            lm_database_data.ID= target_ID;
                            lm_database_msg.data.push_back(lm_database_data);
                            lm_current_data.ID=target_ID;
                            lm_current_msg.data.push_back(lm_current_data);
                            target_first_time=false;
                        }else{
                            lm_database_msg.data[target_ID-1].x= lm_current_data.x;
                            lm_database_msg.data[target_ID-1].y= lm_current_data.y;
                            lm_current_data.ID= target_ID;  
                            lm_current_msg.data.push_back(lm_current_data);  
                        }
                    }else if(label!="target"){
                        for (int j = 0; j < IDlimit; ++j)
                        {
                            Dist=sqrt(pow(lm_database_msg.data[j].x-lm_current_data.x,2)+pow(lm_database_msg.data[j].y-lm_current_data.y,2));
                            if (Dist<1.5 && (j+1)!=target_ID){ //true if the point is close to the one of the others put it in same position
                                lm_database_msg.data[j].x= lm_current_data.x;
                                lm_database_msg.data[j].y= lm_current_data.y;
                                lm_current_data.ID= j+1;
                                lm_current_msg.data.push_back(lm_current_data); 
                                break;
                            }else if(Dist>1.5 && j==IDlimit-1){ //if not then add it at the end
                                int last=lm_database_msg.data.size();
                                lm_database_data.x= lm_current_data.x;
                                lm_database_data.y= lm_current_data.y;
                                lm_database_data.ID=last+1;
                                lm_database_msg.data.push_back(lm_database_data);
                                lm_current_data.ID=last+1;
                                lm_current_msg.data.push_back(lm_current_data);  
                                break;         
                            }
                        }
                    }
                }
                
            }
            SLAM();
        }
        
        lm_current_msg.data.clear();
    }
}


void callback2(const nav_msgs::Odometry::ConstPtr& odom){
    
    odom_first=true;
    robot_x= odom->pose.pose.position.x; //only measurements
    robot_y= odom->pose.pose.position.y;
    msg= odom->pose.pose.orientation;
    tf::quaternionMsgToTF(msg, quat);
    tf::Matrix3x3(quat).getRPY(robot_roll, robot_pitch, robot_yaw);
    odometry.pose.pose.position.x= robot_x;
    odometry.pose.pose.position.y= robot_y;
    odometry.pose.pose.orientation= msg;

    if (first_correction== false){
        rx= robot_x; 
        ry= robot_y;
        ryaw= robot_yaw;
        covariance.block(0,0,3,3)= covariance.block(0,0,3,3)+motion_covariance;
    }else{
        rx= (robot_x-previous_robot_x)+rxc; //final robot states
        ry= (robot_y-previous_robot_y)+ryc;
        ryaw= (robot_yaw-previous_robot_yaw)+ryawc;   
        covariance.block(0,0,3,3)= covariance.block(0,0,3,3)+motion_covariance;
        motion= motion+motion_covariance(0,0);
        //ROS_INFO_STREAM("MOTION"<<motion);
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "SLAM");
    ros::NodeHandle n; 
    I2<< 1,0,
        0,1;
    I3<< 1,0,0,
        0,1,0,
        0,0,1; 

    //motion_noise
    motion_covariance<< 0.005, 0, 0,
                        0, 0.05, 0,
                        0, 0, 0.001;

    //motion_covariance<< 0.005, 0, 0,
    //                    0, 0.005, 0,
    //                    0, 0, 0.01;

    sensor_noise_original<< 1, 0,
                            0, 1;

    //initialization
    covariance<< 0, 0, 0,
                 0, 0, 0,
                 0, 0, 0;

    //publishers
    ros::Subscriber sub1 = n.subscribe("/observations", 1, callback1);
    ros::Subscriber sub2 = n.subscribe("/odometry/filtered", 1, callback2);

    //subscribers
    ros::Publisher pub1 = n.advertise<nav_msgs::Odometry>("/odometry", 1000); 
    ros::Publisher pub2 = n.advertise<visualization_msgs::MarkerArray>("/landmarks", 1000);
    ros::Publisher pub5 = n.advertise<visualization_msgs::MarkerArray>("/robot_covariance", 1000);  
    ros::Publisher pub3 = n.advertise<nav_msgs::Odometry>("/husky", 1000); 
    ros::Publisher pub4 = n.advertise<std_msgs::Float64MultiArray>("/output", 1000); 

    husky.header.frame_id= "map";
    odometry.header.frame_id= "map";

    tf::TransformBroadcaster br;
    tf::StampedTransform transform;
    transform.frame_id_ = "map";
    transform.child_frame_id_ = "odom";

    int cnt=0;
    ros::Rate loop_rate(10);
    while (ros::ok())
    {
        ros::spinOnce();
    
        //robot position estimation
        husky.pose.pose.position.x= rx;
        husky.pose.pose.position.y= ry;
        husky.pose.pose.orientation= tf::createQuaternionMsgFromRollPitchYaw(0,0,ryaw);

        //marker topic for robot
        ComplexEigenSolver<MatrixXcd> ces(covariance.block(0,0,3,3));
        MatrixXcd values= ces.eigenvalues();
        MatrixXcd vectors= ces.eigenvectors();
        robot_marker.markers.resize(cnt+1);
        robot_marker.markers[cnt].header.frame_id = "map";
        robot_marker.markers[cnt].header.stamp = ros::Time::now();
        robot_marker.markers[cnt].ns = "robot";
        robot_marker.markers[cnt].id = cnt;
        robot_marker.markers[cnt].action = visualization_msgs::Marker::ADD;
        robot_marker.markers[cnt].type = visualization_msgs::Marker::CYLINDER;
        robot_marker.markers[cnt].pose.position.x = rx;
        robot_marker.markers[cnt].pose.position.y = ry;
        robot_marker.markers[cnt].pose.position.z = 0;
        robot_marker.markers[cnt].scale.x = sqrt(abs(ces.eigenvalues()(0)));
        robot_marker.markers[cnt].scale.y = sqrt(abs(ces.eigenvalues()(1)));
        robot_marker.markers[cnt].scale.z = 0.02;
        robot_marker.markers[cnt].color.a = 0.8;
        robot_marker.markers[cnt].color.r = 0.0;
        robot_marker.markers[cnt].color.g = 0.0;
        robot_marker.markers[cnt].color.b = 1.0;
        cnt=cnt+1;

        //marker topic for map
        int landm_number= lm_database_msg.data.size();
        if (target_seen==true){
            Markerarr.markers.resize(landm_number+1);
        }else{
            Markerarr.markers.resize(landm_number);
        }
        for (int i = 0; i < landm_number; i++)
        {
            ComplexEigenSolver<MatrixXcd> ces(covariance.block(3+(i)*2,3+(i)*2,2,2));
            MatrixXcd values= ces.eigenvalues();
            MatrixXcd vectors= ces.eigenvectors();
            Markerarr.markers[i].header.frame_id = "map";
            Markerarr.markers[i].header.stamp = ros::Time::now();
            Markerarr.markers[i].ns = "points_and_lines";
            Markerarr.markers[i].id = i+1;
            Markerarr.markers[i].action = visualization_msgs::Marker::ADD;
            Markerarr.markers[i].type = visualization_msgs::Marker::CYLINDER;
            Markerarr.markers[i].pose.position.x = state(3+(i)*2,0);    //target will also eventually be included here
            Markerarr.markers[i].pose.position.y = state(4+(i)*2,0);
            Markerarr.markers[i].pose.position.z = 0;
            Markerarr.markers[i].scale.x = sqrt(abs(ces.eigenvalues()(0)));
            Markerarr.markers[i].scale.y = sqrt(abs(ces.eigenvalues()(1)));
            //ROS_INFO_STREAM(values);
            //ROS_INFO_STREAM(abs(ces.eigenvalues()(0)));
            //ROS_INFO_STREAM(abs(ces.eigenvalues()(1)));
            Markerarr.markers[i].scale.z = 0.02;
            Markerarr.markers[i].color.a = 0.8;
            Markerarr.markers[i].color.r = 1.0;
            Markerarr.markers[i].color.g = 0.1;
            Markerarr.markers[i].color.b = 0.0;
        }
        if (target_seen==true){
            int pos= landm_number;
            ComplexEigenSolver<MatrixXcd> ces(covariance_target.block(0,0,2,2));
            MatrixXcd values= ces.eigenvalues();
            MatrixXcd vectors= ces.eigenvectors();
            Markerarr.markers[landm_number].header.frame_id = "map";
            Markerarr.markers[landm_number].header.stamp = ros::Time::now();
            Markerarr.markers[landm_number].ns = "points_and_lines";
            Markerarr.markers[landm_number].id = landm_number+1;
            Markerarr.markers[landm_number].action = visualization_msgs::Marker::ADD;
            Markerarr.markers[landm_number].type = visualization_msgs::Marker::CYLINDER;
            Markerarr.markers[landm_number].pose.position.x = target(0,0); //the tracked always position, if not now tracked then the last tracked
            Markerarr.markers[landm_number].pose.position.y = target(1,0);
            Markerarr.markers[landm_number].pose.position.z = 0;
            Markerarr.markers[landm_number].scale.x = sqrt(abs(ces.eigenvalues()(0)));
            Markerarr.markers[landm_number].scale.y = sqrt(abs(ces.eigenvalues()(1)));
            Markerarr.markers[landm_number].scale.z = 0.02;
            Markerarr.markers[landm_number].color.a = 0.8;
            Markerarr.markers[landm_number].color.r = 0.0;
            Markerarr.markers[landm_number].color.g = 1.0;
            Markerarr.markers[landm_number].color.b = 0.0;
        }

        //tf publisher
        transform.setOrigin(tf::Vector3(-rx+robot_x, -ry+robot_y, 0));
        transform.setRotation(tf::Quaternion(0, 0, ryaw-robot_yaw));
        transform.stamp_ = ros::Time::now();
        br.sendTransform(transform);

        if (target_include==true){
            target_output(0,0)=state(3+2*(target_ID-1),0);
            target_output(1,0)=state(4+2*(target_ID-1),0);
        }else{
            target_output(0,0)=target(0,0);
            target_output(1,0)=target(1,0);
        }
        //ROS_INFO_STREAM("\n target position \n"<<target_output);
        if (target_include==true){
            target_covariance= covariance.block(3+2*(target_ID-1),3+2*(target_ID-1),2,2).determinant();
        }else{
            target_covariance= covariance_target.determinant();
        }
        //ROS_INFO_STREAM("\n target covariance \n"<<covariance_target);

        //output topic
        output.data.clear();
        output.data.push_back(robot_x);// robot estimated distance travelled
        output.data.push_back(sqrt(pow(rx,2)+pow(ry,2)));// robot estimated distance travelled 
        output.data.push_back(target_output(0,0));// target estimated position x
        output.data.push_back(target_output(1,0));// target estimated position y
        output.data.push_back(0.5*log(4973.1011016823*covariance.block(0,0,3,3).determinant()));//entropy of robot
        output.data.push_back(0.5*log(291.3521265216*target_covariance));// entropy of target
        output.data.push_back(0.5*log(final_cov)); //target info gain

  
        //publish odometry
        pub1.publish(odometry); //odometry
        //publish visualization
        pub2.publish(Markerarr); //landmarks and target
        //publish robot estimation and visualization
        pub5.publish(robot_marker); //for robot covariance
        pub3.publish(husky); // for robot position
        //publish any kind of output
        pub4.publish(output);
        double robot_target= sqrt(pow(target_output(0,0),2)+pow(target_output(1,0),2))-sqrt(pow(rx,2)+pow(ry,2));
        if (robot_target< include_dist && robot_target>3){
            target_include=true;
        }
        loop_rate.sleep();
    }
    return 0;
    

}
































































