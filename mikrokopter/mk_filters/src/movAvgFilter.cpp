#include <ros/ros.h>
#include <iostream>
#include <stdio.h>
#include <vector> 
#include <std_msgs/Float32.h>
#include <vector>
#include <string>
#include <numeric>

using namespace std;
using namespace ros;

int windsize = 5; //size of filter window 
vector<float> window(windsize);
Publisher output_pub;

char* nameNewTopic(char* topic){
	char* newTopic = (char *)malloc(512);
	strtok(topic,"/");
	char* part0 =(char*) "mikrokopter";
	char* part1 = strtok(NULL,"/");
	char* part2 = strtok(NULL,"/");
	part1 = "filter";
	strcpy(newTopic, "/");
	strcat(newTopic, part0);
	strcat(newTopic, "/");
	strcat(newTopic, part1);
	strcat(newTopic, "/");
	strcat(newTopic, part2);
	return newTopic;
}

void callback(const std_msgs::Float32& val){
	std_msgs::Float32 buffer;
	if (val.data!=window.back()){
		window.push_back(val.data);
		if (window.size()>=windsize){
			window.erase(window.begin());
		}
//		for (int i = 0; i<window.size();i++){
//			cout<<window[i]<<",";
//		}
		float sum_of_elems = accumulate(window.begin(),window.end(),0);
		buffer.data = sum_of_elems/windsize;
	}else{
		buffer.data = val.data;
	}	
//	cout<<endl;
//	cout<<val.data<<":"<<buffer.data<<endl;
	output_pub.publish(buffer);
}

int main(int argc, char** argv){
	char* topic=argv[1];
	ROS_INFO("Median Filter");
	init(argc,argv,"medFilter");
	NodeHandle nh;
	cout<<"Filtering topic: "<<topic<<endl;
	Subscriber sub = nh.subscribe(topic, 1, callback);
	char* newTopic = nameNewTopic(topic);
	cout<<"Publishing to a new topic: "<<newTopic<<endl;
	output_pub = nh.advertise<std_msgs::Float32>(newTopic,1000);
	spin();
	return 0;
}
