/* \author Aaron Brown */
// Create simple 3d highway enviroment using PCL
// for exploring self-driving car sensors

//#include "render/render.h"
#include "highway.h"
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	viewer->setBackgroundColor(0, 0, 0);

	// set camera position and angle
	viewer->initCameraParameters();
	float x_pos = 0;
	viewer->setCameraPosition ( x_pos-26, 0, 15.0, x_pos+25, 0, 0, 0, 0, 1);

	Highway highway(viewer);

	//initHighway(viewer);

	int frame_per_sec = 30;
	int sec_interval = 10;
	int frame_count = 0;
	int time_us = 0;

	double egoVelocity = 25;

	while (frame_count < (frame_per_sec*sec_interval))
	{
	    std::cout << "frame_count=" << frame_count << "\n";
	    if (frame_count == 99) {
	        std::cout << "Stopping\n";
	    }
		viewer->removeAllPointClouds();
		viewer->removeAllShapes();

		//stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		highway.stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		viewer->spinOnce(1000/frame_per_sec);
		frame_count++;
		time_us = 1000000*frame_count/frame_per_sec;
		
	}

	std::vector<Eigen::VectorXd>& x_hist = highway.traffic[0].ukf.x_history;
    ofstream x_hist_file;
    x_hist_file.open("/home/thomas/Documents/udacity_exercises/SFND_Unscented_Kalman_Filter/x_hist.txt");
    for (Eigen::VectorXd& x : x_hist) {
        for (int i = 0; i < x.size(); i++) {
            x_hist_file << x(i) << "\t";
        }
        x_hist_file << "\n";
    }
    x_hist_file.close();

}