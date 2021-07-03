#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <json/json.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/videoio.hpp>
#include "mandelbrot.h"

using namespace std;
using namespace cv;
using namespace Json;



int main(int argc, char** argv){
	char* filename; 
	double gloremin;
	double gloremax;
	double gloimmin;
	double gloimmax;
	int divisions;
	string directory;
	string output_filename;
	if(argc > 1){
		ifstream in;
		Value root;
		filename = argv[1];
		in.open(filename);
		in>>root;
		in.close();
		gloremin = root["remin"].asDouble();
		gloremax = root["remax"].asDouble();
		gloimmin = root["immin"].asDouble();
		gloimmax = root["immax"].asDouble();
		divisions = root["divisions"].asDouble();
		directory = root["directory"].asString(); 
		output_filename = root["output_filename"].asString();
	}else{
		gloremin = -2;
		gloremax = 1;
		gloimmin = -1;
		gloimmax = 1;
		output_filename = "output.avi";
	}

	Mat M(500, 500, CV_8UC3, Scalar(1,255,255));
	double remin = gloremin;
	double remax = gloremax;
	double immin = gloimmin;
	double immax = gloimmax;
	double rerange;
	double imrange;
	double reductionRate = 0.01;
	double dleft = 0.01;
	double dright = 0.01;
	double dup = 0.01;
	double ddown = 0.01;
	
	const double fps = 30;
	int steps = 200;
	Size S(500, 500);
	cv::VideoWriter vid(output_filename, cv::VideoWriter::fourcc('P','I','M','1'), fps, S, 1);
	for(int i = 0; i < steps; i++){
		mandelbrotImage(M, remin, remax, immin, immax);
		imshow("Mandelbrot Set!", M);
		waitKey(100);
		vid<<M;
		cout<<"Out"<<endl;
		if(!(i%10)){
			findRegionOfIntrest(M, dleft, dright, dup, ddown);
		}
		rerange = remax - remin;
		imrange = immax - immin;
		remin = dleft*rerange + remin;
		remax = remax - dright*rerange;
		immin = dup * imrange + immin;
		immax = immax - ddown * imrange;
	}
	vid.release();
}
