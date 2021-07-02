#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <json/json.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/videoio.hpp>

using namespace std;
using namespace cv;
using namespace Json;


int mandelbrot(complex<double> c, int max);
void cumilate(int* histogram, int* &cumilative, int size);
void makeHistZero(int* histogram, int size);
void scaleCumilative(int* cumilative, int hueSize);
void printHist(int *cumilative, int huesize);
void mandelbrotImage(Mat M, double remin, double remax, double immin, double immax);
void findRegionOfIntrest(Mat M, double& remin, double& remax, double& immin, double& immmax);
bool isRegionOfIntrest(Mat subset);

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
		imshow("hello", M);
		waitKey(100);
		vid<<M;
		cout<<"Out"<<endl;
		//cout<<"dleft "<<dleft<<" dright "<<dright<<" dup "<<dup<<" ddown "<<ddown<<endl;
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

void findRegionOfIntrest(Mat M, double& dleft, double& dright, double& dup, double& ddown){
	uchar *p;
	int noOfBlocks = 100;
	int lengthOfBlockHz = M.cols/noOfBlocks;
	int lengthOfBlockVr = M.rows/noOfBlocks;
	int noOfMoves = 0;
	bool isReg = 0;
	int startIndx = noOfBlocks/2;
	int blkIndxHz = startIndx;
	int blkIndxVr = startIndx;
	int noOfMovesSpec = 0;
	int move = 1;
	int dir = 0;
	bool firstTime = 1;
	bool doit = 1;
	while(doit){
		if((noOfMoves+1)%2)
			noOfMovesSpec = (noOfMoves + 1)/2 + 1;
		else
			noOfMovesSpec = (noOfMoves + 1)/2;
		noOfMoves++;
		//Now design a function which takes direction and moves
		//Gotta let the whole damn thing work for the first time
		//cout<<"noOfMovesSpec: "<<noOfMovesSpec<<endl;
		while(noOfMovesSpec){
			if(!firstTime){
				if(dir == 0)
					blkIndxHz++;
				else if(dir == 1)
					blkIndxVr++;
				else if(dir == 2)
					blkIndxHz--;
				else if(dir == 3)
					blkIndxVr--;
			}
			firstTime = 0;
			//cout<<noOfMoves<<endl;
			//cout<<"blkIndxHz: "<<blkIndxHz<<"blkIndxVr "<<blkIndxVr<<endl;
			Mat subset(M, Rect(blkIndxHz*lengthOfBlockHz, blkIndxVr*lengthOfBlockVr, lengthOfBlockHz, lengthOfBlockVr));
			isReg = isRegionOfIntrest(subset);
			if(isReg){
				//Here goes the code to calculate the percentage decrement needed at borders.
				dleft = (double)blkIndxHz/(double)noOfBlocks;
				dright = 1 - dleft;
				dup = (double)blkIndxVr/(double)noOfBlocks;
				ddown = 1 - dup;
				dleft = 0.1*dleft;
				dright = 0.1*dright;
				dup = 0.1*dup;
				ddown = 0.1*ddown;
				doit = 0;
				break;
			}
			noOfMovesSpec--;
		}

		dir++;
		if(dir%4 == 0){
			dir = 0;
		}
	}
}
bool isRegionOfIntrest(Mat subset){
	int cols = subset.cols*subset.channels();	
	int noOfSetElems = 0;
	bool value = 0;
	uchar* p;
	for(int i = 0; i < subset.rows; i++){
		p = subset.ptr<uchar>(i);
		for(int j = 0; j < cols; j += 3){
			if(p[j+2] == 255){
				noOfSetElems++;
			}
		}
	}
	double proportion = (double)noOfSetElems/(double)(subset.cols*subset.rows);
	//cout<<"proportion: "<<proportion<<endl;
	if(proportion > 0.30 && proportion < 0.70)
		value = 1;
	else
		value = 0;
	return value;
}

void mandelbrotImage(Mat M, double remin, double remax, double immin, double immax){
	M = Scalar(1, 255, 255);
	int maxIterBlue = 255;
	complex<double> c;
	
	double rewidth = (remax - remin);
	double imwidth = (immax - immin);
	uchar *p;
	int cols = M.cols*M.channels();
	int mandel;
	const int hueSize = 180;
	int histogram[hueSize];

	makeHistZero(histogram, hueSize);
	string filenameImage;
	string prefix;
	Value outValue; 
	int count = 0;
	makeHistZero(histogram, hueSize);
	for(int i = 0; i < M.rows; i++){
		p = M.ptr<uchar>(i);
		for(int j = 0; j < cols; j += 3){
			c = complex<double>(remin + (double)j*(remax-remin)/(3*M.cols), 
							immin + (double)i*(immax-immin)/M.rows);

			if((mandel = mandelbrot(c, maxIterBlue)) == maxIterBlue){
				p[j+2] = 0;
			}
			else{
				p[j] = mandel*hueSize/maxIterBlue;
				histogram[p[j]]++;
			}
		}
	}
	int* cumilative;
	cumilate(histogram, cumilative, hueSize);
	if(cumilative[hueSize-1] == 0){
		printHist(cumilative, hueSize);
		printHist(histogram, hueSize);
	}

	scaleCumilative(cumilative, hueSize);
	for(int i = 0; i < M.rows; i++){
		p = M.ptr<uchar>(i);
		for(int j = 0; j < cols; j += 3){
			if(p[j+2] != 0)
				p[j] = cumilative[p[j]];
		}
	}
	cvtColor(M, M, COLOR_HSV2BGR);
}
/**
 * Calculates the number of iterations it takes for z to be ejected of the magnitude 2
 * The maximum number of iterations which are to be calculated is given by max
 **/
int mandelbrot(complex<double> c, int max){
	complex<double> z = (0.0,0.0);
	int iter = 0;
	while(abs(z) <= 2 && iter < max){
		z = z*z + c;
		iter++;
	}
	return iter;
}
/**
 * Calculates the cumilative histogram
 **/
void cumilate(int* histogram, int* &cumilative, int size){
	cumilative = new int[size];
	cumilative[0] = histogram[0];
	for(int i = 1; i < size; i++){
		cumilative[i] = cumilative[i-1] + histogram[i];
	}
}
/**
 * Sets the array values to zero just to be safe.
 **/
void makeHistZero(int* histogram, int size){
	histogram[0] = 1;
	for(int i = 1; i < size; i++){
		histogram[i] = 0;
	}
}
/**
 * Scales the cumilative histogram to be exactly between 0-179
 **/
void scaleCumilative(int* cumilative, int hueSize){
	for(int i = 0; i < hueSize; i++){
		cumilative[i] = cumilative[i]*179/cumilative[hueSize-1];
	}
}
void printHist(int *cumilative, int huesize){
	for(int i = 0; i < huesize; i++){
		cout<<cumilative[i]<<"\t";
	}
	cout<<endl;
}
