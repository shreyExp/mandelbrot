#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <json/json.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;
using namespace Json;


int mandelbrot(complex<double> c, int max);
void cumilate(int* histogram, int* &cumilative, int size);
void makeHistZero(int* histogram, int size);
void scaleCumilative(int* cumilative, int hueSize);
void printHist(int *cumilative, int huesize);
void mandelbrotImage(Mat M, double remin, double remax, double immin, double immax);

int main(int argc, char** argv){
	char* filename; 
	double gloremin;
	double gloremax;
	double gloimmin;
	double gloimmax;
	int divisions;
	string directory;
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
	}else{
		gloremin = -2;
		gloremax = 1;
		gloimmin = -1;
		gloimmax = 1;
	}

	Mat M(500, 500, CV_8UC3, Scalar(1,255,255));
	int steps = 1000;
	double remin = gloremin;
	double remax = gloremax;
	double immin = gloimmin;
	double immax = gloimmax;
	double rerange;
	double imrange;
	double reductionRate = 0.01;
	
	for(int i = 0; i < steps; i++){
		mandelbrotImage(M, remin, remax, immin, immax);
		imshow("hello", M);
		waitKey(100);
		rerange = remax - remin;
		imrange = immax - immin;
		remin = reductionRate * rerange + remin;
		remax = remax - reductionRate * rerange;
		immin = reductionRate * imrange + immin;
		immax = immax - reductionRate * imrange;
	}
}

void mandelbrotImage(Mat M, double remin, double remax, double immin, double immax){
	M = Scalar(1, 255, 255);
	int maxIterBlue = 80;
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
