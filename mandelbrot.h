#include <complex>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "imageutilities.h"

using namespace std;
using namespace cv;

int mandelbrot(complex<double> c, int max);
void mandelbrotImage(Mat M, double remin, double remax, double immin, double immax);
void findRegionOfIntrest(Mat M, double& remin, double& remax, double& immin, double& immmax);
bool isRegionOfIntrest(Mat subset);


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



