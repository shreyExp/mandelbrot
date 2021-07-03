#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;

void cumilate(int* histogram, int* &cumilative, int size);
void makeHistZero(int* histogram, int size);
void scaleCumilative(int* cumilative, int hueSize);
void printHist(int *cumilative, int huesize);

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
