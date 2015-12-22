#ifndef _HOUGH_UTILITIES_H_
#define _HOUGH_UTILITIES_H_

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "utility.hpp"
using namespace cv;
using namespace std;

/* Constant declaration */
#define MAX_IMAGESIZE   	2048
#define MAX_BRIGHTNESS  	65535 	/* Maximum gray level */
#define MAX_BRIGHTNESS_8  	255 	/* Maximum gray level */
#define GRAYLEVEL       	65535 	/* No. of gray levels */
#define GRAYLEVEL_8			255
#define MAX_FILENAME    	256 	/* Filename length limit */
#define MAX_BUFFERSIZE  	256
#define PORCENTAJE 			5		//Porcentaje de Variacion del Radio del Disco Solar
#define PORCENTAJE_CENTRO	20		//Porcentaje de Variacion maximo de los deplazamientos del centro segun enfoque del disco solar
#define FACTOR_RADIO		0.25
#define RADIO_INICIAL 		960


float xGradient(Mat image, int x, int y);
float yGradient(Mat image, int x, int y);
Mat gradient(const Mat& im_in);
//Mat dilate(Mat& image);
Mat erode(const Mat& im_in);
Mat edge(const Mat& image1, const Mat& image2);
int otsu_th(const Mat& im_in, int xmax,int ymax);
Mat binarizar (const Mat& im_in);
int Min(const int a[]);
Vector<Point> findOnes (const Mat& im_in);
Vector<int> findMaximo (const Mat& im_in, int radio);
Vector<int> findMasVotado (int** matrizMaximos, const int size);
void diferenciaPosiciones(double** centrosMedia);
double** calcularMediaCentros(double** centros, double** centrosRotados);
Vector<Point> Get_Randon_points (const Vector<Point> puntos,int N);
Point Centrodegravedad (const Vector<Point> puntos);


#endif
