#ifndef FLATFIELD_HPP
#define FLATFIELD_HPP

#include <opencv2/highgui/highgui.hpp>

#include <vector>
#include "utility.hpp"

#define ROWS 512
#define COLS 512
#define REL8TO64 72340172838076673

using namespace cv;
using namespace std;

int getImages(vector<Mat>&, Mat&, const double, const double);

//Mat getConst(vector<Mat>&, const Mat&, Mat&, const int[8][2]);		//Descomentar si queremos volver al const int[][]
Mat getConst(vector<Mat>&, const Mat&, Mat&,int**);

void doIteration(const Mat&, Mat&, const Mat&, const Mat&, const int[8][2]);

Mat iterate(const Mat&, Mat&, const Mat&, const Mat&, const int[8][2], const unsigned int);
Mat binarizar (const Mat&);

//NUEVAS FUNCIONES
void getCon(const Mat& mskiq,const Mat& mskir,const Mat& dataiq,const Mat& datair, Mat& pixCnt, Mat& con, int dx, int dy);
Mat getConstNueva(vector<Mat>& data, const Mat& tmp, Mat& pixCnt, int **disp);
void getGainTmp(const Mat& mskiq,const Mat& mskir, Mat& gainTmp, Mat& gain, int dx, int dy);
void calculateStats(const Mat& pixCnt, Mat& gainTmp, Mat& gain);
void doIterationNueva(const Mat& con, Mat& gain, const Mat& tmp, const Mat& pixCnt, const int disp[8][2]);

#endif
