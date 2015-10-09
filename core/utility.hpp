#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;

#define NUMERO_IMAGENES		9

Mat log10(const Mat&);

Mat to16U(const Mat&);

#endif
