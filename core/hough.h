#ifndef _HOUGH_H_
#define _HOUGH_H_

#include "houghUtilities.h"

double** hough(const Mat& img, int radio);
double** calcularHough(Mat &img, int radio);

#endif /* _HOUGH_H_ */
