#include "utility.hpp"

#include <cmath>

Mat log10(const Mat& matrix) {

	Mat ret;

	// Calculo del logaritmo natural
	log(matrix, ret);

	// Criba de los valores no deseados
	ret = max(ret, 0.0);

	// Conversion a logaritmo decimal
	ret = ret / log(10.0);

	return ret;

}

Mat to16U(const Mat& matrix) {

	Mat ret;
	matrix.copyTo(ret);

	double vMin, vMax;

	minMaxLoc(ret, &vMin, &vMax, NULL, NULL);
	ret += abs(vMin);

	minMaxLoc(ret, &vMin, &vMax, NULL, NULL);
	ret /= vMax;
	ret *= 65535;

	ret.convertTo(ret, CV_16UC1);

	return ret;

}
