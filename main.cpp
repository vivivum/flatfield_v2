#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include <opencv2/highgui/highgui.hpp>

#include "core/flatfield.hpp"
#include "core/utility.hpp"
#include "core/hough.h"
#include "core/houghUtilities.h"

using namespace cv;
using namespace std;

#define RMIN 30096.8
#define RMAX 65535.0
#define LOOPS 10
#define STATS true

Mat rotate(Mat src, double angle)
{
    Mat dst;
    Point2f pt(src.cols/2., src.rows/2.);
    Mat r = getRotationMatrix2D(pt, angle, 1.0);
    warpAffine(src, dst, r, Size(src.cols, src.rows));
    return dst;
}

int main(int argc, char** argv) {

	char imageName[] = "./img/im0X.tiff";

	vector<Mat> datacube;

	// Leer imagenes desde fichero, guardandolas en el vector de datos
	for(unsigned int i = 0; i < NUMERO_IMAGENES; i++) {

		imageName[9] = 48 + i;					//El 48 es el 0 en ASCII, con lo que sustituimos la X por el numero de la imagen en cada pasada en la posicion 9
		
		datacube.push_back(imread(imageName, -1));
		//datacube[i].convertTo(datacube[i], CV_64FC1);	//OJO leer un solo canal CV_32FC1!!!!!!!!!!!
		datacube[i].convertTo(datacube[i], CV_32FC1);	//OJO leer un solo canal CV_32FC1!!!!!!!!!!!

	}

	//Creamos la matriz de centros de las imagenes desplazadas [NUMERO_IMAGENES][2 coordenadas]
	double** centros = new double*[NUMERO_IMAGENES];
	for (int i = 0; i < NUMERO_IMAGENES; i++)
		centros[i] = new double[2];

	double** centrosRotados = new double*[NUMERO_IMAGENES];
		for (int i = 0; i < NUMERO_IMAGENES; i++)
			centrosRotados[i] = new double[2];

	cout << "Empezamos a calcular Hough" << endl;
	Vector<double**> matrizResultadosHough;

	//-------------------------------------------------
	//Obtenemos el numero de pasada lmax para Hough dependiento del tamaño de la imagen
	Size s = datacube[0].size();
	int xmax = s.height;
	int lmax = (double)PORCENTAJE/100 * xmax;

	if(lmax%2 == 0)		//Compruebo que sea impar
		lmax++;

	cout << "LMAXXX: " << lmax << endl;
	//-------------------------------------------------

	//****************************************************************************************
	//Con el objetivo de hacer estadistica de los resultados bla bla bla
	ofstream fs("EstadisticasdeCentro.txt");
	for(int i=0; i < NUMERO_IMAGENES; i++){
		matrizResultadosHough.push_back(calcularHough(datacube[i], RADIO_INICIAL));	//Es una matriz de la dimension lmax*4 (votacion, x, y, radio)
		cout << "Pasada: " <<i << endl;

		#ifdef STATS
			for (int j=0; j < lmax; j++){
				for(int h=0; h < 4; h++){
					fs << matrizResultadosHough[i][j][h] << "\t";
				}
				fs << endl;
			}
			fs << endl;
		#endif

		#ifndef STATS
			centros[i][1] = matrizResultadosHough[i][1];					//Cogemos el valor del centro X
			centros[i][0] = matrizResultadosHough[i][2];					//Cogemos el valor del centro Y
			cout << "Imagen Manual " << i << ": " << centros[i][0] << ", " << centros[i][1] << endl;
		#endif

		//---------------------------------------------------------------
	}
	fs.close();

	cout << "Hough Finalizado" << endl;

//
//	//****************************************************************************************
//	//OPERACION ROTADA, VOLVEMOS A HACER LA MISMA OPERACION CON LAS IMAGENES ROTADAS Y EL FLIP
//	Mat fliped;
//	fliped.convertTo(fliped, CV_64F);
//	Vector<double**> matrizResultadosHough2;
//	ofstream fs2("EstadisticasdeCentro_rotados.txt");
//	for(int i=0; i < NUMERO_IMAGENES; i++){
//		flip(rotate(datacube[i], -90),fliped,1);
//		matrizResultadosHough2.push_back(calcularHough(fliped, RADIO_INICIAL)); 	//Es una matriz de la dimension lmax*4 (votacion, x, y, radio)
//
//		#ifdef STATS
//			for (int j=0; j < lmax;j++){
//				for(int h=0; h < 4; h++){
//					fs2 << matrizResultadosHough2[i][j][h] << "\t";
//				}
//				fs2 << endl;
//			}
//			fs2 << endl;
//		#endif
//
//		#ifndef STATS
//			centrosRotados[i][1] = matrizResultadosHough[i][1];					//Cogemos el valor del centro X
//			centrosRotados[i][0] = matrizResultadosHough[i][2];					//Cogemos el valor del centro Y
//			cout << "Imagen Manual " << i << ": " << centrosRotados[i][0] << ", " << centrosRotados[i][1] << endl;
//		#endif
//
//	}
//	fs2.close();
//	//***************************************************
//
//	cout << "Hough ROTADO Finalizado\n\nMostramos las diferencias de los centros" << endl;
//
//	double** centrosMedia = calcularMediaCentros(centros, centrosRotados);
//	diferenciaPosiciones(centrosMedia); //Desplazamientos relativos a la primera desplazada


//	for(int i=0; i < 9; i++){
//		cout << "Imagen " << i << ": " << dispersion[i][0] << ", " << dispersion[i][1] << endl;
//	}

/*
	//Evitamos tener que cargar el hough en los test guardando los resultado por AHORA, BORRAR AL EJECUTAR COMPLETO
	dispersion[0][0] = 0; dispersion[0][1] = 0;
	dispersion[1][0] = 7; dispersion[1][1] = -25;
	dispersion[2][0] = 32; dispersion[2][1] = -41;
	dispersion[3][0] = 63; dispersion[3][1] = -28;
	dispersion[4][0] = 77; dispersion[4][1] = -3;
	dispersion[5][0] = 64; dispersion[5][1] = 26;
	dispersion[6][0] = 35; dispersion[6][1] = 34;
	dispersion[7][0] = 3; dispersion[7][1] = 17;
	*/
	//-----------------------------------------------------------------------------------------------------------------

	/*
	for(int i=0; i < 8; i++){
			cout << "Imagen Manual " << i << ": " << dispersion[i][0] << ", " << dispersion[i][1] << endl;
	}

	// Tabla de desplazamientos - TODO Calcular los desplazamientos de forma automatica
	const int disp[8][2] = {{0, 0}, {8, -26}, {32, -42}, {63, -30}, {79, -3}, {64, 24}, {35, 33}, {3, 18}};

	// Declaracion de la plantilla de pixeles buenos dependiendo del numero de Imagenes, el numero de imagenes maximo es 16, 1 bit por imagen
	Mat tmp(datacube[0].size(), CV_16UC1, 0.0);

	// Calculo de la plantilla de pixeles malos
	int err = getImages(datacube, tmp, RMIN, RMAX);

	if (err)
		return 1;

	// Declaracion de la matriz de conteo de pares de pixeles
	Mat pixCnt(datacube[0].size(), CV_8UC1, 0.0);

	// Calculo del termino constante
	Mat con = getConst(datacube, tmp, pixCnt, dispersion);

#ifdef DEBUG

	namedWindow("Constant term", CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
	imshow("Constant term", con);
	waitKey(0);

#endif

	// Calculo de la ganancia unitaria
	Mat pixCntAux = max(pixCnt, 1.0);
	pixCntAux.convertTo(pixCntAux, CV_64F);
	
	Mat gain = con / pixCntAux;

	// Calculo del flatfield
	Mat flat = iterate(con, gain, tmp, pixCnt, disp, LOOPS);
	flat = to16U(flat);

#ifdef DEBUG

	namedWindow("Flatfield", CV_WINDOW_NORMAL | CV_WINDOW_KEEPRATIO);
	imshow("Flatfield", flat);	
	waitKey(0);

#endif

	imwrite("./ff.png", flat);
	*/

	int g;
	cin >> g;
	return 0;

}
