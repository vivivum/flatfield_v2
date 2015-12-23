#include "hough.h"
#include "houghUtilities.h"

/**
 *	Calculates Maximum matrix in each turn, changing circle radius value
 *
 *  @param puntos		Conjinto de puntos (obtained from Binarized image ) sobre el circulo para estimar coordenadas del centro y radio
 *  @param Centro 		punto central estimado
 *   @param radio 		Approximate radius of the circumference
 *
 *  @return centro  	punto central de salida
 */

double** hough(const Vector<Point> puntos, Point Centro, int radio ){
	//Variables de la funcion
	//cout << "hough" << endl;
	double r;
	double det;
	int aa,bb;
	float a,b;
	int lmax;
	int dimensionAcumulador = 2000;//ojo que es uno mas

	//Contendran los valores y posiciones de los pixeles max y min de la imagen
	double min,max;
	Point min_loc, max_loc;

	int xmax = 2048;

	float paso=0.1;
	float Xmin,Xmax;
	float Ymin,Ymax;
	Xmin=Centro.x-paso*(dimensionAcumulador/2);// Xmin and Xmax are boundaries of coordinates of initial center
	Xmax=Centro.x+paso*(dimensionAcumulador/2);
	Ymin=Centro.y-paso*(dimensionAcumulador/2);// Ymin and Ymax are boundaries of coordinates of initial center
	Ymax=Centro.y+paso*(dimensionAcumulador/2);

	cout << Xmin << " " << Xmax << " "<< Ymin << " "<< Ymax << endl;

	//Calculamos cuantas pasadas realizaremos teniendo en cuenta el tamaño
	//de la imagen y nos aseguramos que se un numero impar
	lmax = (double)PORCENTAJE/100 * xmax;
	if(lmax%2 == 0)
		lmax++;

	cout <<"LMAX: "  << lmax << endl;

	//Creamos la matriz de Maximos, la cual contendra los mejores candidatos de cada pasada
	double** matrizMaximo = new double*[lmax];
	for (int i = 0; i < lmax; ++i)
		matrizMaximo[i] = new double[4];		// Es de dimension 4 por contener: el valor votacion, posX, posY y el radio

	//Recorremos lmax veces la imagenes cambiando el radio
	for(int l = 1; l <= lmax; l++){
		//Calculamos el radio de esta pasada
		r=radio + (l - ((int)(lmax/2)+1))*FACTOR_RADIO;

		//Dimensiones del acumulador 200x200 pixeles
		//Inicializamos el Acumulador
		Mat acu(dimensionAcumulador,dimensionAcumulador, CV_16UC1, Scalar(0));

//int htr;
		for(unsigned int h=0; h < puntos.size();h++){//get points from vector point puntos
			a=Xmin;
			//while (a<=Xmax) {
			for(a=Xmin; a < Xmax; a+=paso){
				//cout << "A: " << a << endl;
				det = r*r - (puntos[h].x-a)*(puntos[h].x-a);
				if(det > 0){
					b=(puntos[h].y-sqrt(det));
					//cout << "B: " << b << endl;
					//cin >> htr;
					if((b > Ymin) && (b < Ymax)){
						aa=round((a-Xmin)/paso);
						bb=round((b-Ymin)/paso);
						acu.at<ushort>(aa,bb) = acu.at<ushort>(aa,bb)+1;
						//cout << "ACU: " << acu.at<ushort>(aa,bb) << endl;
						//cout << "AA: "  << aa << "\tBB: " << bb << endl << endl;
						//cin >> htr;
					}
				}
				//a=a+paso;
			}
		}

		//cout << "Salio del bucle " << endl;
		//Calculo el maximo y su posicion
		minMaxLoc(acu, &min, &max, &min_loc, &max_loc);

		cout << "Max: " << max << ", " << max_loc << "\tMin: " << min << ", " << min_loc << endl;
		//Introducimos los valores del maximo en la matrizMaximo
		matrizMaximo[l-1][0] = max;
		matrizMaximo[l-1][1] = max_loc.x*paso+Xmin;
		matrizMaximo[l-1][2] = max_loc.y*paso+Ymin;
		matrizMaximo[l-1][3] = r;

		//Liberamos la memoria de Acu
		acu.release();
	}
	cout <<"final hough"<<  endl;
	return matrizMaximo;
}

//double** hough(const Mat& img, int radio){
//	//Variables de la funcion
//	//cout << "hough" << endl;
//	double r;
//	int i,j;
//	double det;
//	int a,b;
//	int lmax;
//
//	//Contendran los valores y posiciones de los pixeles max y min de la imagen
//	double min,max;
//	Point min_loc, max_loc;
//
//	//Obtenemos las dimensiones de la imagen
//	Size s = img.size();
//	int xmax = s.height;
//	int ymax = s.width;
//
//	//Calculamos cuantas pasadas realizaremos teniendo en cuenta el tamaño de la imagen y nos aseguramos que se un numero impar
//	lmax = (double)PORCENTAJE/100 * xmax;
//	if(lmax%2 == 0)
//		lmax++;
//
//	//Creamos la matriz de Maximos, la cual contendra los mejores candidatos de cada pasada
//	double** matrizMaximo = new double*[lmax];
//	for (int i = 0; i < lmax; ++i)
//		matrizMaximo[i] = new double[4];		// Es de dimension 4 por contener el valor, posX, posY y el radio
//
//	//Recorremos lmax veces la imagenes cambiando el radio
//	for(int l = 1; l <= lmax; l++){
//		//Calculamos el radio de esta pasada
//		r=radio + (l - ((int)(lmax/2)+1))*FACTOR_RADIO;
//
//		//Inicializamos el Acumulador
//		Mat acu(xmax,ymax, CV_16UC1, Scalar(0));
//		for (i=0; i < xmax; i++){
//			for (j=0; j < ymax; j++){
//				int valorPixel = (img.at<uchar>(j,i));
//				if(valorPixel == MAX_BRIGHTNESS_8){
//					//int a_min = xmax/2 - radio*PORCENTAJE_CENTRO;
//					//int a_max = xmax/2 + radio*PORCENTAJE_CENTRO;
//					// int rango=a_max-a_min+1;		//representa el número de pixels donde se busca la coordenada Xc del centro del disco
//					//for(a=a_min; a < a_max; a++){
//					for(a=0; a < xmax; a++){
//						det = r*r - (i-a)*(i-a);
//						if(det > 0){
//							b=round(j-sqrt(det));
//							//if((b > a_min) && (b < a_max)){
//							if((b > 0) && (b < ymax)){
//								acu.at<ushort>(a,b) = acu.at<ushort>(a,b)+1;
//							}
//						}
//					}//Final For del A
//				}//Final If
//			}
//			//cout << "i = " << xmax << endl;
//		}
//		//cout << "Salio del bucle " << endl;
//		//Calculo el maximo y su posicion
//		minMaxLoc(acu, &min, &max, &min_loc, &max_loc);
//
//		//Introducimos los valores del maximo en la matrizMaximo
//		matrizMaximo[l-1][0] = max;
//		matrizMaximo[l-1][1] = max_loc.x;
//		matrizMaximo[l-1][2] = max_loc.y;
//		matrizMaximo[l-1][3] = r;
//
//		//Liberamos la memoria de Acu
//		acu.release();
//	}
//	//cout <<"final hough"<<  endl;
//	return matrizMaximo;
//}

/**
 *  Pre-procesamos la imagen para pasarla a la funcion Hough y devolvemos la matriz de maximos
 *
 *  @param image	Imagen original a la que calcular la matriz de maximos
 *  @param radio 	Radio aproximado de la circunferencia
 *
 *  @return double** Devolvemos la matriz de maximos
 */
double** calcularHough(Mat &image, int radio){

	//Comprobamos que hemos cargado bien la imagen
	if (!image.data){
		cout << "ERROR. Imagen no cargada correctamente o no encontrada, por favor revise la ruta." << endl;
		return 0;
	}

	//Convertimos la imagen a 32 bits
	//image.convertTo(image, CV_32SC1);

	//Aplicamos los gradientes para hallar los bordes de la imagen
	image = gradient(image);


	//image.convertTo(image, CV_8UC1);
	double min;
	double max;
	minMaxIdx(image, &min, &max);
	Mat adjMap;
	convertScaleAbs(image, adjMap, 255 / max);

	imwrite("gradEscalado.png", adjMap);

	minMaxIdx(adjMap, &min, &max);
	cout << "Minimo: " << min << "\tMaximo: " << max << endl;


	Mat binarizada = binarizar(adjMap);

	Vector<Point> puntosBlancosTotales = findOnes(binarizada);
	Point centroGravedad = Centrodegravedad (puntosBlancosTotales);
	int N = 100;			//Numero de puntos blancos aleatorios
	cout << "N: " << N << endl;
	Vector<Point> randomBlancos =  Get_Randon_points (puntosBlancosTotales, N);

	//Calculamos hough a partir de la imagen con bordes y esa misma imagen binarizada
	return hough( randomBlancos, centroGravedad, radio);
}

//--------------------------------------------------------------------------------------------------------------------
