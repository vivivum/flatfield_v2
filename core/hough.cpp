#include "hough.h"

/**
 *	Calculates Maximum matrix in each turn, changing circle radius value
 *
 *  @param img			Binarized image with its edge calculated
 *  @param radio 		Approximate radius of the circumference
 *
 *  @return matrizMaximo 	Return maximum matrix
 */
double** hough(const Mat& img, int radio){
	//Variables de la funcion
	//cout << "hough" << endl;
	double r;
	int i,j;
	double det;
	int a,b;
	int lmax;

	//Contendran los valores y posiciones de los pixeles max y min de la imagen
	double min,max;
	Point min_loc, max_loc;

	//Obtenemos las dimensiones de la imagen
	Size s = img.size();
	int xmax = s.height;
	int ymax = s.width;

	//Calculamos cuantas pasadas realizaremos teniendo en cuenta el tamaño de la imagen y nos aseguramos que se un numero impar
	lmax = (double)PORCENTAJE/100 * xmax;
	if(lmax%2 == 0)
		lmax++;

	//Creamos la matriz de Maximos, la cual contendra los mejores candidatos de cada pasada
	double** matrizMaximo = new double*[lmax];
	for (int i = 0; i < lmax; ++i)
		matrizMaximo[i] = new double[4];		// Es de dimension 4 por contener el valor, posX, posY y el radio

	//Recorremos lmax veces la imagenes cambiando el radio
	for(int l = 1; l <= lmax; l++){
		//Calculamos el radio de esta pasada
		r=radio + (l - ((int)(lmax/2)+1))*FACTOR_RADIO;

		//Inicializamos el Acumulador
		Mat acu(xmax,ymax, CV_16UC1, Scalar(0));
		for (i=0; i < xmax; i++){
			for (j=0; j < ymax; j++){
				int valorPixel = (img.at<uchar>(j,i));
				if(valorPixel == MAX_BRIGHTNESS_8){
					//int a_min = xmax/2 - xmax*PORCENTAJE_CENTRO;
					//int a_max = xmax/2 + xmax*PORCENTAJE_CENTRO;

					//for(a=a_min; a < a_max; a++){
					for(a=0; a < xmax; a++){
						det = r*r - (i-a)*(i-a);
						if(det > 0){
							b=round(j-sqrt(det));
							//if((b > a_min) && (b < a_max)){
							if((b > 0) && (b < ymax)){
								acu.at<ushort>(a,b) = acu.at<ushort>(a,b)+1;
							}
						}
					}//Final For del A
				}//Final If
			}
			//cout << "i = " << xmax << endl;
		}
		//cout << "Salio del bucle " << endl;
		//Calculo el maximo y su posicion
		minMaxLoc(acu, &min, &max, &min_loc, &max_loc);

		//Introducimos los valores del maximo en la matrizMaximo
		matrizMaximo[l-1][0] = max;
		matrizMaximo[l-1][1] = max_loc.x;
		matrizMaximo[l-1][2] = max_loc.y;
		matrizMaximo[l-1][3] = r;

		//Liberamos la memoria de Acu
		acu.release();
	}
	//cout <<"final hough"<<  endl;
	return matrizMaximo;
}

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
	image.convertTo(image, CV_32SC1);

	//Aplicamos los gradientes para hallar los bordes de la imagen
	image = gradient(image);

	image.convertTo(image, CV_8UC1);

	//Calculamos hough a partir de la imagen con bordes y esa misma imagen binarizada
	return hough( edge( binarizar(image) , image) , radio);
}

//--------------------------------------------------------------------------------------------------------------------
