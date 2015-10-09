#include "houghUtilities.h"

/*
 * Computes the x component of the gradient vector
 * at a given point in a image.
 * returns gradient in the x direction using Sobel Operator
 *
 * @param image 	Contains the original image to calculate Sobel
 * @param x			X position of the pixel that represents the center of Sobel gradien in this step
 * @param y			Y position of the pixel that represents the center of Sobel gradien in this step
 *
 * @return int		Value of Sobel gradient in X axis
 */
int xGradient(Mat image, int x, int y){
    return image.at<int>(y-1, x-1) + 2*image.at<int>(y, x-1) + image.at<int>(y+1, x-1)
    	  -image.at<int>(y-1, x+1) - 2*image.at<int>(y, x+1) - image.at<int>(y+1, x+1);
}
//--------------------------------------------------------------------------------------------------------------------


/*
 * Computes the y component of the gradient vector
 * at a given point in a image
 * returns gradient in the y direction using Sobel Operator
 *
 * @param image 	Contains the original image to calculate Sobel
 * @param x			X position of the pixel that represents the center of Sobel gradien in this step
 * @param y			Y position of the pixel that represents the center of Sobel gradien in this step
 *
 * @return int		Value of Sobel gradient in Y axis
 */

int yGradient(Mat image, int x, int y){
    return image.at<int>(y-1, x-1) + 2*image.at<int>(y-1, x) + image.at<int>(y-1, x+1)
    	  -image.at<int>(y+1, x-1) - 2*image.at<int>(y+1, x) - image.at<int>(y+1, x+1);
}
//--------------------------------------------------------------------------------------------------------------------


/*
 * Calculate gradient of image using sobel kernel
 *
 * @param im_im 	32 bits grayscale Input image
 *
 * @return image	32 bits grayscale Gradient image
 */
Mat gradient(const Mat& im_in){
	Size s = im_in.size();
	int xmax = s.height;
	int ymax = s.width;

	float gx, gy, sum;
	Mat image(xmax,ymax, CV_32SC1, Scalar(0));

	for(int y = 1; y < im_in.rows - 1; y++){
		for(int x = 1; x < im_in.cols - 1; x++){
			gx = xGradient(im_in, x, y);
			gy = yGradient(im_in, x, y);
			//sum = sqrt(gx*gx + gy*gy);
			//sum = abs(gx) + abs(gy);
			sum = gx*gx + gy*gy;
			//valor_Pixel = sum > 255 ? 255:sum;
			//valor_Pixel = sum < 0 ? 0 : sum;
			//image.at<uchar>(y,x) = valor_Pixel;
			image.at<int>(y,x) = sum;
		}
	}
	imwrite( "gradient.tiff", image );
	double min,max;
	Point min_loc, max_loc;
	minMaxLoc(image, &min, &max, &min_loc, &max_loc);
	cout << "Max: " << max << endl;

	imshow("gradient.tiff", (image/4294967296)*255);



	waitKey(0);
	return image;
}
//--------------------------------------------------------------------------------------------------------------------


/*
 * Apply Edge dilation to 8 bits image
 *
 * @param image		Contains the image to which the dilate will be applied
 *
 * @return image 	Contains the image once the dilate is applied
 */
Mat dilate(Mat& image){
	Size s = image.size();
	int xmax = s.height;
	int ymax = s.width;

	//Mark neighbouring pixels of whites to modify them later
    for (int i=0; i < xmax; i++){
        for (int j=0; j < ymax; j++){
            if (image.at<uchar>(i,j) == MAX_BRIGHTNESS_8){
            	if (i>0 && image.at<uchar>(i-1,j)==0) 		image.at<uchar>(i-1,j) = 2;		//left neighbour
				if (j>0 && image.at<uchar>(i,j-1)==0) 		image.at<uchar>(i,j-1) = 2;		//top neighbour
				if (i+1<xmax && image.at<uchar>(i+1,j)==0) 	image.at<uchar>(i+1,j) = 2;		//right neighbour
				if (j+1<ymax && image.at<uchar>(i,j+1)==0)	image.at<uchar>(i,j+1) = 2;		//bottom neighbour
            }
        }
    }

    //Set marked pixels to white
    for (int i=0; i<xmax; i++){
        for (int j=0; j<ymax; j++){
            if (image.at<uchar>(i,j) == 2){
            	image.at<uchar>(i,j) = MAX_BRIGHTNESS_8;
            }
        }
    }

    return image;
}
//--------------------------------------------------------------------------------------------------------------------

/*
 * Apply Edge Erode to 8 bits image
 *
 * @param im_in		Contains the image to which the erode will be applied
 *
 * @return image 	Contains the image once the erode is applied
 */
Mat erode(const Mat& im_in){
	Size s = im_in.size();
	int xmax = s.height;
	int ymax = s.width;

	Mat image = im_in.clone();

    for (int i=0; i < xmax; i++){
        for (int j=0; j < ymax; j++){
            if (image.at<uchar>(i,j) == 0){
            	if (i>0 && image.at<uchar>(i-1,j)==MAX_BRIGHTNESS_8) 		image.at<uchar>(i-1,j) = 2;		//left neighbour
				if (j>0 && image.at<uchar>(i,j-1)==MAX_BRIGHTNESS_8) 		image.at<uchar>(i,j-1) = 2;		//top neighbour
				if (i+1<xmax && image.at<uchar>(i+1,j)==MAX_BRIGHTNESS_8) 	image.at<uchar>(i+1,j) = 2;		//right neighbour
				if (j+1<ymax && image.at<uchar>(i,j+1)==MAX_BRIGHTNESS_8)	image.at<uchar>(i,j+1) = 2;		//bottom neighbour
            }
        }
    }
    for (int i=0; i < xmax; i++){
        for (int j=0; j<ymax; j++){
            if (image.at<uchar>(i,j) == 2){
            	image.at<uchar>(i,j) = 0;
            }
        }
    }

    return image;
}
//--------------------------------------------------------------------------------------------------------------------


/*
 * Edge function compares binarized image with gradien image, throw AND Operator
 * The algorithm checks for white pixels in both images and if true then apply white to final image
 *
 * @param binarized		Contains binarized image
 * @param gradient		Contains gradient image
 *
 * @return edgeImage 	Contains the result image
 */
Mat edge(const Mat& binarized, const Mat& gradient){
	Size s = binarized.size();
	int xmax = s.height;
	int ymax = s.width;

	Mat edgeImage(xmax,ymax, CV_8UC1, Scalar(0));

	for (int i=0; i < xmax; i++){
		for (int j=0; j < ymax; j++){
			//If both pixels are equals and white, then we apply white to the final image (AND Operator)
			if(binarized.at<uchar>(i,j) == gradient.at<uchar>(i,j) && binarized.at<uchar>(i,j) == MAX_BRIGHTNESS_8)
				edgeImage.at<uchar>(i,j) = MAX_BRIGHTNESS_8;
			else
				edgeImage.at<uchar>(i,j) = 0;
		}
	}
/*
	imwrite( "edge.tiff", edgeImage );
	imwrite( "binary.tiff", binarized );
	imwrite( "gradient.tiff", gradient );

	Mat image;
	binarized.convertTo(image, CV_32SC1);
	imwrite( "image.tiff", image );

	imshow("edgeimage", edgeImage);
	waitKey(0);
	imshow("binarizada", binarized);
	waitKey(0);
	imshow("gradiente", gradient);
	waitKey(0);
*/


	return edgeImage;
}
//--------------------------------------------------------------------------------------------------------------------

/**
 *  Calculates threshold for the binarization
 *
 *  @param im_in	Image to which calculate threshold
 *  @param xmax		X size of image
 *  @param ymax		Y size of image
 *
 *  @return int 	Return the value of threshold
 */
int otsu_th(const Mat& im_in, int xmax,int ymax)	{
	double max_sigma;
	int i, x, y; /* Loop variable */
	int threshold; /* threshold for binarization */

	/* Histogram generation */
	int *hist = new int[GRAYLEVEL_8];
	for (i = 0; i < GRAYLEVEL_8; i++) hist[i] = 0;
	for (y = 0; y < ymax; y++)
		for (x = 0; x < xmax; x++) {
			hist[im_in.at<uchar>(x,y)]++;
		}



	/* calculation of probability density */
	int numeroPixel = (xmax * ymax);

	double *prob = new double[GRAYLEVEL_8];
	for ( i = 0; i < GRAYLEVEL_8; i ++ ) {
		prob[i] = (double)hist[i] / numeroPixel;
	}
	delete []hist;

	/* omega & myu generation */
	double *omega = new double[GRAYLEVEL_8]; /* prob of graylevels */
	double *myu = new double[GRAYLEVEL_8];   /* mean value for separation */
	omega[0] = prob[0];
	myu[0] = 0.0;       /* 0.0 times prob[0] equals zero */
	for (i = 1; i < GRAYLEVEL_8; i++) {
		omega[i] = omega[i-1] + prob[i];
		myu[i] = myu[i-1] + i*prob[i];
	}
	delete []prob;

	/* sigma maximization
     sigma stands for inter-class variance
     and determines optimal threshold value */
	threshold = 0;
	max_sigma = 0.0;

	double *sigma = new double[GRAYLEVEL_8]; /* inter-class variance */

	for (i = 0; i < GRAYLEVEL_8-1; i++) {
		if (omega[i] != 0.0 && omega[i] != 1.0)
			sigma[i] = pow(myu[GRAYLEVEL_8-1]*omega[i] - myu[i], 2) /
			(omega[i]*(1.0 - omega[i]));
		else
			sigma[i] = 0.0;
		if (sigma[i] > max_sigma) {
			max_sigma = sigma[i];
			threshold = i;
		}
	}

	//Limpiamos los punteros para liberar memoria

	delete []omega;
	delete []myu;
	delete []sigma;

	return threshold;
}
//--------------------------------------------------------------------------------------------------------------------


/**
 *  Generate a binarized image in black & white
 *
 *  @param im_in Imagen to which apply binarization
 *
 *  @return im_b Return a new binarized image
 */
Mat binarizar (const Mat& im_in){
	Size s = im_in.size();
	int xmax = s.height;
	int ymax = s.width;

	Mat im_b(MAX_IMAGESIZE, MAX_IMAGESIZE, CV_8UC1, Scalar(0));

	int threshold = otsu_th(im_in,xmax,ymax);
	//int threshold = 15000;

	cout << "Threshold: " << threshold << endl;

	/* Binarization output into image */
	for (int y = 0; y < ymax; y++){
		for (int x = 0; x < xmax; x++){
			if (im_in.at<uchar>(y,x) > threshold){
				im_b.at<uchar>(y,x) = MAX_BRIGHTNESS_8;
			}
		}
	}
	return im_b;
}
//--------------------------------------------------------------------------------------------------------------------


/**
 *  Search the smallest value on an array
 *
 *  @param a[]	 	Array to analyze
 *
 *  @return int		Return the smallest found value
 */
int Min(const int a[]){
	int min = a[0];
	for (int i = 1; i < NUMERO_IMAGENES; i++){
		if (a[i] < min){
			min = a[i];
		}
	}
	return min;
}
//--------------------------------------------------------------------------------------------------------------------


/**
 *  Buscamos en una imagen los puntos que tiene valor blanco
 *  y los introducimos en un vector que devolvemos
 *
 *  @param im_in Contiene la imagen a analizar
 *
 *  @return Vector<Point> Devolvemos un vector de puntos que tiene valor blanco
 */
Vector<Point> findOnes (const Mat& im_in){
	Size s = im_in.size();
	int xmax = s.height;
	int ymax = s.width;
	Vector <Point> VectorCoordenadas;

	for (int i=0; i < xmax; i++){
		for (int j=0; j < ymax; j++){
			//Si el valor de la imagen en el pixel X,Y es blanco lo guardamos
			if((im_in.at<ushort>(j,i)) == MAX_BRIGHTNESS){
				VectorCoordenadas.push_back(Point(i, j));
			}
		}
	}

	return VectorCoordenadas;
}
//--------------------------------------------------------------------------------------------------------------------

/**
 *  Buscamos en la matriz Acu de la pasada actual el valor maximo de esta
 *  y lo guardamos con sus coordenadas y radio calculado
 *
 *  @param im_in	Contiene la referencia a la matriz Acu actual
 *  @param radio	Pasamos el radio de la pasada actual para introducirlo en el vector
 *
 *  @return Vector<int> Devolvemos un vector con el valor, posicion y radio del elemento con mayor valor de Acu
 */
Vector<int> findMaximo (const Mat& im_in, int radio){

	Size s = im_in.size();
	int xmax = s.height;
	int ymax = s.width;

	Vector <int> VectorMaximoL;	//Almaceno los datos del maximo local
	int xl, yl;
	ushort valorl = im_in.at<ushort>(0,0);

	for (int i=0; i < xmax; i++){
		for (int j=0; j < ymax; j++){
			//Si el valor es mayor que el guardado actual actualizamos
			if((im_in.at<ushort>(i,j)) < valorl){
				xl = i;
				yl = j;
				valorl = im_in.at<ushort>(i,j);
			}
		}
	}

	VectorMaximoL.push_back((int)valorl);			//Valor Maximo de esta pasada
	VectorMaximoL.push_back(xl);					//Posicion X del valor maximo de esta pasada
	VectorMaximoL.push_back(yl);					//Posicion Y del valor maximo de esta pasada
	VectorMaximoL.push_back(radio);					//Valor del Radio de esta pasada

	return VectorMaximoL;
}
//--------------------------------------------------------------------------------------------------------------------


/**
 *  Buscamos en la matriz de los elementos mas votados de cada pasada
 *  el cual tiene un valor mas votado y devolvemos su fila
 *
 *  @param vec_in 	Contiene la matriz de todos los elementos
 *
 *  @return Vector<int>	Devolvemos el vector del elemento mas votado de todos
 */
Vector<int> findMasVotado (int** matrizMaximos, const int size){
	int maxVal = 0;
	int pos = 0;
	for (int i = 0; i < size; i++){
		if(matrizMaximos[i][0] > maxVal){
			maxVal = matrizMaximos[i][0];
			pos = i;
		}
	}
	Vector <int> elementoMasVotado;
	elementoMasVotado.push_back(matrizMaximos[pos][0]);
	elementoMasVotado.push_back(matrizMaximos[pos][1]);
	elementoMasVotado.push_back(matrizMaximos[pos][2]);
	elementoMasVotado.push_back(matrizMaximos[pos][3]);

	return elementoMasVotado;
}
//--------------------------------------------------------------------------------------------------------------------
/*
 *	Calculamos las diferencias de posiciones con respecto a la primera imagen
 */
void diferenciaPosiciones(double** centrosMedia){
	//Restamos el centro de todas las imagenes al primero, para conocer el desplazamiento con respecto a la primera imagen de la imagen "i"
	for(int i=1; i < NUMERO_IMAGENES; i++){
		centrosMedia[i][0] = centrosMedia[0][0] - centrosMedia[i][0];
		centrosMedia[i][1] = centrosMedia[0][1] - centrosMedia[i][1];
	}
	centrosMedia[0][0] = centrosMedia[0][1] = 0;		//No existe diferencia entre el primero y si mismo
}

//--------------------------------------------------------------------------------------------------------------------


/*
 * Obtenemos la media teniendo en cuenta los dos centros que hemos calculado
 */
double** calcularMediaCentros(double** centros,double** centrosRotados){
	double** centrosMedia = new double*[NUMERO_IMAGENES];
	for (int i = 0; i < NUMERO_IMAGENES; i++)
		centrosMedia[i] = new double[2];

	for (int i = 1; i < NUMERO_IMAGENES; i++){
		centrosMedia[i][0] = (centros[i][0] + centrosRotados[i][1])/2;
		centrosMedia[i][1] = (centros[i][1] + centrosRotados[i][0])/2;
	}

	return centrosMedia;
}

