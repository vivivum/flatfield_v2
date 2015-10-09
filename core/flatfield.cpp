#include "flatfield.hpp"
#include "utility.hpp"

#include <algorithm>
#include <cmath>
#include <math.h>
#include <iostream>
/*
 *  8-bit unsigned integer (uchar)   				CV_8U
 *	8-bit signed integer (schar)					CV_8S
 *	16-bit unsigned integer (ushort)				CV_16U
 *	16-bit signed integer (short)					CV_16S
 *	32-bit signed integer (int)						CV_32S
 *	32-bit floating-point number (float)			CV_32F
 *	64-bit floating-point number (double)			CV_64F
 *
 *
 *	CV_8UC1 (single channel array with 8 bit unsigned integers)
 *	CV_8UC2 (2 channel array with 8 bit unsigned integers)
 *	CV_8UC3 (3 channel array with 8 bit unsigned integers)
 *	CV_8UC4 (4 channel array with 8 bit unsigned integers)
 *	CV_8UC(n) (n channel array with 8 bit unsigned integers (n can be from 1 to 512) )
 *
 */


/**
 *  Filtramos las imagenes de entrada y obtenemos la mascara necesaria para calcular el flatfield
 *
 *  @param data	     N Imagenes de entrada
 *  @param tmp 	     Mascara inicial de 16 bits por pixel inicializada a 0
 *	@param rMin		 Umbral minimo de filtrado
 *	@param rMax		 Umbral maximo de filtrado
 *
 */
int getImages(vector<Mat>& data, \
              Mat& tmp, \
              const double rMin, \
              const double rMax) {

	// Comprobar el numero de imagenes
	unsigned int nImages = data.size();

	if(nImages > NUMERO_IMAGENES) {
		cerr << "Error - El numero de imagenes deben de ser " << NUMERO_IMAGENES << "..." << endl;
		return 1;
	}

	//Obtenemos las dimensiones de una imagen cualquiera (Todas deben ser del mismo tamaño)
	Size s = data[0].size();
	int xmax = s.height;
	int ymax = s.width;

	// Calcular la plantilla de pixeles buenos solo sobre las imagenes desplazadas
	for(unsigned int i = 1; i < NUMERO_IMAGENES; i++) {

		Mat msk(xmax,ymax,CV_16U,Scalar(0));

		for(int x=0; x < xmax; x++){
			for(int y=0; y < ymax; y++){
				if(data[i].at<double>(y,x) >= rMin && data[i].at<double>(y,x) <= rMax){
					msk.at<ushort>(y,x) = 1;
				}
			}
		}

		cout << "Imagen " << i + 1 << ": " << (unsigned int)(sum(msk).val[0]) << " pixeles buenos..." << endl;

		// Acumula la mascara en la plantilla de pixeles buenos
		tmp = tmp | (msk * (1 << i));

		// Emplear la mascara para cribar los pixeles malos
		msk.convertTo(msk, data[i].type());
		//Marcamos los pixeles malos sobre los datos
		data[i] = data[i].mul(msk);
	}
	return 0;
}



Mat getConstOlder(vector<Mat>& data, \
             const Mat& tmp, \
             Mat& pixCnt, \
             int **disp) {

	Mat con(data[0].size(), CV_64F);
	vector<Mat> dat;

	unsigned int loopCnt = 0;

	// Calculo del logaritmo comun (base 10) de la imagen
	dat.push_back(log10(data[0]));

	for(unsigned int iq = 1; iq < 8; iq++) {

		// Calculo del logaritmo comun (base 10) de la imagen
		dat.push_back(log10(data[iq]));

		// Obtencion de la mascara
		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			Mat mskir = (tmp & (1 << ir)) / (1 << ir);

			// Calcula de los desplazamientos relativos
			int dx = disp[iq][0] - disp[ir][0];
			int dy = disp[iq][1] - disp[ir][1];

			// Calculo de los extremos de las ventanas
			unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + ROWS; // FILAS
			unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + COLS; // COLUMNAS
			unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + ROWS; // FILAS
			unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + COLS; // COLUMNAS

			// Calcular ventanas de mascara. MskiqROI y mskirROI son del mismo tamaño, aunque estan desplazadas unas con respecto a la otra una distancia relativa.
			Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
			Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));

			// Calcular la mascara de las ventanas
			Mat mskDouble;
			Mat msk = mskiqROI.mul(mskirROI);

			msk.convertTo(mskDouble, CV_64F);

			cout << "-------------------------------------" << endl;
			cout << "Depth mskiq: " << mskiq.depth() << "\tDepth msk: " << msk.depth() << "\tDepth mskDouble: " << mskDouble.depth() << endl;

			cout << "Depth dat[iq]: " << dat[iq].depth() << "\tDepth dat[ir]: " << dat[ir].depth() << endl;
			//----------------------------------------------------------------------------------------
			cout << "Pasada " << ir << endl << "-------------------------------------" << endl;

			Size size = msk.size();
			int xmax = size.height;
			int ymax = size.width;

			cout << "Msk \tAncho: " << xmax << "\tAlto: " << ymax << endl;
			//----------------------------------------------------------------------------------------

			// Calcular ventanas de datos
			Mat datiqROI(dat[iq], Range(jyl, jyh), Range(jxl, jxh));
			Mat datirROI(dat[ir], Range(iyl, iyh), Range(ixl, ixh));
			//iq=log10(dataiq.at<uchar>(yActualQ, xActualQ))
			//----------------------------------------------------------------------
			size = datiqROI.size();
			xmax = size.height;
			ymax = size.width;
			cout << "datiqROI\tAncho: " << xmax << "\tAlto: " << ymax << endl;

			size = datirROI.size();
			xmax = size.height;
			ymax = size.width;
			cout << "datirROI\tAncho: " << xmax << "\tAlto: " << ymax << endl;
			//-----------------------------------------------------------------------


			// Calcular diferencia
			Mat diff = (datiqROI - datirROI).mul(mskDouble);

			//-----------------------------------------------------------------------
			size = diff.size();
			xmax = size.height;
			ymax = size.width;
			cout << "diff\tAncho: " << xmax << "\tAlto: " << ymax << endl;
			//-----------------------------------------------------------------------
/*
			imshow("diff",diff);

						waitKey( 0 );
*/
//----------------------------------------------
			double min,max;
			Point min_loc,max_loc;

			minMaxLoc(diff, &min, &max, &min_loc, &max_loc);

			Mat B;

			diff.convertTo(B, CV_32F);

			diff.convertTo(B,CV_8U,255.0/(max-min),-255.0/min);

			cout << "Min: " << min << "\tMax: " << max << endl;

//---------------------------------------------------

			// Calcular ventanas del termino constante
			Mat conJROI(con, Range(jyl, jyh), Range(jxl, jxh));
			Mat conIROI(con, Range(iyl, iyh), Range(ixl, ixh));

			// Aplicar la diferencia a las ventanas del termino constante
			conJROI = conJROI + diff;
			conIROI = conIROI - diff;

			// Calcular ventanas de la matriz de conteo de pares de pixeles
			Mat pixCntJROI(pixCnt, Range(jyl, jyh), Range(jxl, jxh));
			Mat pixCntIROI(pixCnt, Range(iyl, iyh), Range(ixl, ixh));

			// Aplicar la mascara a las ventanas de la matriz de pares de pixeles
			pixCntJROI = pixCntJROI + msk;
			pixCntIROI = pixCntIROI + msk;


			//---------------------------------------------------------
			double minVal, maxVal;
			minMaxLoc(mskDouble, &minVal, &maxVal); //find minimum and maximum intensities

			cout << "Depth CON: " << con.depth() << "\tDepth PixCnt: " << pixCnt.depth() << endl;
			cout << "Channels CON: " << con.channels() << "\tChannels PixCnt: " << pixCnt.channels() << endl;
			cout << "Min: " << minVal << "\tMax: " << maxVal << endl;
			cout << "\nValor mskDouble: " << con.at<double>(Point(260,260)) << endl;

			imshow("con",con);
			waitKey( 0 );

#ifdef PROGRESS

			cout << "getConst: Iteracion " << ++loopCnt << " de 28..." << endl;

#endif

		}

	}

	///////////////////////////
	int gh; cin >> gh;
	///////////////////////////
	return con;

}

//*************************************************************************************
//*************************************************************************************
void getCon(const Mat& mskiq,const Mat& mskir,const Mat& dataiq,const Mat& datair, Mat& pixCnt, Mat& con, int dx, int dy){

	// Calculo de los extremos de las ventanas
	unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + ROWS; // FILAS
	unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + COLS; // COLUMNAS
	unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + ROWS; // FILAS
	unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + COLS; // COLUMNAS

	unsigned int maxY = (jyh - jyl), maxX = (jxh - jxl);		//Calculamos el tamaño de las ROI
	unsigned int xActualQ, yActualQ, xActualR, yActualR;

	for (unsigned int y = 0; y < maxY; y++){
		for (unsigned int x = 0; x < maxX; x++){
			//Calculamos los indices de desplazamientos actuales de ambas regiones
			yActualQ = (y + jyl);	yActualR = (y + iyl);
			xActualQ = (x + jxl);	xActualR = (x + ixl);

			uchar msk = mskiq.at<uchar>(yActualQ, xActualQ) * mskir.at<uchar>(yActualR,xActualR);		//COMPROBAR RESULTADO!!!!!!!!!!!!!!!!!!

			//**************************Calculamos el CON**************************
			double mskDouble = msk;

			double logDataiQ = log10(dataiq.at<double>(yActualQ, xActualQ));
			double logDataiR = log10(datair.at<double>(yActualR, xActualR));

			//double diff = (logDataiQ - logDataiR) (AND o *) mskDouble;
			double diff = (logDataiQ - logDataiR) * mskDouble;		//COMPROBAR RESULTADO!!!!!!!!!!!!!!!!!!

			con.at<double>(yActualQ, xActualQ) = con.at<double>(yActualQ, xActualQ) + diff;
			con.at<double>(yActualR, xActualR) = con.at<double>(yActualR, xActualR) - diff;

			//*************************Calculamos el PixCnt************************
			pixCnt.at<uchar>(yActualQ, xActualQ) = pixCnt.at<uchar>(yActualQ, xActualQ) + msk;
			pixCnt.at<uchar>(yActualR, xActualR) = pixCnt.at<uchar>(yActualR, xActualR) + msk;
		}
	}
}



Mat getConst(vector<Mat>& data, \
             const Mat& tmp, \
             Mat& pixCnt, \
             int **disp) {

	Mat con(data[0].size(), CV_64F);
	vector<Mat> dat;

	unsigned int loopCnt = 0;

	for(unsigned int iq = 1; iq < 8; iq++) {

		// Obtencion de la mascara
		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			Mat mskir = (tmp & (1 << ir)) / (1 << ir);

			// Calcula de los desplazamientos relativos
			int dx = disp[iq][0] - disp[ir][0];
			int dy = disp[iq][1] - disp[ir][1];

			getCon(mskiq, mskir, data[iq], data[ir], pixCnt, con, dx, dy);

			cout << "Pasada " << iq << "," << ir << endl;
			imshow("con",con);

			waitKey( 0 );

#ifdef PROGRESS

			cout << "getConst: Iteracion " << ++loopCnt << " de 28..." << endl;

#endif

		}

	}

	/////////////ELIMINAR/////////////
	int gh; cin >> gh;
	//////////////////////////////////
	return con;

}

//*************************************************************************************
//*************************************************************************************


void doIterationNueva(const Mat& con, \
                 Mat& gain, \
                 const Mat& tmp, \
                 const Mat& pixCnt, \
                 const int disp[8][2]) {

	unsigned int loopCnt = 0;

	// Creacion de la ganancia temporal
	Mat gainTmp;
	con.copyTo(gainTmp);

	for(unsigned int iq = 1; iq < 8; iq++) {

		// Obtencion de la mascara
		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			Mat mskir = (tmp & (1 << ir)) / (1 << ir);

			// Calcula de los desplazamientos relativos
			int dx = disp[iq][0] - disp[ir][0];
			int dy = disp[iq][1] - disp[ir][1];

			// Calculo de los extremos de las ventanas
			unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + ROWS; // FILAS
			unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + COLS; // COLUMNAS
			unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + ROWS; // FILAS
			unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + COLS; // COLUMNAS

			// Calcular ventanas de mascara
			Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
			Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));

			// Calcular la mascara de las ventanas
			Mat msk = mskiqROI.mul(mskirROI);
			msk.convertTo(msk, CV_64F);

			// Calcular ventanas de ganancia y ganancia temporal
			Mat gainTmpJROI(gainTmp, Range(jyl, jyh), Range(jxl, jxh));
			Mat gainTmpIROI(gainTmp, Range(iyl, iyh), Range(ixl, ixh));
			Mat gainJROI(gain, Range(jyl, jyh), Range(jxl, jxh));
			Mat gainIROI(gain, Range(iyl, iyh), Range(ixl, ixh));

			// Modificar la ganancia temporal en base a la ganancia y la mascara
			gainTmpJROI = gainTmpJROI + gainIROI.mul(msk);
			gainTmpIROI = gainTmpIROI + gainJROI.mul(msk);

#ifdef PROGRESS

			cout << "doItera : IteraciÃ³n " << ++loopCnt << " de 28..." << endl;

#endif

		}

	}

	// Calcular ganancia unitaria para no dividir por cero
	Mat pixCntAux = max(pixCnt, 1.0);
	pixCntAux.convertTo(pixCntAux, CV_64F);

	gainTmp = gainTmp / pixCntAux;

	// Eliminar elementos a cero (de la matriz de pares de pixeles)
	Mat index = min(pixCnt, 1.0);
	index.convertTo(index, CV_64F);

	gainTmp = gainTmp.mul(index);

	// Calcular sumatorios
	double sum2 = sum(gainTmp)[0];
	double sum3 = sum(gainTmp.mul(gainTmp))[0];
	double nPix = sum(index)[0];

	// Eliminar elementos mas de 5-Sigma veces alejados de la media
	double ave2 = sum2 / nPix; //Ganancia media de la CCD
	double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);

	//Marco los elementos las posiciones de los elementos con un cero que son
	index = (abs(gainTmp - ave2) > fiveSigma) / 255;
	index.convertTo(index, CV_64F);

	sum2 = sum2 - sum(gainTmp.mul(index))[0];
	nPix = nPix - sum(index)[0];

	// Normalizar la tabla de ganancias
	ave2 = sum2 / nPix;

	gainTmp = gainTmp - ave2;

	// Devolver la tabla de ganancias
	gain = gainTmp;

}


//*************************************************************************************
//*************************************************************************************

void getGainTmp(const Mat& mskiq,const Mat& mskir, Mat& gainTmp, Mat& gain, int dx, int dy){

	// Calculo de los extremos de las ventanas
	unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + ROWS; // FILAS
	unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + COLS; // COLUMNAS
	unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + ROWS; // FILAS
	unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + COLS; // COLUMNAS

	unsigned int maxY = (jyh - jyl), maxX = (jxh - jxl);		//Calculamos el tamaño de las ROI
	unsigned int xActualQ, yActualQ, xActualR, yActualR;

	for (unsigned int y = 0; y < maxY; y++){
		for (unsigned int x = 0; x < maxX; x++){
			//Calculamos los indices de desplazamientos actuales de ambas regiones
			yActualQ = (y + jyl);	yActualR = (y + iyl);
			xActualQ = (x + jxl);	xActualR = (x + ixl);

			//Obtenemos el valor de la mascara en el pixel actual

			uchar msk = mskiq.at<uchar>(yActualQ, xActualQ) * mskir.at<uchar>(yActualR,xActualR);	//COMPROBAR RESULTADO!!!!!!!!!!!!!!!!!!

			double mskDouble = msk * REL8TO64;		//COMPROBAR RESULTADO!!!!!!!!!!!!!!!!!!

			//Obtenemos los valores de las ganancias en los pixeles actuales
			double gainTmpJ = gainTmp.at<double>(yActualQ, xActualQ);
			double gainTmpI = gainTmp.at<double>(yActualR, xActualR);
			double gainJ = gain.at<double>(yActualQ, xActualQ);
			double gainI = gain.at<double>(yActualR, xActualR);

			//Modificar la ganancia temporal en base a la ganancia y la mascara
			gainTmp.at<double>(yActualQ, xActualQ) = gainTmpJ + (gainI*mskDouble);
			gainTmp.at<double>(yActualR, xActualR) = gainTmpI + (gainJ*mskDouble);

			/*
			//Modificar la ganancia temporal en base a la ganancia y la mascara DE OTRA FORMA
			gainTmp.at<double>(yActualQ, xActualQ) = gainTmp.at<double>(yActualQ, xActualQ) + ( gain.at<double>(yActualR, xActualR) * mskDouble );
			gainTmp.at<double>(yActualR, xActualR) = gainTmp.at<double>(yActualR, xActualR) + ( gain.at<double>(yActualQ, xActualQ) * mskDouble );
			*/
		}
	}
}

void calculateStats(const Mat& pixCnt, Mat& gainTmp, Mat& gain){
	Mat pixCntAux = max(pixCnt, 1.0);
	pixCntAux.convertTo(pixCntAux, CV_64F);

	Mat index = min(pixCnt, 1.0);
	index.convertTo(index, CV_64F);

	Size s = gainTmp.size();
	int xmax = s.height;
	int ymax = s.width;
	for (int y = 0 ; y < ymax; y++)
		for (int x = 0; x < xmax; x++)
			gainTmp.at<double>(y, x) = gainTmp.at<double>(y, x) / pixCntAux.at<double>(y, x);

	double sum2 = 0, sum3 = 0, nPix = 0;
	for (int y = 0 ; y < ymax; y++){
		for (int x = 0; x < xmax; x++){
			sum2 += gainTmp.at<double>(y,x);												//COMPROBAR VAL[0]
			sum3 += (gainTmp.at<double>(y,x) * gainTmp.at<double>(y,x));					//COMPROBAR VAL[0]
			nPix += index.at<double>(y,x);													//COMPROBAR VAL[0]
		}
	}

	// Eliminar elementos mas de 5-Sigma veces alejados de la media
	double ave2 = sum2 / nPix; //Ganancia media de la CCD
	double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);

	//Marco los elementos las posiciones de los elementos con un cero que son
	for (int y = 0 ; y < ymax; y++){
		for (int x = 0; x < xmax; x++){
			index.at<double>(y,x) = (abs(gainTmp.at<double>(y,x) - ave2) > fiveSigma) / 255;
		}
	}
	index.convertTo(index, CV_64F);

	for (int y = 0 ; y < ymax; y++){
		for (int x = 0; x < xmax; x++){
			sum2 -= (gainTmp.at<double>(y,x) * index.at<double>(y,x));				//COMPROBAR VAL[0] EN AMBOS
			nPix -= (index.at<double>(y,x));										//COMPROBAR VAL[0]
		}
	}

	// Normalizar la tabla de ganancias
	ave2 = sum2 / nPix;
	for (int y = 0 ; y < ymax; y++){
		for (int x = 0; x < xmax; x++){
			gainTmp.at<double>(y,x) -= ave2;
		}
	}

	// Devolver la tabla de ganancias
	gain = gainTmp;
}

void doIteration(const Mat& con, \
					  Mat& gain, \
					  const Mat& tmp, \
					  const Mat& pixCnt, \
					  const int disp[8][2]) {

	unsigned int loopCnt = 0;

	// Creacion de la ganancia temporal
	Mat gainTmp;
	con.copyTo(gainTmp);

	for(unsigned int iq = 1; iq < 8; iq++) {

		// Obtencion de la mascara
		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);

		for(unsigned int ir = 0; ir < iq; ir++) {

			// Obtencion de la mascara
			Mat mskir = (tmp & (1 << ir)) / (1 << ir);

			// Calcula de los desplazamientos relativos
			int dx = disp[iq][0] - disp[ir][0];
			int dy = disp[iq][1] - disp[ir][1];

			getGainTmp(mskiq, mskir, gainTmp, gain, dx, dy);

#ifdef PROGRESS

			cout << "doItera : IteraciÃ³n " << ++loopCnt << " de 28..." << endl;

#endif

		}

	}

	calculateStats(pixCnt, gainTmp, gain);

}

//*************************************************************************************
//*************************************************************************************


Mat iterate(const Mat& con, \
            Mat& gain, \
            const Mat& tmp, \
            const Mat& pixCnt, \
            const int disp[8][2], \
			const unsigned int loops) {

	for(unsigned int i = 0; i < loops; i++) {

		doIteration(con, gain, tmp, pixCnt, disp);

#ifdef PROGRESS

		cout << "iterate : IteraciÃ³n " << i + 1 << " de " << loops << "..." << endl;

#endif

	}

	// Calculo de la imagen de flatfield
	Mat flat = gain * log(10.0);
	exp(flat, flat);

	Mat tmpAux = (tmp > 0) / 255;
	tmpAux.convertTo(tmpAux, CV_64F);

	flat = flat.mul(tmpAux);

	return flat;

}
