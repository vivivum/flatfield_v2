//#include "flatfield.hpp"
//#include "utility.hpp"
//
//#include <algorithm>
//#include <cmath>
//#include <iostream>
//
///*
//Mat binarizar (const Mat& im_in){
//	Mat im_out;
//	//Hacemos una copia de la imagen de entrada
//	im_out=im_in;
//
//	int im_out_p[] = {1,1,1,1,1,1,1,1,1};
//
//	int k[][] = {{1,1,1},{1,1,1},{1,1,1}};
//
//	Size s = im_out.size();
//	int rows = s.height;
//	int cols = s.width;
//
//	int z = 0;
//	for(int i = 1; i < rows-2; i++){
//		for(int i = 1; i < cols-2; i++){
//			z = 1;
//			for(int ik=-1; ik <= 1; ik++){
//				for(int jk=-1; jk <= 1; jk++){
//					im_out_p[z] = im_out.at<uint>(j+jk,i+ik) * k[ik+2][jk+2];
//					z++;
//				}
//			}
//			im_out.at<uint>(j,i) = min(im_out_p);
//		}
//	}
//
//	Mat diff = im_in - im_out;
//	return diff;
//}
//*/
//int getImagesOriginal(vector<Mat>& data, \
//              Mat& tmp, \
//              const double rMin, \
//              const double rMax) {
//
//	// Comprobar el numero de imagenes
//	unsigned int nImages = data.size();
//
//	if(nImages > 8) {
//
//		cerr << "Error - El numero de imagenes deben de ser ocho..." << endl;
//
//		return 1;
//
//	}
//
//	// Calcular la plantilla de pixeles buenos
//	for(unsigned int i = 0; i < 8; i++) {
//
//		// Calcular la mascara, dividiendo por 255 debido a que, en OpenCV, TRUE
//		// equivale a 255
//		Mat msk = ((data[i] >= rMin) & (data[i] <= rMax)) / 255;
//
//		cout << "Imagen " << i + 1 << ": " << (unsigned int)(sum(msk).val[0]) << " pixeles buenos..." << endl;
//
//		// Guardar la mascara en la plantilla de pixeles buenos
//		tmp = tmp | (msk * (1 << i));
//
//		// Emplear la mascara para cribar los pixeles malos
//		msk.convertTo(msk, data[i].type());
//		data[i] = data[i].mul(msk);
//
//	}
//
//	return 0;
//
//}
//
//Mat getConstOriginal(vector<Mat>& data, \
//             const Mat& tmp, \
//             Mat& pixCnt, \
//             const int disp[8][2]) {
//
//	Mat con(data[0].size(), CV_64F);
//	vector<Mat> dat;
//
//	unsigned int loopCnt = 0;
//
//	// Calculo del logaritmo comun (base 10) de la imagen
//	dat.push_back(log10(data[0]));
//
//	for(unsigned int iq = 1; iq < 8; iq++) {
//
//		// Calculo del logaritmo comun (base 10) de la imagen
//		dat.push_back(log10(data[iq]));
//
//		// Obtencion de la mascara
//		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);
//
//		for(unsigned int ir = 0; ir < iq; ir++) {
//
//			// Obtencion de la mascara
//			Mat mskir = (tmp & (1 << ir)) / (1 << ir);
//
//			// Calcula de los desplazamientos relativos
//			int dx = disp[iq][0] - disp[ir][0];
//			int dy = disp[iq][1] - disp[ir][1];
//
//			// Calculo de los extremos de las ventanas
//			unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + data[0].rows; // FILAS
//			unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + data[0].cols; // COLUMNAS
//			unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + data[0].rows; // FILAS
//			unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + data[0].cols; // COLUMNAS
//
//			// Calcular ventanas de mascara
//			Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
//			Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Calcular la mascara de las ventanas
//			Mat mskDouble;
//			Mat msk = mskiqROI.mul(mskirROI);
//			msk.convertTo(mskDouble, CV_64F);
//
//			// Calcular ventanas de datos
//			Mat datiqROI(dat[iq], Range(jyl, jyh), Range(jxl, jxh));
//			Mat datirROI(dat[ir], Range(iyl, iyh), Range(ixl, ixh));
//
//			// Calcular diferencia
//			Mat diff = (datiqROI - datirROI).mul(mskDouble);
//
//			// Calcular ventanas del termino constante
//			Mat conJROI(con, Range(jyl, jyh), Range(jxl, jxh));
//			Mat conIROI(con, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Aplicar la diferencia a las ventanas del termino constante
//			conJROI = conJROI + diff;
//			conIROI = conIROI - diff;
//
//			// Calcular ventanas de la matriz de conteo de pares de pixeles
//			Mat pixCntJROI(pixCnt, Range(jyl, jyh), Range(jxl, jxh));
//			Mat pixCntIROI(pixCnt, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Aplicar la mascara a las ventanas de la matriz de pares de pixeles
//			pixCntJROI = pixCntJROI + msk;
//			pixCntIROI = pixCntIROI + msk;
//
//#ifdef PROGRESS
//
//			cout << "getConst: IteraciÃ³n " << ++loopCnt << " de 28..." << endl;
//
//#endif
//
//		}
//
//	}
//
//	return con;
//
//}
//
//
//
//
//
//
//int getImages(vector<Mat>& data, \
//              Mat& tmp, \
//              const double rMin, \
//              const double rMax) {
//
//	// Comprobar el numero de imagenes
//	unsigned int nImages = data.size();
//
//	int h;//Borrar SOLO DEBUG
//
//	if(nImages > 8) {
//
//		cerr << "Error - El numero de imagenes deben de ser ocho..." << endl;
//
//		return 1;
//
//	}
//
//	//Obtenemos las dimensiones de una imagen cualquiera (Todas deben ser del mismo tamaño)
//	Size s = data[0].size();
//	int xmax = s.height;
//	int ymax = s.width;
//
//	// Calcular la plantilla de pixeles buenos
//	for(unsigned int i = 0; i < 8; i++) {
//
//		Mat msk(xmax,ymax,CV_8U,Scalar(0));
//
//		for(int x=0; x < xmax; x++){
//			for(int y=0; y < ymax; y++){
//				if(data[i].at<double>(y,x) >= rMin && data[i].at<double>(y,x) <= rMax){
//					msk.at<uchar>(y,x) = 1;
//				}
//			}
//		}
//
//		cout << "Imagen " << i + 1 << ": " << (unsigned int)(sum(msk).val[0]) << " pixeles buenos..." << endl;
//
//		// Guardar la mascara en la plantilla de pixeles buenos
//		tmp = tmp | (msk * (1 << i));
//
//		// Emplear la mascara para cribar los pixeles malos
//		msk.convertTo(msk, data[i].type());
//		data[i] = data[i].mul(msk);
//	}
//
//	return(0);
//}
//
//
//
//Mat getConst(vector<Mat>& data, \
//             const Mat& tmp, \
//             Mat& pixCnt, \
//             int **disp) {
//
//	Mat con(data[0].size(), CV_64F);
//	vector<Mat> dat;
//
//	unsigned int loopCnt = 0;
//
//	// Calculo del logaritmo comun (base 10) de la imagen
//	dat.push_back(log10(data[0]));
//
//	for(unsigned int iq = 1; iq < 8; iq++) {
//
//		// Calculo del logaritmo comun (base 10) de la imagen
//		dat.push_back(log10(data[iq]));
//
//		// Obtencion de la mascara
//		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);
//
//		for(unsigned int ir = 0; ir < iq; ir++) {
//
//			// Obtencion de la mascara
//			Mat mskir = (tmp & (1 << ir)) / (1 << ir);
//
//			// Calcula de los desplazamientos relativos
//			int dx = disp[iq][0] - disp[ir][0];
//			int dy = disp[iq][1] - disp[ir][1];
//
//			// Calculo de los extremos de las ventanas
//			unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + data[0].rows; // FILAS
//			unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + data[0].cols; // COLUMNAS
//			unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + data[0].rows; // FILAS
//			unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + data[0].cols; // COLUMNAS
//
//			// Calcular ventanas de mascara. MskiqROI y mskirROI son del mismo tamaño, aunque estan desplazadas unas con respecto a la otra una distancia relativa.
//			Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
//			Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));
//
//
//			//Mostramos las ROI
//			Mat regiones = Mat::zeros( mskiq.size(), CV_8UC3 );
//			rectangle( regiones, Point( jyl, jxl ), Point( jyh, jxh), Scalar( 0, 55, 255 ), 1,  4);
//			rectangle( regiones, Point( iyl, ixl ), Point( iyh, ixh), Scalar( 255, 55, 255 ), 1,  4);
//
//
//
//
//
//
//			// Calcular la mascara de las ventanas
//			Mat mskDouble;
//			Mat msk = mskiqROI.mul(mskirROI);
//			msk.convertTo(mskDouble, CV_64F);
//
//			//----------------------------------------------------------------------------------------
//			cout << "Pasada " << ir << endl << "-------------------------------------" << endl;
//
//			Size size = msk.size();
//			int xmax = size.height;
//			int ymax = size.width;
//
//			cout << "Msk \tAncho: " << xmax << "\tAlto: " << xmax << endl;
//			//----------------------------------------------------------------------------------------
//
//			// Calcular ventanas de datos
//			Mat datiqROI(dat[iq], Range(jyl, jyh), Range(jxl, jxh));
//			Mat datirROI(dat[ir], Range(iyl, iyh), Range(ixl, ixh));
//
//			//----------------------------------------------------------------------
//			size = datiqROI.size();
//			xmax = size.height;
//			ymax = size.width;
//			cout << "datiqROI\tAncho: " << xmax << "\tAlto: " << xmax << endl;
//
//			size = datirROI.size();
//			xmax = size.height;
//			ymax = size.width;
//			cout << "datirROI\tAncho: " << xmax << "\tAlto: " << xmax << endl;
//			//-----------------------------------------------------------------------
//
//
//			// Calcular diferencia
//			Mat diff = (datiqROI - datirROI).mul(mskDouble);
//
//			//-----------------------------------------------------------------------
//			size = diff.size();
//			xmax = size.height;
//			ymax = size.width;
//			cout << "diff\tAncho: " << xmax << "\tAlto: " << xmax << endl;
//			//-----------------------------------------------------------------------
///*
//			imshow("diff",diff);
//
//						waitKey( 0 );
//*/
//
//			double min,max;
//			Point min_loc,max_loc;
//
//			minMaxLoc(diff, &min, &max, &min_loc, &max_loc);
//
//			Mat B;
//
//			diff.convertTo(B, CV_32F);
//
//			diff.convertTo(B,CV_8U,255.0/(max-min),-255.0/min);
//
//			cout << "Min: " << min << "\tMax: " << max << endl;
//
//
//
//			// Calcular ventanas del termino constante
//			Mat conJROI(con, Range(jyl, jyh), Range(jxl, jxh));
//			Mat conIROI(con, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Aplicar la diferencia a las ventanas del termino constante
//			conJROI = conJROI + diff;
//			conIROI = conIROI - diff;
//
//			// Calcular ventanas de la matriz de conteo de pares de pixeles
//			Mat pixCntJROI(pixCnt, Range(jyl, jyh), Range(jxl, jxh));
//			Mat pixCntIROI(pixCnt, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Aplicar la mascara a las ventanas de la matriz de pares de pixeles
//			pixCntJROI = pixCntJROI + msk;
//			pixCntIROI = pixCntIROI + msk;
//
//			imshow("con",con);
//
//			waitKey( 0 );
//
//#ifdef PROGRESS
//
//			cout << "getConst: Iteracion " << ++loopCnt << " de 28..." << endl;
//
//#endif
//
//		}
//
//	}
//
//	///////////////////////////
//	int gh; cin >> gh;
//	///////////////////////////
//	return con;
//
//}
//
//
//void doIteration(const Mat& con, \
//                 Mat& gain, \
//                 const Mat& tmp, \
//                 const Mat& pixCnt, \
//                 const int disp[8][2]) {
//
//	unsigned int loopCnt = 0;
//
//	// Creacion de la ganancia temporal
//	Mat gainTmp;
//	con.copyTo(gainTmp);
//
//	for(unsigned int iq = 1; iq < 8; iq++) {
//
//		// Obtencion de la mascara
//		Mat mskiq = (tmp & (1 << iq)) / (1 << iq);
//
//		for(unsigned int ir = 0; ir < iq; ir++) {
//
//			// Obtencion de la mascara
//			Mat mskir = (tmp & (1 << ir)) / (1 << ir);
//
//			// Calcula de los desplazamientos relativos
//			int dx = disp[iq][0] - disp[ir][0];
//			int dy = disp[iq][1] - disp[ir][1];
//
//			// Calculo de los extremos de las ventanas
//			unsigned int jyl = max(0, -dy), jyh = min(0, -dy) + con.rows; // FILAS
//			unsigned int jxl = max(0, -dx), jxh = min(0, -dx) + con.cols; // COLUMNAS
//			unsigned int iyl = max(0,  dy), iyh = min(0,  dy) + con.rows; // FILAS
//			unsigned int ixl = max(0,  dx), ixh = min(0,  dx) + con.cols; // COLUMNAS
//
//			// Calcular ventanas de mascara
//			Mat mskiqROI(mskiq, Range(jyl, jyh), Range(jxl, jxh));
//			Mat mskirROI(mskir, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Calcular la mascara de las ventanas
//			Mat msk = mskiqROI.mul(mskirROI);
//			msk.convertTo(msk, CV_64F);
//
//			// Calcular ventanas de ganancia y ganancia temporal
//			Mat gainTmpJROI(gainTmp, Range(jyl, jyh), Range(jxl, jxh));
//			Mat gainTmpIROI(gainTmp, Range(iyl, iyh), Range(ixl, ixh));
//			Mat gainJROI(gain, Range(jyl, jyh), Range(jxl, jxh));
//			Mat gainIROI(gain, Range(iyl, iyh), Range(ixl, ixh));
//
//			// Modificar la ganancia temporal en base a la ganancia y la mascara
//			gainTmpJROI = gainTmpJROI + gainIROI.mul(msk);
//			gainTmpIROI = gainTmpIROI + gainJROI.mul(msk);
//
//#ifdef PROGRESS
//
//			cout << "doItera : IteraciÃ³n " << ++loopCnt << " de 28..." << endl;
//
//#endif
//
//		}
//
//	}
//
//	// Calcular ganancia unitaria
//	Mat pixCntAux = max(pixCnt, 1.0);
//	pixCntAux.convertTo(pixCntAux, CV_64F);
//
//	gainTmp = gainTmp / pixCntAux;
//
//	// Eliminar elementos a cero (de la matriz de pares de pixeles)
//	Mat index = min(pixCnt, 1.0);
//	index.convertTo(index, CV_64F);
//
//	gainTmp = gainTmp.mul(index);
//
//	// Calcular sumatorios
//	double sum2 = sum(gainTmp)[0];
//	double sum3 = sum(gainTmp.mul(gainTmp))[0];
//	double nPix = sum(index)[0];
//
//	// Eliminar elementos mas de 5-Sigma veces alejados de la media
//	double ave2 = sum2 / nPix;
//	double fiveSigma = 5 * sqrt((sum3 / nPix) - ave2 * ave2);
//
//	index = (abs(gainTmp - ave2) > fiveSigma) / 255;
//	index.convertTo(index, CV_64F);
//
//	sum2 = sum2 - sum(gainTmp.mul(index))[0];
//	nPix = nPix - sum(index)[0];
//
//	// Normalizar la tabla de ganancias
//	ave2 = sum2 / nPix;
//
//	gainTmp = gainTmp - ave2;
//
//	// Devolver la tabla de ganancias
//	gain = gainTmp;
//
//}
//
//
//Mat iterate(const Mat& con, \
//            Mat& gain, \
//            const Mat& tmp, \
//            const Mat& pixCnt, \
//            const int disp[8][2], \
//			const unsigned int loops) {
//
//	for(unsigned int i = 0; i < loops; i++) {
//
//		doIteration(con, gain, tmp, pixCnt, disp);
//
//#ifdef PROGRESS
//
//		cout << "iterate : IteraciÃ³n " << i + 1 << " de " << loops << "..." << endl;
//
//#endif
//
//	}
//
//	// Calculo de la imagen de flatfield
//	Mat flat = gain * log(10.0);
//	exp(flat, flat);
//
//	Mat tmpAux = (tmp > 0) / 255;
//	tmpAux.convertTo(tmpAux, CV_64F);
//
//	flat = flat.mul(tmpAux);
//
//	return flat;
//
//}
