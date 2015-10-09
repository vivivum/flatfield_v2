 M = dlmread('EstadisticasdeCentro.txt', '\t');
 MR = dlmread('EstadisticasdeCentro_rotados.txt', '\t');
 
 M=M(:,1:4);
 MR=MR(:,1:4);
 
 centrosBuenos = [
           1022     1021 
           990.683  1315.75      
           1200.68  1251.75 
           1330.68  1055.75
           1230.68  807.750
           1012.68  681.750 
           800.683  801.750   
           730.683  1031.75 
           850.683  1291.75   
      ];
  
%  centrosBuenos = zeros(9,2);
%  centrosBuenos(:,1) = centrosBuenosOr(:,2);
%  centrosBuenos(:,2) = centrosBuenosOr(:,1);
 
 numeroImagenes = 9;
 porcentaje = 0.8;
 
 
 im = zeros(103,4,numeroImagenes);
 umbral = zeros(numeroImagenes);
 
 %para imagenes rotadas
 imR = zeros(103,4,numeroImagenes);
 umbralR = zeros(numeroImagenes);
 
 centros = zeros(numeroImagenes,3);%coordenadas del centro x u y radio
 cmax = zeros(numeroImagenes,3);
 cPonderados = zeros(numeroImagenes,3);%coordenadas del centro x u y radio
 erroresMedia = zeros(numeroImagenes, 2);
 erroresPico = zeros(numeroImagenes, 2);
 erroresPond = zeros(numeroImagenes, 2);
 
 
 centrosR = zeros(numeroImagenes,3);
 cmaxR = zeros(numeroImagenes,3);
 cPonderadosR = zeros(numeroImagenes,3);
 erroresMediaR = zeros(numeroImagenes, 2);
 erroresPicoR = zeros(numeroImagenes, 2);
 erroresPondR = zeros(numeroImagenes, 2);
 
 
 for i = 0:numeroImagenes-1
     im(:,:,i+1) = M((1+103*i):((1+103*i)+102),:);
     umbral(i+1) = porcentaje * max(im(:,1,i+1)) + (1-porcentaje)*min(im(:,1,i+1));
     A=im(:,:,i+1);
     A0=A( find( A(:,1) > umbral(i+1) ),:); 
%      size(A0)
%      A0
%      pause;
     umbralMax = max(im(:,1,i+1));
     
     F=A( find( A(:,1) == umbralMax ),:); 
     
     centros(i+1,1) = mean(A0(:,2));
     centros(i+1,2) = mean(A0(:,3));
     centros(i+1,3) = mean(A0(:,4));
     
     maxPond = sum(A0(:,1));
     %media ponderada
     cPonderados(i+1,1) = sum(A0(:,1).*A0(:,2))/maxPond;
     cPonderados(i+1,2)  = sum(A0(:,1).*A0(:,3))/maxPond;
     cPonderados(i+1,3)  = sum(A0(:,1).*A0(:,4))/maxPond;
  
     
     cmax(i+1,1) =  F(2);
     cmax(i+1,2) =  F(3);
     cmax(i+1,3) =  F(4);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     imR(:,:,i+1) = MR((1+103*i):((1+103*i)+102),:);
     umbralR(i+1) = porcentaje * max(imR(:,1,i+1)) + (1-porcentaje)*min(imR(:,1,i+1));
     AR=imR(:,:,i+1);
     A0R=AR( find( AR(:,1) > umbralR(i+1) ),:); 
%      size(A0)
%      A0
%      pause;
     umbralMaxR = max(imR(:,1,i+1));
     
     FR=AR( find( AR(:,1) == umbralMaxR ),:); 
     
     centrosR(i+1,1) = mean(A0R(:,2));
     centrosR(i+1,2) = mean(A0R(:,3));
     centrosR(i+1,3) = mean(A0R(:,4));
     
     maxPondR = sum(A0R(:,1));
     
     cPonderadosR(i+1,1) = sum(A0R(:,1).*A0R(:,2))/maxPondR;
     cPonderadosR(i+1,2)  = sum(A0R(:,1).*A0R(:,3))/maxPondR;
     cPonderadosR(i+1,3)  = sum(A0R(:,1).*A0R(:,4))/maxPondR;
  
     
     cmaxR(i+1,1) =  FR(2);
     cmaxR(i+1,2) =  FR(3);
     cmaxR(i+1,3) =  FR(4);
     
 end    
     
 
 temp = centrosR(:,1);
 centrosR(:,1) = centrosR(:,2);
 centrosR(:,2) = temp;
 
 temp = cmaxR(:,1);
 cmaxR(:,1) = cmaxR(:,2);
 cmaxR(:,2) = temp;
 
 temp = cPonderadosR(:,1);
 cPonderadosR(:,1) = cPonderadosR(:,2);
 cPonderadosR(:,2) = temp;
 
 erroresMedia = centrosBuenos - centros(:,1:2)
 erroresPico = centrosBuenos - cmax(:,1:2)
 erroresPond = centrosBuenos - cPonderados(:,1:2)
 
 
 erroresMedia_C = centrosBuenos - (centros(:,1:2) + centrosR(:,1:2))/2
 erroresPico_C = centrosBuenos - (cmax(:,1:2) + cmaxR(:,1:2))/2
 erroresPond_C = centrosBuenos - (cPonderados(:,1:2) + cPonderadosR(:,1:2))/2
 
 
 
 
 
 
 
 
 
 
 
 
 
 