imsize=2048;
paso=0.1;
C_x=1024;
C_y=1020;
dim=10;
Xmi=C_x-paso*(dim/2)
Xmx=C_x+paso*(dim/2)
Ymi=C_y-paso*(dim/2)
Ymx=C_y+paso*(dim/2)
difx=Xmx-Xmi
dify=Ymx-Ymi
a=[Xmi:paso:Xmx];

aa=round((a-Xmi)/paso)

max(aa)*paso+Xmi
i=Xmi
Xmx
h=1
    while i<=Xmx
    det=i*i
    h=h+1
    i=i+paso    
    end
    %lmax = (double)PORCENTAJE/100 * xmax;
    lmax = round(15/100 * 2048) % numero de radios
	if mod(lmax,2) == 0
		lmax=lmax+1
    end    
    
  %  //Recorremos lmax veces la imagenes cambiando el radio
%	for(int l = 1; l <= lmax; l++){
%		//Calculamos el radio de esta pasada
%		r=radio + (l - ((int)(lmax/2)+1))*FACTOR_RADIO;
r=zeros(1,lmax)
radio=100
FACTOR_RADIO=0.25
for i=1:lmax
    r(i)=radio+(i - (round(lmax/2)+1))*FACTOR_RADIO;
end