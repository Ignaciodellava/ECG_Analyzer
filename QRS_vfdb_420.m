% PRÁCTICA 4 DETECTOR DE LATIDOS DE ECG

%  IGNACIO DE LLANO VARELA & CESAR CARRASCO





clear
clc

%----------------------------Actividad 1-----------------------------------
% Descargamos la señal
%Si la frecuencia de muestreo son 250 muestras por segundo, un minuto seran
%15000 muestras
[signal,Fs,tm]=rdsamp('vfdb/420',[],15000,0);
signal=signal(:,1);
figure(1),plot(tm,signal);
%Guardamos un solo canal




%----------------------------Actividad 2-----------------------------------
%T.F
tf=fft(signal);
figure(2),plot(abs(tf).^2)%abs es valor absoluto porque hay valores negativos

% en vez de entre -x y x esta pasado el cero y luego por ahi, hay que
% pasarlo ese trozo a -x transformada de fourier desplazada
tfshift=fftshift(tf);
ejef=linspace(-Fs/2, Fs/2, length(tfshift));
figure(3),plot(ejef, abs(tfshift/length(signal)));

%La centramos en cero
continua=mean(signal);
lineacont=ones(1,length(signal))*continua;
figure(4),plot(tm,signal,tm,lineacont,'r');
ssc=signal-continua;
figure(5),plot(tm,signal,tm,ssc,'r');





%----------------------------Actividad 3-----------------------------------

%Filtro Paso Bajo para eliminar altas frecuencias
%fc=11; 
fc=25; %Escogemos una frecuencia de corte de 25 Hz para obtener un resultado mejor en vez de entre 5 y 15
w=fc/(Fs/2);
[b,a]=butter(6,w);
freqz(b,a);
x=filter(b,a,ssc);
figure(5),plot(x,'b') ,hold on, plot(signal,'r')
%Observamos que hemos eliminado el ruido de alta frecuencia, los piquitos
%pequeños


tf2=fft(x);%T.F
figure(6),plot(abs(tf2).^2)%abs es valor absoluto porque hay valores negativos
tfshift2=fftshift(tf2);
ejef2=linspace(-Fs/2, Fs/2, length(tfshift2));
figure(7),plot(ejef2, abs(tfshift2/length(x))), hold on,plot(ejef, abs(tfshift/length(signal)),'r')
%En la representacion en frecuencia vemos que efectivamente el ruido se ha
%reducido y se ha eliminado el de valores mas altos

%Filtro Paso Alto para eliminar las bajas frecuencias no deseadas
fc2=5;
w2=fc2/(Fs/2);
[b2,a2]=butter(6,w2,'high');
freqz(b2,a2);
x2=filter(b2,a2,x);
figure(8),plot(x,'b') ,hold on, plot(x2,'r')
%Aplicando el filtro alteramos un poco la señal pero conservamos las curvas
%mas pronunciadas las del qrs.


%----------------------------Actividad 4-----------------------------------
%Derivacion
a4=8;
b4=[2 1 0 -1 -2];
x4=filter(b4,a4,x2);
figure(9),plot(x4,'b') ,hold on, plot(x2,'r')
%Resaltamos los cambios de amplitud

%Lo elevaos al cuadrado para tener solo valores positivos y aumentar la
%amplitud de los picos que queremos detectar
x5=x4.^2;
figure(10),plot(x5,'b')


%Transformamos los picos en una sola onda mas fácil de detectar
vit=0.150*Fs;
vit=round(vit);
b5=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
a5=vit;
x6=filter(b5,a5,x5);
figure(11),plot(x6,'b')




%----------------------------Actividad 5-----------------------------------

help rdann
help findpeaks

%Normalizamos la señal para que sea mas sencilla y universal de tratar
x6=x6/max(x6);
figure(11),plot(x6,'b')

%Detectamos sobre la señal obtenida los picos de la señal con un mínimo de
%0.3 de amplitud ya que los qrs de menor valor que hemos observado en
%nuestras pruebas despues de normalizarlos tenian esos valores. Añadimos
%una distancia minima de 50 entre deteccion de latido y latido y  un 0.25
%de caida abrupta para asegurar ue es un pico de la señal que sube y baja.
[valor, localizacion]=findpeaks(x6, 'MinPeakHeight',0.3,'MinPeakProminence',0.25,'MinPeakDistance',50);

%     [...] = findpeaks(...,'MinPeakHeight',MPH) finds only those peaks that
%     are greater than the minimum peak height, MPH. MPH is a real-valued
%     scalar. The default value of MPH is -Inf.
 
%     [...] = findpeaks(...,'MinPeakProminence',MPP) finds peaks guaranteed
%     to have a vertical drop of more than MPP from the peak on both sides
%     without encountering either the end of the signal or a larger
%     intervening peak. The default value of MPP is zero.
 
%     [...] = findpeaks(...,'MinPeakDistance',MPD) finds peaks separated by
%     more than the minimum peak distance, MPD. This parameter may be
%     specified to ignore smaller peaks that may occur in close proximity to
%     a large local peak. For example, if a large local peak occurs at LOC,
%     then all smaller peaks in the range [N-MPD, N+MPD] are ignored. If not
%     specified, MPD is assigned a value of zero.  





%Creación de eje de tiempo para observar si escogemos bien los picos y
%eliminamos otros que no queramos luego del umbral.
eje=linspace(0,100,length(x6));

%Representamos la señal final con los picos que ha detectado
plot(eje,x6), hold on, plot(eje(localizacion),x6(localizacion),'r*')
%Vemos a que posicion corresponden de la señal original
plot(eje,signal), hold on, plot(eje(localizacion),signal(localizacion),'r*')

%----------------------------Actividad 5b-----------------------------------
%Rudimentario
%Hacemos una estimacion de la longitud del qrs
finalqrs=localizacion;
principioqrs=finalqrs-43;
%NEGRO=inicio
%AMARILLO=final
plot(eje,signal), hold on, plot(eje(finalqrs),signal(finalqrs),'y*'),hold on, plot(eje(principioqrs),signal(principioqrs),'k*')


%-------------------------------Actividad 6-------------------------------------------

%Utilizamos la función rdann para obtener los datos de las anotaciones de
%la base de datos physionet que son los supuestamente acertados
 ann=rdann('vfdb/420','atr',[],15000,0);

eje=linspace(0,100,length(signal));

%Verdes detectados con ann, rojos, los nuestros detectados con el procesamiento de la
%señal. Representamos los rojos para comparar si es mas o menos acertado
%el algoritmo. Si es así procedems a calcular los valores de la detección de
%cada latido. 
plot(eje,signal),hold on, plot(eje(ann),signal(ann), 'g*'),hold on, plot(eje(localizacion),x6(localizacion),'r*');

%Por alguna razon esta señal solo tiene una anotacion en physionet



%Un ligero desplazamiento también es válido
%Ahora intentaremos calcular la sensibilidad y el valor predictivo positivo
%de este algoritmo de detección de latidos. Lo haremos creano un rango de
%muestras mas y menos que la anotación de physionet, segun nuestro algoritmo, 
% hay una variacion de aproximadamente 60 muestras. si nuestra detección
%de latidos, la roja, esta dentro del rango de muestras sumadas y restadas
%al valor verde detectado por physionet, será un Verdadero Positivo, si
%está fuera será un Falso Positivo, y si no lo hay en ningun caso, pero
%physionet si tiene uno indicado será un Falso negativo

%Creamos los valores vacíos
 FN=[];
 detec=[];
 FP=[];
 VP=[];

for i=1:length(localizacion)
    if ann(find(ann>=(localizacion(i)-60) & ann<= (localizacion(i)+60)))
        detec=[detec localizacion(i)];
        if isempty(detec)
        FP=[FP localizacion(i)];
    else
        VP =[VP localizacion(i)];
        end
    else 
    end
end

for j=1:length(ann)
    if localizacion(find(localizacion>=(ann(j)-60)&localizacion<=(ann(j)+60)))
 detec=[detec ann(j)];
        if isempty(detec)
        end
    else
        FN =[FN ann(j)];
        end
end

 Sensibilidad= length(VP)/(length(VP)+length(FN));
 VPP=length(VP)/(length(VP)+length(FP));









%--------------------------------------------------------------------------------
% i=1;
% while i<length(picos)
%      FP=[];
%     VP=[];
%     
% detec=find(ann>picos(i)-15);
% detec =detec(detec<picos(i)+15);
% 
%   if isempty(detec)
%         FP=[FP picos(i)];
%     else
%         VP =[VP detec];
%     end
% end


%--------------------------------------------------------------------------------
% while 1i == 1:length(picos)
%     FP=[];
%     VP=[];
% 
%     detec =ann(ann>picos(1i)-10);            %15 //Margen
%     detec=detec(detec<picos(1i)+10);
% 
%     if isempty(detec)
%         FP=[FP,picos(1i)];
%     else
%         VP =[VP,detec];
%     end
%  
% end
% 
% 
% while 2i == 1:length(ann)
%     FN=[];
%     detec =ann(picos>ann(2i)-10);  %15 //Margen
%     detec=detec(detec<ann(2i)+10);
% 
%     if isempty(detec)
%         FN=[FN,ann(1i)];
%     end
% end





