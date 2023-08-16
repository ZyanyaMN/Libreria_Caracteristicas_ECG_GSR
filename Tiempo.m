function [] = Tiempo(ecg,Fs)
%G=2, %Ganancia
%ecg = ecg/G;
ecg = (ecg - mean(ecg))/std(ecg);
t = (1:1:length(ecg))*(1/Fs);


umbral_y=6*mean(abs(ecg));
umbral_x=0.3*Fs;


%FinPeaks
[pks,locs]=findpeaks(ecg,'MinPeakHeight',umbral_y,'MinPeakDistance',umbral_x);

%%HRV
HRV=diff(t(locs));
new_t=downsample(t,125);
HRV_f=interp1(t(locs(1:end-1)),HRV,new_t);

%BPM
BPM=(1./HRV_f)*60;
BPM1=BPM(~isnan(BPM));
new_t1=new_t(~isnan(BPM));


%% MÈTRICAS

%RMSSD (Raíz cuadrada de la media de las diferencias de la suma de los cuadrados entre intervalos RR adyacentes)
%PRIMER PASO
resta_RR=diff(HRV);
%SEGUNDO
resta_RR2=resta_RR.^2;
%Tercer
suma_resta_RR2=sum(resta_RR2);
%Cuarto
norm_resta_RR2=suma_resta_RR2/length(resta_RR);
%Quinto
RMSSD=sqrt(norm_resta_RR2);


%(SDNN:Desviación estándar todos los intervalos R-R)
SDNN=std(HRV);
SDNN=SDNN*1000;

%Pasarlo a ms 
RMSSD_ms=RMSSD*1000;

%%Calcular PRR50 (PNN50:Número de intervalos adyacentes que varían por más de 50ms expresado en porcentaje)
abs_resta_RR=abs(diff(HRV));

%Logica booleana
num_max_50=abs_resta_RR>50/1000;

%Calculo final ('Porcentaje')
PRR50=sum(num_max_50)/length(abs_resta_RR);
PRR50=PRR50*100;


%AVNN (mean RR)
AVNN= mean(HRV);
AVNN=AVNN*1000;

%Promedio FC (FRECUENCIA CARDIACA)
HR= mean (BPM1);
%Mínimo BPM
BPM_MIN=min(BPM1);
%Máximo BPM
BPM_MAX=max(BPM1);


% Mostrar los resultados
disp(['RMSDD:' num2str(RMSSD_ms) 'ms']);
disp(['PRR50: ' num2str(PRR50) '%']);
disp(['AVNN:' num2str(AVNN) 'ms']);
disp(['PROMEDIO BPM:' num2str(HR) '1/min']);
disp(['SDNN:' num2str(SDNN) 'ms']);

end