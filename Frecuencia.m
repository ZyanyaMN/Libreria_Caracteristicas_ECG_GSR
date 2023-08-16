function [] = Frecuencia(ecg,Fs)
f_vlf = [0.000, 0.04];   % Frecuencia muy baja (VLF) en Hz
f_lf = [0.04, 0.15];     % Frecuencia baja (LF) en Hz
f_hf = [0.15, 0.4];      % Frecuencia alta (HF) en Hz
window_minutes = 4.27; %tiempo de ventana (256 segundos) 
overlap = .5; % 50% solapamiento
resample = 4; %Interpolación es de 4 hz 

%FinPeaks
t = (1:1:length(ecg))*(1/Fs);
umbral_y=6*mean(abs(ecg));
umbral_x=0.3*Fs;
[pks,locs]=findpeaks(ecg,'MinPeakHeight',umbral_y,'MinPeakDistance',umbral_x);

%HRV
HRV=diff(t(locs));

tnn=cumsum(HRV); %Calculo de vector de tiempo de la señal
RRA= HRV-mean(HRV); % eliminar ruido DC

t_max = tnn(end); %Tiempo máximo de la señal
f_max=1; %1 Hz máximo de frecuencia


t_win =round(60 *window_minutes); % Ventana de 4.27 minutos*60=256 segundos 

num_windows = round(t_max / t_win);
if (num_windows < 1)
    num_windows = 1;
    t_win = floor(tnn(end)-tnn(1));
end

%Reconstrucción 4 hz 
fs_uni = resample;   %Hz
steps=1/resample; %espacios 

tnn_uni = tnn(1) : steps : tnn(end); %nuevo vector tiempo
n_win_uni = round(t_win / steps); % muestras por ventana 
num_windows_uni = round(length(tnn_uni) / n_win_uni); %número de ventanas por muestras


ts = t_win / (n_win_uni-1);   % Tiempo de muestreo 
f_res = 1 / (n_win_uni * ts); % Resolución espectral   
f_res = f_res / resample;     % Interpolación de muestras 
% 
% f_start = 0; % Frecuencia inicial
% f_end = 100; % Frecuencia final
% num_points = n_win_uni; % Número de puntos en el vector
% 
% f_axis= (linspace(f_start, f_res, Fs)');
f_axis = (f_res : f_res :f_max)'; %frecuencias ciclos por unidad de tiempo

nni_uni = interp1(tnn, RRA, tnn_uni, 'spline')'; %Interpolación 4 Hz
%window = hamming(n_win_uni, 'periodic');
window =hamming(n_win_uni); %Ventana
[Pxx,frequencies] = pwelch(nni_uni, window, overlap,f_axis, fs_uni); %PWelch
pxx  = Pxx * (1/mean(window)); % gain correction

% Encontrar los índices correspondientes a las frecuencias de interés
index_vlf = find(frequencies >= f_vlf(1) & frequencies <= f_vlf(2));
index_lf = find(frequencies >= f_lf(1) & frequencies <= f_lf(2));
index_hf = find(frequencies >= f_hf(1) & frequencies <= f_hf(2));

% Extraer las amplitudes de la PSD para cada banda de frecuencia
psd_vlf = pxx(index_vlf);
psd_lf = pxx(index_lf);
psd_hf = pxx(index_hf);

% Encontrar el índice del pico más alto en cada banda de frecuencia
[~, index_max_vlf] = max(psd_vlf);
[~, index_max_lf] = max(psd_lf);
[~, index_max_hf] = max(psd_hf);

% Obtener la frecuencia del pico más alto en cada banda
peak_frequency_vlf = frequencies(index_vlf(index_max_vlf));
peak_frequency_lf = frequencies(index_lf(index_max_lf));
peak_frequency_hf = frequencies(index_hf(index_max_hf));

% Calcular la potencia en cada banda de frecuencia en ms^2
power_vlf = sum(psd_vlf) * (frequencies(2) - frequencies(1));
power_lf = sum(psd_lf) * (frequencies(2) - frequencies(1));
power_hf = sum(psd_hf) * (frequencies(2) - frequencies(1));

% Calcular la potencia en cada banda de frecuencia en ms^2
power_vlf_n = trapz(psd_vlf);
power_lf_n= trapz(psd_lf);
power_hf_n = trapz(psd_hf);

% Calcular el porcentaje de potencia en cada banda de frecuencia
total_power = trapz(pxx);
percent_power_vlf = (power_vlf_n / total_power) * 100;
percent_power_lf = (power_lf_n/ total_power) * 100;
percent_power_hf = (power_hf_n / total_power) * 100;

% Mostrar los resultados
disp(['Pico VLF: ' num2str(peak_frequency_vlf) ' Hz']);
disp(['Pico LF: ' num2str(peak_frequency_lf) ' Hz']);
disp(['Pico HF: ' num2str(peak_frequency_hf) ' Hz']);
disp(['Potencia VLF: ' num2str(power_vlf*1000000) ' ms^2']);
disp(['Potencia LF: ' num2str(power_lf*1000000) ' ms^2']);
disp(['Potencia HF: ' num2str(power_hf*1000000) ' ms^2']);
disp(['Porcentaje de potencia VLF: ' num2str(percent_power_vlf) '%']);
disp(['Porcentaje de potencia LF: ' num2str(percent_power_lf) '%']);
disp(['Porcentaje de potencia HF: ' num2str(percent_power_hf) '%']);
end