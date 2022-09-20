% -------------------------------------------------------------------------
% Simulacao de um Sistema de Transmissao Digital em banda Base
% Sera transmitida uma imagem usando modulacao 2-PAM
% -------------------------------------------------------------------------

clear; clc; close all;

% -------------------------------------------------------------------------
% Par?metros
% -------------------------------------------------------------------------

Fs=44100;           % Frequencia de Amostragem do Sinal Continuo (Hz)
Ts=1/Fs;            % Periodo de Amostragem (s)

oversampling=10;    % Fator Fs/R

R=Fs/oversampling;  % Taxa de Transmissao em simbolos/s (baud rate)
T=1/R;              % Periodo de Simbolo (s)

del=25;             % Resposta do filtro formatador se estende por (2*del) periodos de simbolo
                    % Numero de amostras do filtro formatador: N=2*(del*oversampling)+1

rolloff=0.5;        % Fator de rolloff dos filtros Tx e Rx

% -------------------------------------------------------------------------
% Leitura da Imagem e mapeamento 2-PAM
%--------------------------------------------------------------------------

im_in=imread('shuttle_80x60.tif');  %Leitura da imagem a ser transmitida
%im_in=imread('lenna512.tif');

L=8;
[size_r,size_c]=size(im_in);
im_size=size_r*size_c;
im_vec=reshape(im_in,1,im_size);

bit_matrix=de2bi(im_vec);
bit_per_symbol=2;               %constelaco 2^bit_per_symbol-PAM, neste caso 2-PAM
bit_symbols=reshape(bit_matrix, im_size*L/bit_per_symbol, bit_per_symbol);
symbols=bi2de(bit_symbols);     %sequencia de bits a ser transmitida
alfabeto=[-3 -1 1 3];                       %alfabeto 2-PAM
symbols=symbols+1;
pam=alfabeto(symbols);


% -------------------------------------------------------------------------
% Filtro de Tx+Rx
% Formatacao de Pulso - Tipo Raiz Quadrada do Cosseno Levantado no Tx e Rx         
% -------------------------------------------------------------------------
filtro=rcosfir(rolloff,del,oversampling,1,'sqrt');

%Formatar e transmitir os simbolos PAM

sinal_tx=upsample(pam,oversampling);           % Realiza Upsampling
sinal_tx_filtrado=conv(sinal_tx,filtro);       % Sinal Filtrado de Transmissao

% -------------------------------------------------------------------------
% Canal Passa Baixas Simulado usando um Filtro de Butterworth
% -------------------------------------------------------------------------

fc=R/2*(1+rolloff);                 % Largura de banda do sinal transmitido
[bn,an]=butter(5,2*fc/Fs);          % Filtro passa-baixas de largura de banda fc

sinal_rx=filtfilt(bn,an,sinal_tx_filtrado);

sigPower = sum(sinal_rx(:).^2)/length(sinal_rx);
sigPower = 10*log10(sigPower);
reqSNR = 5;                         %Relacao sinal-ruido (Eb/N0) em dB
noisePower = sigPower-reqSNR+10*log10(oversampling/2);

wn = wgn(1,length(sinal_rx),noisePower);

sinal_rx_ruido=sinal_rx+wn;

figure(1);

subplot(2,2,1);                    

% Densidade Espectral de Potencia do Sinal Transmitido
[Pxx,F]=pwelch(sinal_tx_filtrado/max(abs(sinal_tx_filtrado)),[],[],[],Fs,'twosided');
plot((F-Fs/2)/1000,10*log10(fftshift(Pxx))); 
grid; xlabel('Frequ?ncia (kHz)'); ylabel('dBm/Hz');
xlim([-(Fs/4)/1000 (Fs/4)/1000]);
title('Sinal Transmitido')

subplot(2,2,2);                    

% Densidade Espectral de Potencia do sinal recebido
[Pxx,F]=pwelch(sinal_rx/max(abs(sinal_rx)),[],[],[],Fs,'twosided');
plot((F-Fs/2)/1000,10*log10(fftshift(Pxx))); 
grid; xlabel('Frequ?ncia (kHz)'); ylabel('dBm/Hz');
xlim([-(Fs/4)/1000 (Fs/4)/1000]);
title('Sinal na Sa?da do Canal')

% -------------------------------------------------------------------------
% Receptor (Filtro Casado)
% -------------------------------------------------------------------------

sinal_rx_casado=conv(sinal_rx_ruido,filtro);       %Filtro casado

subplot(2,2,3);
[H,F]=freqz(bn,an,2048,'whole',Fs);
gain=20*log10(fftshift(abs(H)));
plot((F-Fs/2)/1000,gain); grid;
axis([-Fs/(4*1000) Fs/(4*1000) -50 10]);   %Grafico do m?dulo da funcao de transferencia
xlabel('Frequ?ncia (kHz)');                %do filtro que simula o canal
ylabel('Ganho (dB)');
title('Resposta de Amplitude do Canal');


subplot(2,2,4);                    

% Densidade Espectral de Potencia do sinal recebido
[Pxx,F]=pwelch(sinal_rx_casado/max(abs(sinal_rx)),[],[],[],Fs,'twosided');
plot((F-Fs/2)/1000,10*log10(fftshift(Pxx))); 
grid; xlabel('Frequ?ncia (kHz)'); ylabel('dBm/Hz');
xlim([-(Fs/4)/1000 (Fs/4)/1000]);
title('Sinal na Sa?da do Filtro Casado')

pam_rx=downsample(sinal_rx_casado,oversampling);     
pam_rx=pam_rx(del*2+1:length(pam_rx)-del*2);        

%Estimacao dos simbolos PAM recebidos

symbols_rx_quant=quantalph(pam_rx,alfabeto);

bit_symbols_rx=(symbols_rx_quant+1)/2;

bit_matrix_rx=reshape(bit_symbols_rx, im_size, L);   %Reconstrucao da imagem no receptor
im_vec_rx=bi2de(bit_matrix_rx);
im_in_rx=reshape(im_vec_rx,size_r,size_c);     

figure(2);

subplot(2,1,1);                   %Visualizacao da imagem transmitida
colormap(gray);
h=image(im_in);
set(h,'CDataMapping','scaled')
axis('equal');
title('Imagem Transmitida');
hold;

subplot(2,1,2);                   %Visualizacaoo da imagem recebida
colormap(gray);
h=image(im_in_rx);
set(h,'CDataMapping','scaled')
axis('equal');
title('Imagem Recebida');

eyediagram(sinal_rx_casado(1,1:5000),2*oversampling)     %Diagrama de olho

[num_erros,BER]=symerr(bit_matrix_rx,bit_matrix);

disp ('N?mero de erros = ') 
disp(num_erros)
disp ('BER simulada = ')
disp(BER)
BER_teorica = 0.5*erfc(sqrt(10^(reqSNR/10)));
disp ('BER teorica = ')
disp(BER_teorica)
