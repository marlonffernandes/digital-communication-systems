%##########################################################################
%                   Comunicação Digital
%                   Exercício de Simulação Nº 5
%               Marlon Ferraz Fernandes 
%##########################################################################

clear all;
close all;

% Parâmetros da codificação BCH(63, 57, 1)
n = 63;
k = 57;
t = 1;
gp = bchgenpoly(n, k);
bchEncoder = comm.BCHEncoder(n, k, gp, t);
bchDecoder = comm.BCHDecoder(n, k, gp, t);

[X, Y] = meshgrid(0:k, 0:n-k);
T = (( (X+Y) > t ) + 0) .* X;

tamanhoDados = 1e6;
i = 1;
for EbN0 = 0:1:8    % Conforme Tabela 1: Eb/N0 variando de 0 a 8 dB
    
    % Variáveis auxiliares
    EbN0_linear = 10^(EbN0/10);
    sigma = sqrt(10^(-EbN0/10)/2);  % Calcula a variância do ruído a ser gerado         
    dados = randi([0 1], tamanhoDados, 1);   % Gera os dados a serem transmitidos 
    coded_raw_ber(i) = 0.5*erfc(sqrt( (k/n) * EbN0_linear ));    
    BER_AWGN_teorica_cod(i) = 0*coded_raw_ber(i);

    %----------------------------------------------------------------------  
    % Simulação do desempenho de BPSK no canal AWGN
    
    % Calcula a BER Teórica sem codificação
    BER_AWGN_teorica(i) = 0.5*erfc(sqrt(10^(EbN0/10)));
    
    % Calcula a BER Teórica com codificação
    msg_err = binopdf(0:k, k, coded_raw_ber(i));
    par_err = binopdf(0:(n-k), n-k, coded_raw_ber(i))';
    joint = msg_err(ones(1, n-k+1), :) .* par_err(:, ones(1, k+1));
    BER_AWGN_teorica_cod(i) = sum( sum(T .* joint) ) / k;
    
    % Calcula a BER Simulada
    total_erros = 0;    % Reinicia a contagem de erros de bits
    total_amostras = 0; % Reinicia a contagem de bits transmitidos
    while total_erros < 100 && total_amostras < 10*tamanhoDados
        encodedData = bchEncoder(dados);            % BCH encode data
        modSignal = pskmod(encodedData, 2);         % BPSK modulate
        receivedSignal = awgn(modSignal, (k/n)*EbN0);     % Pass through AWGN channel
        demodSignal = pskdemod(receivedSignal, 2);  % BSPK demodulate
        receivedBits = bchDecoder(demodSignal);     % BCH decode data
        num_erros = biterr(dados, receivedBits);
        
        total_erros = total_erros + num_erros;
        total_amostras = total_amostras + length(dados);
    end
    BER_AWGN_simulada(i) = total_erros/total_amostras;
    
    %----------------------------------------------------------------------
    % Simulação do desempenho BPSK em canais com desvanecimento
    
    % Calcula a BER Teórica sem codificação
    a = ( gamma(1.5)) / (2*sqrt(pi)*gamma(2) );
    b = (1+EbN0_linear).^(-1);
    c = hypergeom([1, 0.5], 2, (1/(1+EbN0_linear)));
    BER_Rayl_teorica(i) = a*b*c;
    
    % Calcula a BER Simulada com codificação
    total_erros = 0;    % Reinicia a contagem de erros de bits
    total_amostras = 0; % Reinicia a contagem de bits transmitidos
    while total_erros < 100 && total_amostras < 10*tamanhoDados
        encodedData = bchEncoder(dados);            % BCH encode data
        modSignal = pskmod(encodedData, 2);         % BPSK modulate
        
        fading = random('Nakagami', 1, 1, [length(encodedData), 1]); % Gera fading Rayleigh
        fadingSignal = fading .* modSignal;                          % TX*fading
        receivedSignal = normrnd(fadingSignal, sigma);               % RX = TX*fading + ruído
        
        demodSignal = pskdemod(receivedSignal, 2);  % BSPK demodulate
        receivedBits = bchDecoder(demodSignal);     % BCH decode data
        
        num_erros = biterr(dados, receivedBits);
        total_erros = total_erros + num_erros;
        total_amostras = total_amostras + length(dados);
    end
    BER_Rayl_simulada(i) = total_erros/total_amostras;
    
    SNR(i) = EbN0;
    i = i+1;
end

% Plota as curvas para o canal AWGN
figure('NumberTitle', 'off', 'Name', 'AWGN');
semilogy(SNR, BER_AWGN_teorica(1,:), 'b-', 'DisplayName', 'BER Teórica sem Código');
hold on;
semilogy(SNR, BER_AWGN_teorica_cod(1,:), 'r-', 'DisplayName', 'BER Teórica com Código');
hold on;
semilogy(SNR, BER_AWGN_simulada(1,:), 'g-', 'DisplayName', 'BER Simulada com Código');
hold on;
grid on;
ylabel('BER'); 
xlabel('Eb/N0 (dB)');
legend;
title('Desempenho BPSK com Cod. BCH(63,57,1) - AWGN');
set(gca, 'YLim', [1e-7 1]);

% Plota as curvas para o canal Rayleigh
figure('NumberTitle', 'off', 'Name', 'Rayleigh');
semilogy(SNR, BER_Rayl_teorica(1,:), 'b-', 'DisplayName', 'BER Teórica sem Código');
hold on;
semilogy(SNR, BER_Rayl_simulada(1,:), 'g-', 'DisplayName', 'BER Simulada com Código');
hold on;
grid on;
ylabel('BER'); 
xlabel('Eb/N0 (dB)');
legend;
title('Desempenho BPSK com Cod. BCH(63,57,1) - Rayleigh');
set(gca, 'YLim', [1e-7 1]);
