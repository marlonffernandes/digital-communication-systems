close all;

% Simula��o do desempenho de BPSK no canal AWGN
i=1;
for ebno = 0:1:10
    numErros = 0;
    sigma = sqrt(10^(-ebno/10)/2);           
    amostras = 2500000;
    RX = normrnd(-1, sigma, 1, amostras);    % Consideramos a transmiss�o de bits '0' com energia unit�ria (s�mbolos -1) 
    TX_est = RX>0;                           % Decisor no receptor               
    numErros = sum(TX_est);                  % Contagem de erros de bit
    ber_teorica(i) = 0.5*erfc(sqrt(10^(ebno/10)));
    ber_simulada(i) = numErros/amostras;
    snr(i) = ebno;
    i = i+1;
end

% Plota as curvas da BER Te�rica e Simulada no canal AWGN
semilogy(snr, ber_teorica(1,:),'DisplayName','BER Te�rica - AWGN');
hold on;
semilogy(snr, ber_simulada(1,:),'*','DisplayName','BER Simulada - AWGN');
hold on;

% Simula��o do desempenho BPSK em canais com desvanecimento
m = [1];   % Valores de m para a distribui��o de Nakagami-m
for n = 1:length(m)
    i = 1;
    for ebno = 0:1:10
        % Calcula a BER Te�rica
        ebno_linear = 10^(ebno/10);
        a = (gamma(m(n)+0.5))/(2*sqrt(pi)*gamma(m(n)+1));
        b = (1+(ebno_linear/m(n))).^(-m(n));
        c = hypergeom([m(n),0.5],(m(n)+1),(1/(1+(ebno_linear/m(n)))));
        ber_teorica(i) = a*b*c;
        
        % Calcula a BER Simulada
        numErros = 0;                   % Reinicia a contagem de erros
        amostras = 2500000;             % Consideramos a transmiss�o de bits '0'
        novas_amostras = [1, amostras];
        sigma = sqrt(10^(-ebno/10)/2);  % Calcula a vari�ncia do ru�do a ser gerado         
        
        % Gera v�rias amostras de h (fading)
        fading = random('Nakagami', m(n), 1, novas_amostras);
        RX = normrnd(-fading, sigma);       % RX = TX*fading + ruido
        TX_est = (RX./fading) > 0;          % Estima o valor de x dividindo pelo 
        
        numErros = sum(TX_est);             % Contagem de erros de bit
        ber_simulada(i) = numErros/amostras;
        snr(i) = ebno;
          
        i = i+1;
    end
    
    % Plota as curvas da BER Te�rica e Simulada considerando desvanecimento
    % dado pela distribui��o de Nakagami-m e Rayleigh (m = 1)
    if m(n) == 1
        legendaTeorica = strcat('BER Te�rica - Rayleigh');
        legendaSimulada = strcat('BER Simulada - Rayleigh');
    else
        legendaTeorica = strcat('BER Te�rica - Nakagami-m (m=', num2str(m(n)), ')');
        legendaSimulada = strcat('BER Simulada - Nakagami-m (m=', num2str(m(n)), ')');
    end
    semilogy(snr, ber_teorica(1,:), 'DisplayName', legendaTeorica);
    hold on;
    semilogy(snr, ber_simulada(1,:), '*', 'DisplayName', legendaSimulada);
    hold on;
end

grid on;
ylabel('BER'); 
xlabel('Eb/N0 (dB)');
legend;
title('Desempenho da Modula��o BPSK em Canais com Desvanecimento');
