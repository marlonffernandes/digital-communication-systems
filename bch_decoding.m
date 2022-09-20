n = 63;
k = 57;
t = 1;
rate = k/n;

ebno = 0:0.5:10;

ebno_lin = 10.^(ebno/10);
bpsk_ber = 0.5*erfc(sqrt(ebno_lin));
coded_raw_ber = 0.5*erfc(sqrt(rate*ebno_lin));
coded_ber = 0*coded_raw_ber;

% Theory
[X,Y] = meshgrid(0:k,0:n-k);
T = (((X+Y)>t)+0).*X;
for i=1:length(ebno)
  msg_err = binopdf(0:k,k,coded_raw_ber(i));
  par_err = binopdf(0:(n-k),n-k,coded_raw_ber(i))';
  joint = msg_err(ones(1,n-k+1),:) .* par_err(:,ones(1,k+1));
  coded_ber(i) = sum(sum(T.*joint))/k;
end

% Simulation (all zeros codeword)
M = 500;
sim_ber = zeros(1,length(ebno));
for i=1:2:15
  y = -1 + randn(M,n)/sqrt(2*rate*ebno_lin(i));
  ybit = (y>0)+0;
  xhat = bchdec(gf(ybit,1),n,k);
  sim_ber(i) = sum(sum(xhat~=0))/M/k;
end

semilogy(ebno,bpsk_ber,'r-',ebno,coded_raw_ber,'g-',ebno,coded_ber,'b-',ebno,sim_ber,'bo');
title('(63,36) BCH with bounded distance decoding');
legend('BPSK','Input to decoder','After decoder','Simulation',3);
xlabel('E_n / N_0');
ylabel('BER');
grid on;
set(gca,'YLim',[1e-7 1]);