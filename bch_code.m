N = 63;
K = 57;
S = 15;
gp = bchgenpoly(63,57);
bchEncoder = comm.BCHEncoder(N,K,gp,S);
bchDecoder = comm.BCHDecoder(N,K,gp,S);
errorRate = comm.ErrorRate('ComputationDelay',3); 
for counter = 1:20
  data = randi([0 0],630,1);                 % Generate binary data
  encodedData = bchEncoder(data);            % BCH encode data
  modSignal = pskmod(encodedData,2);         % BPSK modulate
  receivedSignal = awgn(modSignal,5);        % Pass through AWGN channel
  demodSignal = pskdemod(receivedSignal,2);  % BSPK demodulate
  receivedBits = bchDecoder(demodSignal);    % BCH decode data
  errorStats = errorRate(data,receivedBits); % Compute error statistics
end
fprintf('Error rate = %f\nNumber of errors = %d\n', ...
  errorStats(1), errorStats(2))