%% Paper
%Simulation Study of Double Threshold Energy
%Detection Method for Cognitive Radios
%Pankaj Verma ?, Brahmjit Singh
%
%
format long

N = [500,1000];  
snr_dB =-9; %dB
snr = 10.^(snr_dB./10);
Pf = 1/100:1/100:1;
%-------------------------------------

%% BPSK Signal

L=1500;
data = round(rand(1,L));                               % Data sequence
uni2bip=2*data-1;                                      % Convert unipolar to bipolar
T=1;                                                   % Bit duration
Eb=T/2;                                                % This will result in unit amplitude waveforms
fc=3/T;                                                % Carrier frequency
t=linspace(0,5,1500);                                  % discrete time sequence between 0 and 5*T (15000 samples)
K=length(t);                                           % Number of samples
Nsb=K/length(data);                                    % Number of samples per bit
dd=repmat(data',1,Nsb);                                % replicate each bit Nsb times
bb=repmat(uni2bip',1,Nsb); dw=dd';                     % Transpose the rows and columns
dw=dw(:)'; 

%------ Convert dw to a column vector (colum by column) and convert to a row vector
bw=bb';
bw=bw(:)';                                             % Data sequence samples
w=sqrt(2*Eb/T)*cos(2*pi*fc*t);                         % carrier waveform
bpsk_w=bw.*w;                                          % modulated waveform

%-----AWGN noise with mean 0 and variance -----%
        Noise1 = randn(1,N(1)); 
        Noise2 = randn(1,N(2));
        vn1=var(Noise1);
        vn2=var(Noise2);
        
        %-----Real valued Gaussian Primary User Signal------% 
         Signal1 = sqrt(snr).*bpsk_w(1:500); 
         Signal2 = sqrt(snr).*bpsk_w(1:1000);

         vs1=var(Signal1);
         vs2=var(Signal2);
 
%------- Threshold-----------

Threshold_1 = N(1)*vn1 + qfuncinv(Pf)*sqrt(2*N(1)*(vn1)^2);
Threshold_2 = N(2)*vn2 + qfuncinv(Pf)*sqrt(2*N(2)*(vn2)^2);
%------------------------------------
%% Probabilty of detection theory
Pd_the1 = qfunc((Threshold_1 -N(1)*(vn1+vs1))./(sqrt(2*N(1)*(vn1+vs1)^2)));
Pd_the2 = qfunc((Threshold_2 -N(2)*(vn2+vs2))./(sqrt(2*N(2)*(vn2+vs2)^2)));


%% Probabilty of detection simulated

hwait = waitbar(0,'Please wait ....');
for i=1:length(Pf)
    D1=0;
    D2=0;
    for j=1:10000    %Monte Carlo simulation
        
        %-----AWGN noise with mean 0 and variance 1-----%
        Noise1 = randn(1,N(1)); 
        Noise2 = randn(1,N(2));
        v_n1=var(Noise1);
        v_n2=var(Noise2);
        
        %-----BPSK Signal ------% 
        
        Signal1 = sqrt(snr).*bpsk_w(1:500); 
        Signal2 = sqrt(snr).*bpsk_w(1:1000);
        v_s1=var(Signal1);
        v_s2=var(Signal2);
        
        Recv_Sig1 = Signal1 + Noise1; % Received signal at SU 1
        Recv_Sig2 = Signal2 + Noise2; % Received signal at SU 2
        
        Energy1 = abs(Recv_Sig1).^2; % Energy of received signal over N samples
        Energy2 = abs(Recv_Sig2).^2;
        
        %-----Computation of Test statistic for energy detection-----%
        Test_Statistic1 =sum(Energy1);
        Test_Statistic2 =sum(Energy2);
        
        %-----Theoretical value of Threshold-----%
        Threshold1(i) = N(1)*v_n1 + qfuncinv(Pf(i))*sqrt(2*N(1)*(v_n1)^2);
        Threshold2(i) = N(2)*v_n2 + qfuncinv(Pf(i))*sqrt(2*N(2)*(v_n2)^2);
        
        if(Test_Statistic1>= Threshold1(i))  % Check whether the received energy is greater than threshold, if so,(Probability of detection) counter by 1
           D1 = D1+1;
           
        end
        if(Test_Statistic2 >= Threshold2(i))  % Check whether the received energy is greater than threshold, if so,(Probability of detection) counter by 1
           D2 = D2+1;
            
        end
    end

      Pd_1(i)=D1/j;
      Pd_2(i)=D2/j;
      waitbar(i/length(Pf),hwait);
end
close(hwait);
plot(Pf,Pd_the1,'b')
grid on
hold on
plot(Pf,Pd_the2,'r--s')
plot(Pf,Pd_1,'g o')
plot(Pf,Pd_2,'c *')
axis([0.0001,1,0.0001,1]);

legend('Theoritical, N=500','Stimulated, N=500',...
    'Theoritical, N=1000','Stimulated, N=1000')










