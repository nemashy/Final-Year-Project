clear all;
close all;

N = 24; 
Pfa = 10^-3;        
Iterations = 1e6;   
alpha_ca = N*((Pfa)^(-1/N)-1);
Reference_cells = zeros(N,Iterations); 

% Noise matrix
for i = 1:N
    I = randn(1,Iterations);        
    Q = randn(1,Iterations);
    Reference_cells(i,:) = (I+j.*Q)*(1/sqrt(2));   
end

sld = abs(Reference_cells).^2;  % Signal After Square Law Detector 
Sum_Ref_cells = sum(sld,1);

g = Sum_Ref_cells/N; 

% CA Threshold
Threshold = g.*alpha_ca;

I_test = randn(Iterations,1);
Q_test = randn(Iterations,1);
Test_Signal = Reference_cells(1,:)/sqrt(2); 
Test_AD = abs(Test_Signal).^2;

False_Alarms = sum((Test_AD-Threshold)>0);
Sim_Pfa = False_Alarms/Iterations;

Pfa_Error = 100*(Sim_Pfa-Pfa)/Pfa;

% Calculate Pd

SNR_dB_Range = 0:1:30;
Pd_Simulated = [];
Pd_Calculated = [];

I = randn(1,Iterations);
Q = randn(1,Iterations);

for SNR_dB = SNR_dB_Range
    SNR_linear = 10^(SNR_dB/10);
    cut = sqrt(SNR_linear)*(I+ j*Q)*(1/sqrt(2)) + Reference_cells(5,:);
    cut_AD = abs(cut).^2;
    Detections = sum((cut_AD-Threshold)>0);
    Pd_Simulated = [Pd_Simulated, Detections/Iterations];
    Pd_Calculated = [Pd_Calculated, 1/((1+alpha_ca/(N*(1+SNR_linear)))^N)]; %Pd equation
end 

plot(SNR_dB_Range,Pd_Simulated);
hold on;
plot(SNR_dB_Range,Pd_Calculated);
xlabel('SNR [dB]');
ylabel('Pd');
title(strcat('CA-CFAR Pd vs SNR for Pfa of : ',num2str(Pfa)));
legend('Simulated Pd','Expected Pd')