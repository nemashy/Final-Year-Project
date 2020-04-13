close all;
clear all;

N = 24;
Pfa = 1e-3;
alpha_go = alpha_go(Pfa,N);
SNR_dB = 0:1:30;

Iterations = 1e6; % number of columns
Reference_cells = zeros(N,Iterations);

for i = 1:N
    I = randn(1,Iterations);
    Q = randn(1,Iterations);
    Reference_cells(i,:) = (I + 1j*Q)/sqrt(2);
end
 
Reference_cells_AD = abs(Reference_cells).^2; % using a square law detector
Window_Front = Reference_cells_AD(1:N/2,:);
Window_Back = Reference_cells_AD(N/2+1:end,:);	
Sum_Reference_cells = max(sum(Window_Front),sum(Window_Back));
g = Sum_Reference_cells;
T = g.*alpha_go; % threshold values

I_test = randn(Iterations,1);
Q_test = randn(Iterations,1);
noise = (I_test + 1j*Q_test)/sqrt(2);
noise_AD = abs(noise).^2; % noise

False_Alarms = sum((noise_AD.'-T)>0);  % detection decision
Simulated_Pfa = False_Alarms/Iterations;

Pfa_Error = 100*(Simulated_Pfa-Pfa)/Pfa;

Pd = [];
Pd_sim = [];
for index = 1:length(SNR_dB)
    SNR_linear = 10^(SNR_dB(index)/10);
    
    I_test = randn(Iterations,1);
    Q_test = randn(Iterations,1);
    noise_target = sqrt(SNR_linear)*(I_test + 1j*Q_test)/sqrt(2)+ Reference_cells(10,:).'; %noise + target returns
    noise_target_AD = abs(noise_target).^2;
    
    
    Window_Front = Reference_cells_AD(1:N/2,:);
	Window_Back = Reference_cells_AD(N/2+1:end,:);
	
	Sum_Reference_cells = max(sum(Window_Front),sum(Window_Back));
	g = Sum_Reference_cells;
	T = g.*alpha_go;
    
    
    Detections = ( noise_target_AD.'-T)>0;
    
    Pd_sim = [Pd_sim; sum(Detections(:))/Iterations];
    
    
    % Math for expected values
    
    sum1 = 0;
    for k = 0:1:(N/2-1)
        sum1 = sum1 + nchoosek((N/2 -1 +k),k)*(2+alpha_go/(1+SNR_linear))^(-k);
    end
    
    
    Pd = [Pd ; 2*((1+alpha_go/(1+SNR_linear))^(-N/2) - ((2+ alpha_go/(1+SNR_linear))^(-N/2) * sum1))];% Pd equation

end

plot(SNR_dB,Pd.');
hold on;
plot(SNR_dB,Pd_sim.');
title(strcat('Pd vs SNR for GO-CA-CFAR with Pfa of :',num2str(Pfa)));
xlabel('SNR [dB]');
ylabel('Pd');
legend('Expected Pd','Simulated Pd')
