close all;
clear all;

N = 24;
Pfa = 1e-3;
alpha_so = alpha_so(Pfa,N);
SNR_dB = 0:1:30;
Iterations = 1e6;

Reference_cells = zeros(N,Iterations);
for i = 1:N
    I = randn(1,Iterations);
    Q = randn(1,Iterations);
    Reference_cells(i,:) = (I + 1j*Q)/sqrt(2);
end

Reference_cells_AD = abs(Reference_cells).^2;
Window_Front = Reference_cells_AD(1:N/2,:);
Window_Back = Reference_cells_AD(N/2+1:end,:);
	
Sum_Reference_cells = min(sum(Window_Front),sum(Window_Back));
g = Sum_Reference_cells;

T = g.*alpha_so;

I_test = randn(Iterations,1);
Q_test = randn(Iterations,1);
Test_Signal = (I_test + 1j*Q_test)/sqrt(2);
Test_AD = abs(Test_Signal).^2;

False_Alarms = sum((Test_AD.'-T)>0);
Simulated_Pfa = False_Alarms/Iterations;

Pfa_Error = 100*(Simulated_Pfa-Pfa)/Pfa;

Pd = [];
Pd_sim = [];
for index = 1:length(SNR_dB)
    SNR_linear = 10^(SNR_dB(index)/10); 
    I_test = randn(Iterations,1);
    Q_test = randn(Iterations,1);
    Test_Signal = sqrt(SNR_linear)*(I_test + 1j*Q_test)/sqrt(2)+ Reference_cells(10,:).';
    Test_AD = abs(Test_Signal).^2;
    
    Window_Front = Reference_cells_AD(1:N/2,:);
	Window_Back = Reference_cells_AD(N/2+1:end,:);
	
	Sum_Reference_cells = min(sum(Window_Front),sum(Window_Back));
	g = Sum_Reference_cells;
	T = g.*alpha_so;
   
    Detections = ( Test_AD.'-T )>0;
    Pd_sim = [Pd_sim; sum(Detections(:))/Iterations];
    
    sum1 = 0;
    for k = 0:1:(N/2-1)
        sum1 = sum1 + nchoosek((N/2 -1 +k),k)*(2+alpha_so/(1+SNR_linear))^(-k);
    end
    Pd = [Pd ; 2*(2+alpha_so/(1+SNR_linear))^(-N/2)* sum1];

end
plot(SNR_dB,Pd.');
hold on;
plot(SNR_dB,Pd_sim.');
title(strcat('Pd vs SNR for SO-CA-CFAR with Pfa of :',num2str(Pfa)));
xlabel('SNR [dB]');
ylabel('Pd');
legend('Expected Pd','Simulated Pd')
