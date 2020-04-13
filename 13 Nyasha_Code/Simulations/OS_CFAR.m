clear all;
close all;

Pfa = 1e-3;
Iterations = 1e6;
N = 24;
SNR_dB_vector = 0:1:30;
index = 3*N/4; % kth value

alpha_values = 1:1:100;
Pfa_values = factorial(N)*factorial(alpha_values+N-index)./(factorial(N-index)*factorial(alpha_values+N));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*Pfa));
alpha_os = alpha_values(ind);

Reference_cells = zeros(N,Iterations); 

for i = 1:N
    I = randn(1,Iterations);        % Noise matrix
    Q = randn(1,Iterations);
    Reference_cells(i,:) = (I+j.*Q)*(1/sqrt(2));      
end


%Reference_cells = randn(N,Iterations)./sqrt(2);
Reference_cells_square_law = abs(Reference_cells).^2;
sorted_ref = sort(Reference_cells_square_law,1);
x_k = sorted_ref(index,:);
T_oscfar = x_k .* alpha_os;

Pd_cacfar = [];
I_test = randn(1,Iterations);
Q_test = randn(1,Iterations);
Data_noise = (I_test + j*Q_test)./sqrt(2);
Data_noise_square_detector = abs(Data_noise).^2;
Pd_calculated = [];


I_t = randn(1,Iterations);
Q_t = randn(1,Iterations);

Number_false_alarms = length(find(Data_noise_square_detector>T_oscfar));
Pfa_simulated_cacfar = Number_false_alarms/Iterations;

for SNR_dB = SNR_dB_vector
    
    SNR_linear = 10^(SNR_dB/10);
    %Calculating Pd
    Target_voltage = sqrt(SNR_linear);
    Target_signal = Target_voltage*(I_t + j*Q_t)./sqrt(2);
    Total = Reference_cells(10,:) + Target_signal;
    Total_signal_square_detector = abs(Total).^2;
    Number_detections = length(find(Total_signal_square_detector>T_oscfar));
    P_detection = Number_detections/Iterations;
    
    Pd_cacfar = [Pd_cacfar;P_detection];
    
    sum1 = 1;
    for i = 0:1:(index-1)
        sum1 = sum1 *(N-i)/(N-i+(alpha_os/(1+SNR_linear))); % Pd equation
    end
    Pd_calculated = [Pd_calculated, sum1];
                                 
end
plot(SNR_dB_vector,Pd_cacfar);
hold on;
plot(SNR_dB_vector,Pd_calculated);
xlabel('SNR [dB]');
ylabel('Pd');
title(strcat('OS-CFAR Pd vs SNR for Pfa of : ',num2str(Pfa)));
legend('Simulated Pd','Expected Pd')
