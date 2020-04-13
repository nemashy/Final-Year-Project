clear all;
close all;

Pfa = 1e-3;
Iterations = 1e6;
N = 24;
SNR_dB_vector = 0:1:30;

alpha = N*(Pfa^(-1/N)-1);


Ref_cells = zeros(N,Iterations); 

for i = 1:N
    I = randn(1,Iterations);        % Noise matrix
    Q = randn(1,Iterations);
    Ref_cells(i,:) = (I+j.*Q)*(1/sqrt(2));      % Linear Detector
end


%Reference_cells = randn(N,Iterations)./sqrt(2);
Reference_cells = Ref_cells;
Reference_cells_square_law = abs(Reference_cells).^2;
Sum_column = sum(Reference_cells_square_law,1);
Avg = Sum_column./N;
T_cacfar = Avg .* alpha;

Pd_cacfar = [];
I_test = randn(1,Iterations);
Q_test = randn(1,Iterations);
Data_noise = (I_test + j*Q_test)./sqrt(2);
Data_noise_square_detector = abs(Data_noise).^2;
Pd_calculated = [];


I_t = randn(1,Iterations);
Q_t = randn(1,Iterations);

Number_false_alarms = length(find(Data_noise_square_detector>T_cacfar));
Pfa_simulated_cacfar = Number_false_alarms/Iterations;

for SNR_dB = SNR_dB_vector
    
    SNR_linear = 10^(SNR_dB/10);
    %Calculating Pd
    Target_voltage = sqrt(SNR_linear);
    Target_signal = Target_voltage*(I_t + j*Q_t)./sqrt(2);
    Total = Reference_cells(10,:) + Target_signal;
    Total_signal_square_detector = abs(Total).^2;
    Number_detections = length(find(Total_signal_square_detector>T_cacfar));
    P_detection = Number_detections/Iterations;
    
    Pd_cacfar = [Pd_cacfar;P_detection];
     
    Pd_calculated = [Pd_calculated, 1/((1+alpha/(N*(1+SNR_linear)))^N)];
                                 
end
plot(SNR_dB_vector,Pd_cacfar);
hold on;
plot(SNR_dB_vector,Pd_calculated);


