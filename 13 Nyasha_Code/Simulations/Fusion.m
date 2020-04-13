close all;
clear all;

N = 24;
Pfa = 1e-4;
SNR_dB = 0:1:30;


Iterations = 1e6;

Reference_cells = zeros(N,Iterations);

for i = 1:N
    I = randn(1,Iterations);
    Q = randn(1,Iterations);
    Reference_cells(i,:) = (I + 1j*Q)/sqrt(2);
end

Reference_cells_AD = abs(Reference_cells).^2;

%% GO-CA-CFAR
alpha_go = GO_CA_CFAR_Statistic(Pfa,N);
Window_Front = Reference_cells_AD(1:N/2,:);
Window_Back = Reference_cells_AD(N/2+1:end,:);
	
Sum_Reference_cells = max(sum(Window_Front),sum(Window_Back));
g = Sum_Reference_cells;
T_go = g.*alpha_go;
%% OS-CFAR
index = 3*N/4;
alpha_values = 1:1:500;
Pfa_values = factorial(N)*factorial(alpha_values+N-index)./(factorial(N-index)*factorial(alpha_values+N));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*Pfa));
alpha_os = alpha_values(ind);

sorted_ref = sort(Reference_cells_AD,1);
x_k = sorted_ref(index,:);
T_os = x_k .* alpha_os;

%% CA
alpha_ca = N*(Pfa^(-1/N)-1);
Sum_column = sum(Reference_cells_AD,1);
Avg = Sum_column./N;
T_ca = Avg * alpha_ca;

I_test = randn(Iterations,1);
Q_test = randn(Iterations,1);
Test_Signal = (I_test + 1j*Q_test)/sqrt(2);
Test_AD = abs(Test_Signal).^2;
Test_AD = Test_AD.';


%% Finding positions where there are detections
detection_go = find((Test_AD - T_go)>0);
detection_os = find((Test_AD - T_os)>0);
detection_ca = find((Test_AD - T_ca)>0);

% binary distiction
A_go = zeros(Iterations,1);
A_os = zeros(Iterations,1);
A_ca = zeros(Iterations,1);

for i = 1:length(detection_go)
    A_go(detection_go(i)) = 1;
end
for i = 1:length(detection_os)
    A_os(detection_os(i)) = 1;
end
for i = 1:length(detection_ca)
    A_ca(detection_ca(i)) = 1;
end

% adding the binary value
A = A_go + A_os + A_ca;
%%
False_Alarms = length(find(A>1));
Simulated_Pfa = False_Alarms/Iterations;
Pfa_Error = 100*(Simulated_Pfa-Pfa)/Pfa;
%%
Pd = [];
Pd_sim = [];
for index = 1:length(SNR_dB)
    SNR = 10^(SNR_dB(index)/10);
    %Math for simulated values
    
    I_test = randn(Iterations,1);
    Q_test = randn(Iterations,1);
    Test_Signal = sqrt(SNR)*(I_test + 1j*Q_test)/sqrt(2)+ Reference_cells(10,:).';
    Test_AD = abs(Test_Signal).^2;
    Test_AD = Test_AD.';
    
    %% GO-CA-CFAR
    alpha_go = GO_CA_CFAR_Statistic(Pfa,N);
    Window_Front = Reference_cells_AD(1:N/2,:);
    Window_Back = Reference_cells_AD(N/2+1:end,:);

    Sum_Reference_cells = max(sum(Window_Front),sum(Window_Back));
    g = Sum_Reference_cells;
    T_go = g.*alpha_go;
    %% OS-CFAR
    index = 3*N/4;
    alpha_values = 1:1:500;
    Pfa_values = factorial(N)*factorial(alpha_values+N-index)./(factorial(N-index)*factorial(alpha_values+N));
    [val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*Pfa));
    alpha_os = alpha_values(ind);

    sorted_ref = sort(Reference_cells_AD,1);
    x_k = sorted_ref(index,:);
    T_os = x_k .* alpha_os;

    %% CA
    alpha_ca = N*(Pfa^(-1/N)-1);
    Sum_column = sum(Reference_cells_AD,1);
    Avg = Sum_column./N;
    T_ca = Avg * alpha_ca;

    %% Finding positions where there are detections
    detection_go = find((Test_AD - T_go)>0);
    detection_os = find((Test_AD - T_os)>0);
    detection_ca = find((Test_AD - T_ca)>0);

    % binary distiction
    A_go = zeros(Iterations,1);
    A_os = zeros(Iterations,1);
    A_ca = zeros(Iterations,1);

    for i = 1:length(detection_go)
        A_go(detection_go(i)) = 1;
    end
    for i = 1:length(detection_os)
        A_os(detection_os(i)) = 1;
    end
    for i = 1:length(detection_ca)
        A_ca(detection_ca(i)) = 1;
    end

    % adding the binary value
    A = A_go + A_os + A_ca;
    Detections = length(find(A>1));
    Pd_sim = [Pd_sim; Detections/Iterations];

end

plot(SNR_dB,Pd_sim.');
title(strcat('Pd vs SNR for Fusion_CFAR with Pfa of :',num2str(Pfa)));
xlabel('SNR [dB]');
ylabel('Pd');
hold on
