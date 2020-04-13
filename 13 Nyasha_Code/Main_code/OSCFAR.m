function [Detections] = OSCFAR(Pfa,window_size,g_cells,sif1)

sif = sif1.';
B = sif(:);
target = 0;
signal = abs(B).^2;
N = 2*(window_size);
index = 3*N/4;
Detections = [];

alpha_values = 1:1:500;
Pfa_values = factorial(N)*factorial(alpha_values+N-index)./(factorial(N-index)*factorial(alpha_values+N));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*Pfa));
alpha = alpha_values(ind);

columns = length(signal);

for q = window_size+g_cells/2+1:columns-window_size-g_cells/2
    cut_power = signal(q);
    window_back = signal(q-window_size-(g_cells)/2:q-(g_cells/2)-1);
    window_front = signal(q+(g_cells/2)+1:q+(g_cells/2)+window_size);
    reference_cells = [window_back; window_front];
    reference_cells_sorted = sort(reference_cells);
    x_k = reference_cells_sorted(index);
    threshold = alpha*x_k; 
    if (threshold<cut_power)
         Detections = [Detections; q];
         target = target + 1;
    end
end
%close all;