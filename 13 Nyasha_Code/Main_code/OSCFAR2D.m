function [Detections] = OSCFAR2D(Pfa,window_size,g_cells,sif1)

target = 0;
N = 2*(window_size);
index = 3*N/4;
Detections = [];

alpha_values = 1:1:500;
Pfa_values = factorial(N)*factorial(alpha_values+N-index)./(factorial(N-index)*factorial(alpha_values+N));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*Pfa));
alpha = alpha_values(ind);

signal = abs(sif1).^2;
[rows, columns]= size(sif1);
for d = 1:rows
    for q = window_size+g_cells/2+1:columns-window_size-g_cells/2
        cut_power = signal(d,q);
        window_back = [];
        window_front = [];
        for wind = q-window_size-(g_cells)/2:q-(g_cells/2)-1
            window_back = [window_back; signal(d,wind)];
        end
        for wind1 = q+(g_cells/2)+1:q+(g_cells/2)+window_size
             window_front = [window_front; signal(d,wind1)];
        end
        reference_cells = [window_back; window_front];
        reference_cells_sorted = sort(reference_cells);
        x_k = reference_cells_sorted(index);
        threshold = alpha*x_k; 
        if (threshold<cut_power)
             Detections =  [Detections; q+(d-1)*columns];
             target = target + 1;
        end
    end
end
%close all;