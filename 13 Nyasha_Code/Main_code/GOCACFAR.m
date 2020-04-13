function [Detections] = GOCACFAR(Pfa,window_size,g_cells,sif1)

sif = sif1.';
B = sif(:);

Detections = [];
target = 0;
N = window_size *2;
alpha = GO_CA_CFAR_Statistic(Pfa,N);

signal = abs(B).^2;

columns = length(signal);

% applying a cell averaging cfar
for q = window_size+g_cells+1:columns-window_size-g_cells
    cut_power = signal(q);
    window_back = signal(q-window_size-(g_cells)/2:q-(g_cells/2)-1);
    window_front = signal(q+(g_cells/2)+1:q+(g_cells/2)+window_size);
    sum_reference_cells = max(sum(window_front),sum(window_back));
    g = sum_reference_cells;
    threshold = alpha*g; 
    if (threshold<cut_power)
         Detections = [Detections; q];
         target = target + 1;
    end
end
end
%close all;