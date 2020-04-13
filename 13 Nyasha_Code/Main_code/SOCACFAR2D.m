function [Detections] = SOCACFAR2D(Pfa,window_size,g_cells,sif1)

Detections = [];
target = 0;
N = window_size *2;
alpha = SO_CA_CFAR_Statistic(Pfa,N);

signal = abs(sif1).^2;

[rows, columns]= size(sif1);

% applying a cell averaging cfar
for d =1:rows
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
        sum_reference_cells = min(sum(window_front),sum(window_back));
        g = sum_reference_cells;
        threshold = alpha*g; 
        if (threshold<cut_power)
             Detections = [Detections; q+(d-1)*columns];
             target = target + 1;
        end
    end
end
end
%close all;