clear all;
close all;
File = 'Dataset4.wav';
Pfa = 1e-8;
window_size = 8;
g_cells = 18;

[sif1,dataRange,numPulses,Tp] = cantenna_rti_v3_yunus(File);

Detections_1 = OSCFAR(Pfa,window_size,g_cells,sif1);
Detections_2 = CACFAR(Pfa,window_size,g_cells,sif1);
Detections_3 = GOCACFAR(Pfa,window_size,g_cells,sif1);
Detections_4 = SOCACFAR(Pfa,window_size,g_cells,sif1);
%Detections_6 = CACFAR2D(Pfa,window_size,g_cells,sif1);
%Detections_7 = GOCACFAR2D(Pfa,window_size,g_cells,sif1);
%Detections_8 = OSCFAR2D(Pfa,window_size,g_cells,sif1);
%Detections_9 = SOCACFAR2D(Pfa,window_size,g_cells,sif1);
% Setup constants and parameters

sif = abs(sif1).^2;

hold on;

%locate position of cell on the image and place an 

[num_rows, num_columns] = size(sif);
%choose the detection algorithm to use
Detections = Detections_4;


for k = Detections
    t = floor(k/num_columns)+1;
    time = t*Tp*2; 
    position_range_vector = mod(k,num_columns)+1;
end

%discrading unwanted noise cells and take region of interest in terms of range cells

%{
noise = find((position_range_vector<8)|(position_range_vector>44));
for h = noise
    position_range_vector(h) = [];
    time(h) = [];
end
noise = find((position_range_vector<12)&(time<7.5));
for h = noise
    position_range_vector(h) = [];
    time(h) = [];
end
noise = find((position_range_vector<12)&(time>12));
for h = noise
    position_range_vector(h) = [];
    time(h) = [];
end
noise = find((time<8)&(position_range_vector<15));
for h = noise
    position_range_vector(h) = [];
    time(h) = [];
end
%}

%Defining the time limit for region of interest
min_time = 0; 
%min_time = 1;
%max_time = 7.5;
max_time = numPulses*Tp*2;
t_value_lower = round(min_time/(Tp*2));
t_value_upper = round(max_time/(Tp*2));
t_values = t_value_lower:t_value_upper;

% find available times
times = unique(time,'rows');
detected_columns = length(find(round(times/(Tp*2))==t_values));
P_detection = detected_columns/length(t_values);

%print 'x' on detected position
text(dataRange(position_range_vector),time,'X','FontSize',5);

hold off;
