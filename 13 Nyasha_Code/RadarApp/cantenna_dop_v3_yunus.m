function [] = cantenna_dop_v3_yunus(wavFile)
% function [] = cantenna_dop_v3(wavFile)
% 
% Produces a DTI (Doppler x time intensity) image of the
% cantenna recording. Assumes the data is collected with
% the radar in a continuous wave (CW) mode. See inside
% the function for more parameters.
%
%    wavFile = the filename of the .WAV file to process
%    
% Lincoln Laboratory MIT PEP Radar Course
% (c) 2012 Massachusetts Institute of Technology

% ---- setup constants and parameters ----
c = 299e6; % (m/s) speed of light
cpi = 0.10; % (s) coherent processing interval 
fc = 2590e6; % (Hz) Center frequency (connect VCO Vtune to +5)
maxSpeed = 30; % (m/s) maximum speed to display
overlapFactor = 8; % (unitless) number of overlapped pulse windows (1 for no overlap)
ovsDop = 4; % (unitless) oversample factor for Doppler axis
% ----- end constants and parameters -----

% use a default filename if none is given
if ~exist('wavFile','var')
    wavFile = 'radar_test2.wav';
end

% read the raw wave data
fprintf('Loading WAV file...\n');
[Y,Fs] = audioread(wavFile);

% derived parameters
N = round(cpi*Fs); % # of samples per pulse

% the input appears to be inverted
x = -Y(:,2); % Received signal at baseband
clear Y;

% grab an integer number of overlapped frames
M = floor(numel(x) / N * overlapFactor) - (overlapFactor) + 1;

% compute axes parameters for the plot
% Note: the Doppler data is oversampled by ovsDop
delta_f = (0:ovsDop*N/2-1).' / (ovsDop*N) * Fs; % Doppler freq. (Hz)
lambda = c / fc; % wavelength (m)
speed = delta_f * lambda / 2; % Doppler -> speed
time = (1:M) * cpi / overlapFactor; % collection time (sec)

% limit the speed axis to a reasonable range
speed = speed(speed <= maxSpeed);
nSpeed = numel(speed);

% compute a Doppler window
dopWin = hann_window(N);

% compute the Doppler vs. time plot
fprintf('Processing...\n');
dti = zeros(nSpeed, M);
for mIdx = 1:M
    xInds = (1:N).' + (mIdx-1)*floor(N/overlapFactor); % compute indices
    tmp = double(x(xInds)) .* dopWin; % apply Doppler window
    tmp = tmp - mean(tmp); % remove DC component if it exists
    tmp = fft(tmp, ovsDop*N); % compute oversampled Doppler spectrum
    dti(:,mIdx) = 20*log10( abs(tmp(1:nSpeed)) ); % grab result in dB
end
clear x;
dti = dti.';

% make Doppler vs. time plot
figure;
imagesc(time, speed, dti');
%imagesc(speed,time,dti);
colormap(jet(256));
caxis(max(dti(:)) + [-60 0]); % show 60 dB dynamic range
colorbar;
ylabel('Speed (m/sec)');
xlabel('Time (sec)');
axis xy;


% ---- standard DSP helper functions below ----

function [w] = hann_window(N)
% create a hann (cosine squared) window
w = .5 + .5*cos(2*pi*((0:N-1).'/(N-1) - .5));