classdef RadarRecordAndProcessApp_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        UIAxes                         matlab.ui.control.UIAxes
        RECORDINGFEATUREPanel          matlab.ui.container.Panel
        STARTRECORDINGButton           matlab.ui.control.Button
        RecordinglengthEditFieldLabel  matlab.ui.control.Label
        RecordinglengthEditField       matlab.ui.control.NumericEditField
        RecordingSuccessfulCheckBox    matlab.ui.control.CheckBox
        FilenameEditFieldLabel         matlab.ui.control.Label
        FilenameEditField              matlab.ui.control.EditField
        sLabel                         matlab.ui.control.Label
        Lamp                           matlab.ui.control.Lamp
        PROCESSINGFEATUREPanel         matlab.ui.container.Panel
        WaveFilesFoundListBoxLabel     matlab.ui.control.Label
        WaveFilesFoundListBox          matlab.ui.control.ListBox
        FolderSelectedNameTextAreaLabel  matlab.ui.control.Label
        FolderSelectedNameTextArea     matlab.ui.control.TextArea
        SelectFolderButton             matlab.ui.control.Button
        RampTimeEditFieldLabel         matlab.ui.control.Label
        RampTimeEditField              matlab.ui.control.NumericEditField
        BandwidthEditFieldLabel        matlab.ui.control.Label
        BandwidthEditField             matlab.ui.control.NumericEditField
        msLabel                        matlab.ui.control.Label
        MHzLabel                       matlab.ui.control.Label
    end



    methods (Access = private)

        % Callback function
        function PROCESSButtonPushed(app, event)
     
        end

        % Callback function
        function ModeSelectionButtonGroupSelectionChanged(app, event)
        
        end

        % Button pushed function: SelectFolderButton
        function SelectFolderButtonPushed(app, event)
           global arr
           global X
           
           % uncheck box when new recording made
           app.RecordingSuccessfulCheckBox.Value=0;
           %get the folder which contains wave files
           folder = uigetdir([]); 
           %extract only wave files if folder contains other unwanted files
           d = dir([folder, '/*.wav']);
           n=length(d);
           %create a string array with number of elements equal to number of wave files
           X=strings([1,n]);
           %fill the array
           for k= 1:n
               X(:,k)=d(k).name;
           end
           
           %set items in listbox to array of waves files
           app.WaveFilesFoundListBox.Items=X;
           arr=n; % orginal array length variable used in recording feature
           app.FolderSelectedNameTextArea.Value=folder; % show folder name that was selected
        end

        % Callback function
        function NEWPROCESSButtonValueChanged(app, event)
       
        end

        % Callback function
        function NEWPROCESSButtonPushed(app, event)
          
        end

        % Value changed function: WaveFilesFoundListBox
        function WaveFilesFoundListBoxValueChanged(app, event)
% uncheck box when new recording made
app.RecordingSuccessfulCheckBox.Value=0;

value = app.WaveFilesFoundListBox.Value; %value=wave file selected

            % ---- standard DSP helper functions below ----

function [w] = hann_window(N)
% create a hann (cosine squared) window
w = .5 + .5*cos(2*pi*((0:N-1).'/(N-1) - .5));
end    
function [y] = fft_interp(x,M)
% perform approximate bandlimited interpolation of x by a factor of M
L = 4;
winInds = (-L*M : L*M).'/M * pi;

% get the ideal antialiasing filter's impulse response of length 2*M + 1 
winInds(L*M + 1) = 1;
myWin = sin(winInds) ./ winInds;
myWin(L*M + 1) = 1;

% use the window method; apply a hann window
myWin = myWin .* hann_window(2*L*M + 1);

% insert zeros in data and apply antialias filter via FFT
nFFT = numel(x) * M;
if isreal(x)
    y = ifft( fft(myWin,nFFT) .* repmat(fft(x),[M 1]), 'symmetric');
else
    y = ifft( fft(myWin,nFFT) .* repmat(fft(x),[M 1]) );
end
y = y([L*M+1:end 1:L*M]);
end      


% intelligent way of knowing if CW or FMCW mode
                    %read file (want to read it once only)
                    [Y,Fs] = audioread(value);
                    %extract channel 1 data
                    Channel1= -Y(:,1);
                    %obtain max value
                    maxValue= abs(max(Channel1));
                    %obtain mena value
                    meanValue=abs(mean(Channel1));
                    %create control variable to check which mode 
                    control=(maxValue-meanValue);
                             if control>0.1 
                                FMCWmode=true; %file is FMCW
                             else
                                FMCWmode=false;%file is CW
                             end     
                             
                             
%THIS IS THE FMCW CODE
      if FMCWmode
          
%% Input parameters by user
%Yunus comment: Get from measuring actual transmit signal
%fStart = 2400e6 (Hz) LFM start frequency 
%fStop = 2480e6 (Hz) LFM stop frequency
 
%% Setup constants and parameters
c = 299e6; % (m/s) speed of light

Tp = (app.RampTimeEditField.Value)*1e-3; % 10ms is the default :(s) minimum pulse length
            % Yunus comment: Used to provide a better estimate the pulse length 

numPad = 64; % number of samples to pad for bandlimited interpolation & shifting
             % Yunus comment: Number of extra samples to take around the trigger rising edge signal
             
ovsTrg = 16; % oversampling factor applied when interpolating the trigger signal

ovsRng = 2; % oversampling factor applied when interpolating the range data
            % Yunus comment: Oversampling the IF signal by a factor of two -
            % Yunus comment:  1. Why is oversampling needed? To reduce straddle loss?
            % Yunus comment:  2. Why a factor of two?

nPulseCancel = 3; % number of pulses to use for canceller 
                  % Yunus comment: two pulse canceller
            
maxRange = 100; % (m) maximum range to display


% ----- end constants and parameters -----

%% Read the raw wave data
[Y,Fs] = audioread(value,'native'); %FMCW uses a datype to read file

%% Derived parameters
Np = round(Tp * Fs); % Number of samples in the chirp signal (or per pulse)
BW = (app.BandwidthEditField.Value)*1e6; %fStop - fStart i.e. 80MHz is the default:(Hz) transmit bandwidth
delta_r = c/(2*BW);  % (m) range resolution

%% Change the signs of the measured parameters to be more accurate
% change sign of the input because it appears inverted in practice
trig = -Y(:,1); % the trigger signal is in the first channel
                % on rising edge, the chirp signal increases 
s = -Y(:,2); % the raw mixer output is in the second channel
clear Y;

%% Estimate the actual chirp pulse length from the measured data and store in variable Tp 
% Estimate the index of the rising edge of the trigger and store in variable pulseStarts
% parse the trigger signal (look for threshold crossings)

pulseTrig = (trig > 0);
pulseSum = sum(pulseTrig(1:Np));
pulseStarts = [];
for ii = Np+1:numel(trig)-Np-1        
    if (pulseTrig(ii) && pulseSum==0)
        pulseStarts = [pulseStarts; ii]; %#ok<*AGROW>
    end
   
    % update the running sum
    pulseSum = pulseSum + pulseTrig(ii) - pulseTrig(ii - Np);
end
clear pulseTrig;

% refine using measured parameters: the pulse width
Np = round(min(diff(pulseStarts))/2); % Np is the number of samples of the chirp pulse
                                      % pulseStarts is the index of the rising edge in samples
                                      % Want to only keep the rising edge of triangular wave
                                      % So pulsewidth is estimated by pulsestarts/2*1/fs [in seconds]
Tp = Np / Fs;

%% Pre-compute some windows and other vectors 
Nrange = floor(ovsRng*Np/2); % number of output range samples
dataRange = (0:Nrange-1).' * (delta_r/ovsRng); % labelling the range axes in meters
dataRange = dataRange(dataRange <= maxRange); % labelling the range axes in meters upto maxRange in meters
Nrange_keep = numel(dataRange); % number of range bins which are less than maxRange in meters
%
% Setup windows to be used later
rngWin = hann_window(Np); % the window applied to reduce range sidelobes
padWin = sin( (1:numPad).'/(numPad+1) * pi/2) .^2; % the window applied to the padded data
                                    % Yunus comment: why this custom window. Why not hamming, taylor, hanning, etc
trgWin = hann_window(numPad*2+1); % the window applied to the trigger data

%  Obtain the number of pulses in the data
nSamples = numel(s);
pulseStarts = pulseStarts(pulseStarts+Np+numPad <= nSamples); % pulseStarts is the index of the rising edge in samples
numPulses = numel(pulseStarts); 

% process pulses into a data matrix
sif = zeros(Nrange_keep,numPulses); % sif - Range lines in matrix form 
                                                                        
for pIdx = 1:numPulses
    % bandlimited interpolate the trigger signal
    % Yunus comment: Estimate the fraction of a bin offset needed to align pulseStarts to the zero-crossing

    tmp = double(trig(pulseStarts(pIdx) + (-numPad:numPad))) .* trgWin; 
                  % Yunus comment: pulseStarts(pIdx) contains estimate of indx of zero crossing
                  % Yunus comment: need to estimate fraction offset needed to make it exactly the zero crossing
                  % Yunus comment: need this offset to align the received data                  
    interpTmp = fft_interp(tmp,ovsTrg);
    interpTmp = interpTmp( (numPad*ovsTrg + 1) + (-ovsTrg:ovsTrg) );
    interpOffs = (-ovsTrg:ovsTrg)/ovsTrg;
    myIdx = find(diff(sign(interpTmp))==2)+1;
    tmp2 = interpTmp( myIdx + (-1:0) );
    % linear interpolate to find the zero crossing
    fracOffset = -(interpOffs(myIdx) - tmp2(2)/(tmp2(2)-tmp2(1)) / ovsTrg);
    
    % time-align the data to the trigger event (the zero crossing) 
    % Apply non-integer time-shift in the frequency domain by multiplying by a phase ramp signal
    cInds = pulseStarts(pIdx) + (-numPad:(Np+numPad-1));
    tmp = double(s(cInds)); % Yunus comment: tmp = data received after chirp pulse
                            % Yunus comment: get samples from data from (rising edge-NumPad: NumSamplesChirpPulse + NumPad) 
    tmp(1:numPad) = tmp(1:numPad) .* padWin;
    tmp(end:-1:(end-numPad+1)) = tmp(end:-1:(end-numPad+1)) .* padWin;
   
    % time delay applied in the frequency domain below
    tmp = fft(tmp);
    tmp = tmp .* exp( -1j*(0:(Np+2*numPad-1)).'/(Np+2*numPad)*2*pi*fracOffset );
    tmp = ifft(tmp,'symmetric');
    
    % compute & scale range data from the time-aligned mixer output
    tmp = ifft(tmp(numPad + (1:Np)) .* rngWin, 2*Nrange); % Need to do an IFFT to get range line
    sif(:,pIdx) = tmp(1:Nrange_keep); % sif - Range lines in matrix form
                                      % sif - each column is a range profile
end
%
clear s trig;
%
sif = sif.';

% apply the N-pulse canceller
mti_filter = -ones(nPulseCancel,1)/nPulseCancel;
midIdx = round((nPulseCancel+1)/2);
mti_filter(midIdx) = mti_filter(midIdx) + 1;
sif = convn(sif,mti_filter,'same');

% apply the median CFAR normalization
sif_dB = 20*log10(abs(sif));
% Normalising or making the maximum value 0dB
Max_sif_dB = max(max(sif_dB)); 
sif_dB = sif_dB - Max_sif_dB; 

% plot the MTI normalized results
imagesc(app.UIAxes,dataRange,(1:numPulses)*Tp*2,sif_dB);
app.UIAxes.XLabel.String= 'Range (m)';
app.UIAxes.YLabel.String= 'Time (s)';
app.UIAxes.Title.String='RTI with MTI FMCW Results'; 
colormap(app.UIAxes,jet(256)); 
caxis(app.UIAxes,[-50 0]); %creates suitable colour axes 
colorbar(app.UIAxes); 
axis(app.UIAxes, 'xy'); %flips the axes into correct orientation
axis(app.UIAxes, 'tight'); %allows plot to expand the whole axis area


        
      else 

            %THIS IS THE CW CODE
                    % ---- setup constants and parameters ----
                    c = 299e6; % (m/s) speed of light
                    cpi = 0.10; % (s) coherent processing interval
                    fc = 2590e6; % (Hz) Center frequency (connect VCO Vtune to +5)
                    maxSpeed = 30; % (m/s) maximum speed to display
                    overlapFactor = 8; % (unitless) number of overlapped pulse windows (1 for no overlap)
                    ovsDop = 4; % (unitless) oversample factor for Doppler axis
                    % ----- end constants and parameters -----
                    
                    % already read the raw wave data
                    
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
                    dti = zeros(nSpeed, M);
                    for mIdx = 1:M
                        xInds = (1:N).' + (mIdx-1)*floor(N/overlapFactor); % compute indices
                        tmp = double(x(xInds)) .* dopWin; % apply Doppler window
                        tmp = tmp - mean(tmp); % remove DC component if it exists
                        tmp = fft(tmp, ovsDop*N); % compute oversampled Doppler spectrum
                        dti(:,mIdx) = 20*log10( abs(tmp(1:nSpeed)) ); % grab result in dB
                        
                    end
                    clear x;
                    dti = dti.'; %transpose matrix
                    
                    % make Doppler vs. time plot
                    imagesc(app.UIAxes,speed,time,dti);
                    imagesc(app.UIAxes,time, speed, dti');
                    colormap(app.UIAxes, jet(256));
                    caxis(app.UIAxes,max(dti(:)) + [-60 0]); % show 60 dB dynamic range
                    colorbar(app.UIAxes);
                    app.UIAxes.XLabel.String= 'Time (s)';
                    app.UIAxes.YLabel.String= 'Speed (m/s)';
                    app.UIAxes.Title.String='Doppler CW Results';
                    axis(app.UIAxes, 'xy');%flips the axes into correct orientation
                    axis(app.UIAxes, 'tight'); %allows plot to expand the whole axis area
                    
      end    

        end

        % Value changed function: FolderSelectedNameTextArea
        function FolderSelectedNameTextAreaValueChanged(app, event)
   
        end

        % Button pushed function: STARTRECORDINGButton
        function STARTRECORDINGButtonPushed(app, event)
global arr
global X
            % uncheck box when new recording made
            app.RecordingSuccessfulCheckBox.Value=0;
            
            %set parameters for recording
            Fs=44100;
            NumChannels=2;
            nBits=16;
            ID=-1; %default value
            recorder = audiorecorder(Fs,nBits,NumChannels,ID); %creates recording object
            length=app.RecordinglengthEditField.Value;
            app.Lamp.Color='r'; %change lamp colour to red to show busy recording
            recordblocking(recorder, length); %record for length of time specified without interruptions
            app.Lamp.Color='g'; %change lamp colour to green to show ready to record again and end of last recording
            
            %the progress bar could only displayed after recording done- can't seem to implement progress bar and record at same time
 
            output=getaudiodata(recorder);  %gets array of data from channel 1 and channel 2 
            
            %save the recorded file as a wav file with the filename entered by user
            filename= app.FilenameEditField.Value;
            filename=strcat(filename,'.wav');
            audiowrite(filename,output,Fs);
            
            app.RecordingSuccessfulCheckBox.Value=1; %show the recording was successful
            
           %%update/refresh folder contents after new recording was made
           arr=arr+1;%update index if multiple recoridngs made
           X(:,arr)=filename;
           %set items in listbox to array of waves files
           app.WaveFilesFoundListBox.Items=X;
          
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'UI Figure';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'RADAR RESULTS')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.FontName = 'Gill Sans MT';
            app.UIAxes.FontSize = 18;
            app.UIAxes.Position = [233 11 398 466];

            % Create RECORDINGFEATUREPanel
            app.RECORDINGFEATUREPanel = uipanel(app.UIFigure);
            app.RECORDINGFEATUREPanel.TitlePosition = 'centertop';
            app.RECORDINGFEATUREPanel.Title = 'RECORDING FEATURE';
            app.RECORDINGFEATUREPanel.FontName = 'Bahnschrift';
            app.RECORDINGFEATUREPanel.FontWeight = 'bold';
            app.RECORDINGFEATUREPanel.FontSize = 14;
            app.RECORDINGFEATUREPanel.Position = [17 11 208 188];

            % Create STARTRECORDINGButton
            app.STARTRECORDINGButton = uibutton(app.RECORDINGFEATUREPanel, 'push');
            app.STARTRECORDINGButton.ButtonPushedFcn = createCallbackFcn(app, @STARTRECORDINGButtonPushed, true);
            app.STARTRECORDINGButton.FontName = 'Bahnschrift';
            app.STARTRECORDINGButton.FontSize = 14;
            app.STARTRECORDINGButton.Position = [16 52 143 25];
            app.STARTRECORDINGButton.Text = 'START RECORDING';

            % Create RecordinglengthEditFieldLabel
            app.RecordinglengthEditFieldLabel = uilabel(app.RECORDINGFEATUREPanel);
            app.RecordinglengthEditFieldLabel.HorizontalAlignment = 'right';
            app.RecordinglengthEditFieldLabel.FontName = 'Bahnschrift';
            app.RecordinglengthEditFieldLabel.FontSize = 14;
            app.RecordinglengthEditFieldLabel.Position = [20 94 113 22];
            app.RecordinglengthEditFieldLabel.Text = 'Recording length ';

            % Create RecordinglengthEditField
            app.RecordinglengthEditField = uieditfield(app.RECORDINGFEATUREPanel, 'numeric');
            app.RecordinglengthEditField.HorizontalAlignment = 'center';
            app.RecordinglengthEditField.FontName = 'Bahnschrift';
            app.RecordinglengthEditField.FontSize = 14;
            app.RecordinglengthEditField.Position = [140 94 29 22];

            % Create RecordingSuccessfulCheckBox
            app.RecordingSuccessfulCheckBox = uicheckbox(app.RECORDINGFEATUREPanel);
            app.RecordingSuccessfulCheckBox.Text = 'Recording Successful?';
            app.RecordingSuccessfulCheckBox.FontName = 'Bahnschrift';
            app.RecordingSuccessfulCheckBox.FontSize = 14;
            app.RecordingSuccessfulCheckBox.Position = [22 16 162 22];

            % Create FilenameEditFieldLabel
            app.FilenameEditFieldLabel = uilabel(app.RECORDINGFEATUREPanel);
            app.FilenameEditFieldLabel.HorizontalAlignment = 'right';
            app.FilenameEditFieldLabel.FontName = 'Bahnschrift';
            app.FilenameEditFieldLabel.FontSize = 14;
            app.FilenameEditFieldLabel.Position = [8 130 63 22];
            app.FilenameEditFieldLabel.Text = 'Filename';

            % Create FilenameEditField
            app.FilenameEditField = uieditfield(app.RECORDINGFEATUREPanel, 'text');
            app.FilenameEditField.HorizontalAlignment = 'center';
            app.FilenameEditField.FontName = 'Bahnschrift';
            app.FilenameEditField.FontSize = 14;
            app.FilenameEditField.Position = [79 130 118 22];

            % Create sLabel
            app.sLabel = uilabel(app.RECORDINGFEATUREPanel);
            app.sLabel.FontName = 'Bahnschrift';
            app.sLabel.FontSize = 14;
            app.sLabel.Position = [171 94 25 22];
            app.sLabel.Text = 's';

            % Create Lamp
            app.Lamp = uilamp(app.RECORDINGFEATUREPanel);
            app.Lamp.Position = [168 54 20 20];

            % Create PROCESSINGFEATUREPanel
            app.PROCESSINGFEATUREPanel = uipanel(app.UIFigure);
            app.PROCESSINGFEATUREPanel.TitlePosition = 'centertop';
            app.PROCESSINGFEATUREPanel.Title = 'PROCESSING FEATURE';
            app.PROCESSINGFEATUREPanel.FontName = 'Candara';
            app.PROCESSINGFEATUREPanel.FontWeight = 'bold';
            app.PROCESSINGFEATUREPanel.FontSize = 15;
            app.PROCESSINGFEATUREPanel.Position = [16 207 207 270];

            % Create WaveFilesFoundListBoxLabel
            app.WaveFilesFoundListBoxLabel = uilabel(app.PROCESSINGFEATUREPanel);
            app.WaveFilesFoundListBoxLabel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.WaveFilesFoundListBoxLabel.HorizontalAlignment = 'center';
            app.WaveFilesFoundListBoxLabel.FontName = 'Pontano Sans';
            app.WaveFilesFoundListBoxLabel.FontSize = 14;
            app.WaveFilesFoundListBoxLabel.FontWeight = 'bold';
            app.WaveFilesFoundListBoxLabel.Position = [35 63 112 22];
            app.WaveFilesFoundListBoxLabel.Text = 'Wave Files Found';

            % Create WaveFilesFoundListBox
            app.WaveFilesFoundListBox = uilistbox(app.PROCESSINGFEATUREPanel);
            app.WaveFilesFoundListBox.Items = {'wave file 1', 'wave file 2', 'wave file 3', 'wave file 4'};
            app.WaveFilesFoundListBox.ValueChangedFcn = createCallbackFcn(app, @WaveFilesFoundListBoxValueChanged, true);
            app.WaveFilesFoundListBox.FontName = 'Pontano Sans';
            app.WaveFilesFoundListBox.FontSize = 14;
            app.WaveFilesFoundListBox.BackgroundColor = [0.9412 0.9412 0.9412];
            app.WaveFilesFoundListBox.Position = [21 7 168 55];
            app.WaveFilesFoundListBox.Value = 'wave file 1';

            % Create FolderSelectedNameTextAreaLabel
            app.FolderSelectedNameTextAreaLabel = uilabel(app.PROCESSINGFEATUREPanel);
            app.FolderSelectedNameTextAreaLabel.HorizontalAlignment = 'center';
            app.FolderSelectedNameTextAreaLabel.FontName = 'Bahnschrift';
            app.FolderSelectedNameTextAreaLabel.FontSize = 14;
            app.FolderSelectedNameTextAreaLabel.Position = [28 110 145 22];
            app.FolderSelectedNameTextAreaLabel.Text = 'Folder Selected Name:';

            % Create FolderSelectedNameTextArea
            app.FolderSelectedNameTextArea = uitextarea(app.PROCESSINGFEATUREPanel);
            app.FolderSelectedNameTextArea.ValueChangedFcn = createCallbackFcn(app, @FolderSelectedNameTextAreaValueChanged, true);
            app.FolderSelectedNameTextArea.HorizontalAlignment = 'center';
            app.FolderSelectedNameTextArea.FontName = 'Bahnschrift';
            app.FolderSelectedNameTextArea.FontSize = 14;
            app.FolderSelectedNameTextArea.Position = [42 92 121 19];

            % Create SelectFolderButton
            app.SelectFolderButton = uibutton(app.PROCESSINGFEATUREPanel, 'push');
            app.SelectFolderButton.ButtonPushedFcn = createCallbackFcn(app, @SelectFolderButtonPushed, true);
            app.SelectFolderButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.SelectFolderButton.FontName = 'OCR A Extended';
            app.SelectFolderButton.FontSize = 14;
            app.SelectFolderButton.FontWeight = 'bold';
            app.SelectFolderButton.Position = [33 140 134 30];
            app.SelectFolderButton.Text = 'Select Folder';

            % Create RampTimeEditFieldLabel
            app.RampTimeEditFieldLabel = uilabel(app.PROCESSINGFEATUREPanel);
            app.RampTimeEditFieldLabel.HorizontalAlignment = 'right';
            app.RampTimeEditFieldLabel.FontName = 'Courier New';
            app.RampTimeEditFieldLabel.FontSize = 14;
            app.RampTimeEditFieldLabel.Position = [39 182 81 22];
            app.RampTimeEditFieldLabel.Text = 'Ramp Time';

            % Create RampTimeEditField
            app.RampTimeEditField = uieditfield(app.PROCESSINGFEATUREPanel, 'numeric');
            app.RampTimeEditField.HorizontalAlignment = 'center';
            app.RampTimeEditField.FontName = 'Courier New';
            app.RampTimeEditField.FontSize = 14;
            app.RampTimeEditField.Position = [130 182 30 22];

            % Create BandwidthEditFieldLabel
            app.BandwidthEditFieldLabel = uilabel(app.PROCESSINGFEATUREPanel);
            app.BandwidthEditFieldLabel.HorizontalAlignment = 'center';
            app.BandwidthEditFieldLabel.FontName = 'Courier New';
            app.BandwidthEditFieldLabel.FontSize = 14;
            app.BandwidthEditFieldLabel.Position = [34 214 81 22];
            app.BandwidthEditFieldLabel.Text = 'Bandwidth';

            % Create BandwidthEditField
            app.BandwidthEditField = uieditfield(app.PROCESSINGFEATUREPanel, 'numeric');
            app.BandwidthEditField.HorizontalAlignment = 'center';
            app.BandwidthEditField.FontName = 'Courier New';
            app.BandwidthEditField.FontSize = 14;
            app.BandwidthEditField.Position = [119 214 32 22];

            % Create msLabel
            app.msLabel = uilabel(app.PROCESSINGFEATUREPanel);
            app.msLabel.FontName = 'Courier New';
            app.msLabel.FontSize = 14;
            app.msLabel.Position = [164 182 25 22];
            app.msLabel.Text = 'ms';

            % Create MHzLabel
            app.MHzLabel = uilabel(app.PROCESSINGFEATUREPanel);
            app.MHzLabel.FontName = 'Courier New';
            app.MHzLabel.FontSize = 14;
            app.MHzLabel.Position = [154 214 31 22];
            app.MHzLabel.Text = 'MHz';
        end
    end

    methods (Access = public)

        % Construct app
        function app = RadarRecordAndProcessApp_exported

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end