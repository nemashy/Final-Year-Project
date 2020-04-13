%% Processing recordings

clear all;
close all;
 

FMCW_mode = 1; % = 1 : FMCW (Range and Doppler)
               % = 0: CW    (Doppler only)
               
RecordingNo2Process = 2;
 
wavFile_CW_All = {'20130603 01 test 1 cw fahad to door.wav'; 
                  '20130603 02 test 2 cw fahad walking from door to pole.wav';
                  '20130603 03 test 3 cw fahad walking to door yunus running to and back.wav'};
    
wavFile_FMCW = '20130602_ASED_lab_test1'; 

wavFile_FMCW_All = { '20130603 04 test 1 fmcw fahad arms swinging.wav';
                     '20130603 05 test 2 fmcw fahad walking away yunus running towards.wav';
                     '20130603 06 test 3 fmcw yunus running away fahad walking away.wav'};
    % FMCW files
    % '20130603 04 test 1 fmcw fahad arms swinging'
    % '20130603 05 test 2 fmcw fahad walking away yunus running towards'
    % '20130603 06 test 3 fmcw yunus running away fahad walking away'
               
 if FMCW_mode == 1
     % FMCW             
        wavFile_FMCW = wavFile_FMCW_All{RecordingNo2Process};
        cantenna_rti_v3_yunus(wavFile_FMCW);
 else
     % CW
       wavFile_CW = wavFile_CW_All{RecordingNo2Process};
       cantenna_dop_v3_yunus(wavFile_CW);
 end

