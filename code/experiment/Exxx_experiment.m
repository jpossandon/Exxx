%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exxx_experiment.m 
% Rapid serial tactile presentation + TOJ task
% This is an EEG/Eye-tracking/Tactile-stimulation experiment
% The experiment consists in the presentation of tactile stimulus to the 
% hands in uncrossed or crossed positions, meanwhile gaze is directed to the 
% center of the screen, at the middle between the position of the hands.
% Tactile stimulation occurs in a sequence of variable length according to a
% flat hazard function. During the sequence stimulus occur randomly to the
% left or the right hand at different SOA, with a minimmum of 750 ms (to
% prevent overlap of ERPs or TFRs). At the end of the sequence, the last
% two stimuli are always to different hands and at asynchronies that permit
% the calibration of a TOJ performance curve. Subject respond which side was 
% stimulated by an eye movement to a left or right box displayed in the screen.
%
% Temporal profile
%
%     ts(l/r)--------ts(l/r)---...---ts(l/r)-----ts(!previous)-------response
% ms:       750+rnd(500)   750+rnd(500)     constant
%                                           stimulus
%                                           distribution:(20 30 50 90 170 330 650 1290 1800)
%
%     
%     ts: tactile stimulus (25 ms-200 Hz)
%     response: eye movement to left-right boxes
%
%
% If the experiment script is interrupted, run the following lines in the
% command windows so a new experiment can be initiated:
%       sca                                         % close PTB windows
%       Eyelink('Shutdown');                        % closes the link to the eye-tracker
%       fclose(obj);                                % closes the serial port
%       PsychPortAudio('Close', pahandle);          % closes the audio port
%
% This experiment script is a modification of Touch_experiment.m, also by 
% JPO and used in a previous experiment without EEG in NBP lab (Osnabr?ck)
% JPO, January-2016, Hamburg
%
% TOCHECK: 
% time delay to start audio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all                                                                 % we clear parameters?
% this is for debugging
win.DoDummyMode             = 0;                                            % (1) is for debugging without an eye-tracker, (0) is for running the experiment
win.stim_test               = 1;                                            % (1) for testing the stimulators (always when doing the experiment), (0) to skip
%  PsychDebugWindowConfiguration(0,1);                                       % this is for debugging with a single screen

% Screen parameters
win.whichScreen             = 0;                                            % (CHANGE?) here we define the screen to use for the experiment, it depend on which computer we are using and how the screens are conected so it might need to be changed if the experiment starts in the wrong screen
win.FontSZ                  = 20;                                           % font size
win.bkgcolor                = 0;                                            % screen background color, 127 gray
win.Vdst                    = 66;                                           % (!CHANGE!) viewer's distance from screen [cm]         
win.res                     = [1920 1080];                                  %  horizontal x vertical resolution [pixels]
win.wdth                    = 51;                                           %  51X28.7 cms is teh size of Samsung Syncmaster P2370 in BPN lab EEG rechts
win.hght                    = 28.7;                                         % 
win.pixxdeg                 = win.res(1)/(2*180/pi*atan(win.wdth/2/win.Vdst));% 
win.center_thershold        = 3*win.pixxdeg;                                % distance from the midline threshold for gaze contingent end of trial
win.fixpos                  = [win.res(1)/2, win.res(2)/2];
win.tgtpos                  = [win.fixpos(1,1)-win.fixpos(1,1)/2 win.fixpos(1,2);...
                                win.fixpos(1,1)+win.fixpos(1,1)/2 win.fixpos(1,2)];
win.fixcolor                = 127;
win.fixrad                  = 10;

% Tactile stimulation settings 
win.tact_freq               = 200;                                          % frequency of stimulation in Hz
win.stim_dur                = .015;                                         % duration of tactile stimulation. The vibrator takes some time to stop its motion so for around 50 ms we use 25 ms of stimulation time (ask Tobias for the exact latencies they have measured)
win.tact_minlat             = .8;                                           % ISI durign the sequence goes from tact_minlat to tact_max_lat, according to a flat hazard function
win.halflife                = .1;
win.decay                   = log(2)./win.halflife;                         % lambda
win.tact_max_lat            = 1.8;                                          % to avoid getting by chance a very large number, we set a limit on 99% of the exponential distribution [win.tact_minlat+win.halflife:win.halflife:1.8;cumsum(1./[2.^[1:10]])]

% Blocks and trials
win.total_tact              = 4500;                                         % 4500 total tact and 450 test gives 4500-2*450=3600 stimuli (900 per condition) beside the two last ones of the TOJ task 
win.total_test              = 450;
win.exp_trials              = win.total_test;
win.t_perblock              = 9;
win.calib_every             = 5; 

% Audio, white noise parameters
win.wn_vol                  = 1;                                           % (CHANGE?) adjust to subject comfort

% Device input during the experiment
win.in_dev                  = 1;                                            % (1) - keyboard  (2) - mouse  (3) - pedal (?)    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENT START
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (!CHANGE!) adjust this to the appropiate screen
display(sprintf('\n\n\n\n\n\nPlease check that the display screen resolution is set to 1920x1080 and the screen borders are OK. \nIf screen resolution is changes matlab (and the terminal) need to be re-started.\n'))
if ismac                                                                    % this bit is just so I can run the experiment in my mac without a problem
    exp_path    	= '/Users/jossando/trabajo/Exxx/';                      % path in my mac
else
    exp_path    	= '/home/th/Experiments/Exxx/';
end

win.s_n             = input('Subject number: ','s');                        % subject id number, this number is used to open the randomization file
win.fnameEDF     	= sprintf('s%02d.EDF',str2num(win.s_n));                % EDF name can be only 8 letters long, so we can have numbers only between 01 and 99
pathEDF            	= [exp_path 'data/' sprintf('s%02d/',str2num(win.s_n))];% where the EDF files are going to be saved

if exist([pathEDF win.fnameEDF],'file')                                     % checks whether there is a file with the same name
    rp = input(sprintf('!Filename %s already exist, do you want to overwrite it (y/n)?',win.fnameEDF),'s');
    if (strcmpi(rp,'n') || strcmpi(rp,'no'))
        error('filename already exist')
    end
end
if ~isdir(pathEDF)
    mkdir(pathEDF)
end
win.s_age         	= input('Subject age: ','s');
win.s_hand         	= input('Subject handedness for writing (l/r): ','s');
win.s_gender     	= input('Subject gender (m/f): ','s');
win.setStr      	= sprintf('Subject %s\nAge %s\nHandedness %s\nGender %s\n',...
                        win.s_n,win.s_age,win.s_hand,win.s_gender);         % setting summary
fprintf(win.setStr); 

AssertOpenGL();                                                             % check if Psychtoolbox is working (with OpenGL) TODO: is this needed?
ClockRandSeed();                                                            % this changes the random seed

[IsConnected]       = EyelinkInit(win.DoDummyMode);                        	% open the link with the eyetracker
assert(IsConnected==1, 'Failed to initialize EyeLink!')


prevVerbos = Screen('Preference','Verbosity', 2);                           % this two lines it to set how much we want the PTB to output in the command and display window 
prevVisDbg = Screen('Preference','VisualDebugLevel',3);                     % verbosity-1 (default 3); vdbg-2 (default 4)
Screen('Preference', 'SkipSyncTests', 2)                                    % for maximum accuracy and reliability

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START PTB SCREEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[win.hndl, win.rect]        = Screen('OpenWindow',...                       % starts PTB screen
                                win.whichScreen,win.bkgcolor);  
[win.cntr(1), win.cntr(2)]  = WindowCenter(win.hndl);                       % get where is the display screen center
Screen('BlendFunction',win.hndl, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     % enable alpha blending for smooth drawing
HideCursor(win.hndl);                                                     % this to hide the mouse
Screen('TextSize', win.hndl, win.FontSZ);                                   % sets teh font size of the text to be diplayed
KbName('UnifyKeyNames');                                                    
win.start_time = clock;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START AUDIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Audio signal control in one of the stereo channels the white noise that
% is present to the loudspeaker and in the other the sinusoidal signal that
% is presented to the tactor stimulator amplifier (this signal is
% continouous, the actual stimulation is controlled by the parallel port
% signal)

freq                        = 48000;                                        % noise frequency in Hz
dur                         = 5;                                            % noise duration in seconds (noise goes in a loop so duration is not so relevant)
wavedata1                   = rand(1,freq*dur)*2-1;                         % the noise
wavedata2                   = sin(2.*pi.*win.tact_freq.*...                 % the stimulus
                                [0:1/freq:dur-1/freq]);      
wavedata                    = [wavedata1 ; wavedata2];     
nrchannels                  = 2;                                           

InitializePsychSound;                                                       % Perform basic initialization of the sound driver:

try                                                                         % Try with the 'freq'uency we wanted:
    pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
catch                                                                       % Failed. Retry with default frequency as suggested by device:
    fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
    fprintf('Sound may sound a bit out of tune, ...\n\n');
    psychlasterror('reset');
    pahandle = PsychPortAudio('Open', [], [], 0, [], nrchannels);
end
PsychPortAudio('FillBuffer', pahandle, wavedata);                           % Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('Volume', pahandle,win.wn_vol);                              % Sets the volume (between 0 and 1)
% s = PsychPortAudio('GetStatus', pahandle);                                % Status of the port for testing


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUCTIONS IN GERMAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if win.in_dev == 1
    texts.txtdev = ['Leertaste dr' 252 'cken (press Space)'];
elseif win.in_dev == 2
    texts.txtdev = 'Maustate klicken';
end

% This is the safest way so umlaut are corretly displayed UTF. 
% code for umlauts 228(a),252(u)
texts.txt1     = double(['Die Funktionsf' 228 'higkeit der Ger' 228 'te muss '...    
             252 'berpr' 252 'ft werden\n' texts.txtdev ' zum Fortfahren..']);             
texts.txt2     = double(['Stimulatoren werden auf den Handr' 252 'cken besfestigt ...\n' texts.txtdev]);
texts.txt3     = double(['Der linke Stimulator wird getestet (vibriert drei mal), ...\n danach ' texts.txtdev]);
texts.txt4     = double(['Der rechte Stimulator wird getestet (vibriert drei mal), ...\n danach ' texts.txtdev]);
texts.txt6     = double(['Beginn des Experiments \n Test Block \n Die H' 228 ...
            'nde bitte parallel positionieren (parallel). \n Zum Beginnen die ' texts.txtdev]);
texts.txt7     = double(['Beginn des Experiments\n Test Block \n  F' 252 'r den n' 228 ...
            'chsten Block bitte die H' 228 'nde ' 252 'berkreuzen (crossed). \n Zum Fortfahren die ' texts.txtdev]);      
texts.txt10    = double(['F' 252 'r den n' 228 'chsten Block die H' 228 ...
        'nde parallel positionieren. Zum Fortfahren die ' texts.txtdev]);
texts.txt11    = double(['F' 252 'r den n' 228 'chsten Block die H' 228 ...
        'nde ' 252 'berkreuzen. Zum Fortfahren die '  texts.txtdev]);
% texts.txt12    = double(['Links oder Rechts']);
%these are for debugging
% handstr  = {'Left','Right','','','Left','Right'};
% crossstr = {'Uncrossed','Crossed'};
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EYE-TRACKER SETUP, OPEN THE EDF FILE AND ADDS INFO TO ITS HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[win.el, win.elcmds] = setup_eyetrackerxxx(win, 1);                            % Setups the eye-tracker with the before mentioned parameters

DrawFormattedText(win.hndl,texts.txt1,'center','center',255,55);                  % This 'draws' the text, nothing is displayed until Screen flip
Screen('Flip', win.hndl);                                                   % This is the command that changes the PTB display screen. Here it present the first instructions.

if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
    waitForKB_linux({'space'});                                             % press the space key in the keyboard
elseif win.in_dev == 2
    [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
end

OpenError = Eyelink('OpenFile', win.fnameEDF);                                  % opens the eye-tracking file. It can only be done after setting-up the eye-tracker 
if OpenError                                                                % error in case it is not possible, never happened that I know, but maybe if the small hard-drive aprtition of the eye=tracker is full
    error('EyeLink OpenFile failed (Error: %d)!', OpenError), 
end
Eyelink('Command', sprintf('add_file_preamble_text ''%s''', win.setStr ));       % this adds the information about subject to the end of the header  
wrect = Screen('Rect', win.hndl);                                           
Eyelink('Message','DISPLAY_COORDS %d %d %d %d', 0, 0, wrect(1), wrect(2));  % write display resolution to EDF file

ListenChar(2)                                                               % disable MATLAB windows' keyboard listen (no unwanted edits)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIMULATOR TEST
% ALWAYS: left-hand is # 1 and right-hand is # 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if win.stim_test
    DrawFormattedText(win.hndl,texts.txt2,'center','center',255,55);
    Screen('Flip', win.hndl);
    if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
        waitForKB_linux({'space'});                                             % press the space key in the keyboard
    elseif win.in_dev == 2      
        [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
    end

    DrawFormattedText(win.hndl,texts.txt3,'center','center',255,55);                  % TEST LEFT STIMULATOR (three times)
    Screen('Flip', win.hndl);
    t1              = PsychPortAudio('Start', pahandle, 0, 0, 0);               % Start the sounds

    WaitSecs(1);
    for e=1:3
        Eyelink('command', '!*write_ioport 0x378 0');                  
        WaitSecs(1+rand(1));      
        Eyelink('command', '!*write_ioport 0x378 1');                           % start stimulation by sending a signal through the parallel port (a number that was set by E275_define_tact_states)
        WaitSecs(win.stim_dur);                                                 % for the specified duration
        Eyelink('command', '!*write_ioport 0x378 0');                          % stop stimulation
        WaitSecs(win.stim_dur);
    end
    if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
        waitForKB_linux({'space'});                                             % press the space key in the keyboard
    elseif win.in_dev == 2
        [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
    end

    DrawFormattedText(win.hndl,texts.txt4,'center','center',255,55);                  % TEST RIGHT STIMULATOR (three times)
    Screen('Flip', win.hndl);
    WaitSecs(1);
    for e=1:3
        Eyelink('command', '!*write_ioport 0x378 0');
        WaitSecs(1+rand(1));      
        Eyelink('command', '!*write_ioport 0x378 2');                           % start stimulation by sending a signal through the parallel port (a number that was set by E275_define_tact_states)
        WaitSecs(win.stim_dur);       
        Eyelink('command', '!*write_ioport 0x378 0');                          % stop stimulation
        WaitSecs(win.stim_dur);
    end
    
    PsychPortAudio('Stop', pahandle);

    if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
        waitForKB_linux({'space'});                                             % press the space key in the keyboard
    elseif win.in_dev == 2
        [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
    end
    fprintf('\nTest ready!');
else
    fprintf('\nSkipping stimulators Test.\n Only during experiment debugging');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EYE-TRACKER CALIBRATION AND VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EyelinkDoTrackerSetup(win.el);                                              % calibration/validation (keyboard control)
win.el.eye_used = Eyelink('EyeAvailable');
% if win.el.eye_used==win.el.BINOCULAR,                           % (!TODO!) this I do not know yet
%                 win.el.eye_used = win.el.LEFT_EYE;
% end

Eyelink('WaitForModeReady', 500);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rvals               = [randsample(2:1:win.total_tact-1-win.total_test,...   % here we sample -win.total_test so we can add 1 in the next line (and so we can sample the space of trial stimulation every 1 instead of two)
                        win.total_test-1),win.total_tact-win.total_test];                         
rvals               = cumsum(1+diff([0 sort(rvals)]));                      % length of every trial, stimulation occurs with a flat hazard function, with win.total_test tests, and win.total_tact stimulations. sum(diff([0 sort(rvals)])) is then equal to  win.total_tact
win.tnseq           = diff([0 sort(rvals)]);    

nBlocks             = win.exp_trials./win.t_perblock;                       % # experimental block counting the first test one
win.block_cross     = repmat(randsample([1 2],2),1,nBlocks/2);              % 1 uncrossed 2 crossed
nTrials             = win.exp_trials;

win.block_start     = repmat([1,zeros(1,win.t_perblock-1)],1,nBlocks);      % Trials that are block start
trial_hand          = randsample(repmat([1 2],1,(win.total_tact-win.total_test*2)/2),...
                    win.total_tact-win.total_test*2);
win.trial_trigger   = nan(1,win.total_tact);

nT                  = 1;
win.seqstart        = [1 rvals(1:end-1)+1];
seqstart_hand       = [1 rvals(1:end-1)+1-[2:2:2*449]];
for bb = 1:nBlocks
    for tt = 1:win.t_perblock
        if win.block_cross(bb) == 1 % uncrossed
            win.trial_trigger(win.seqstart(nT):win.seqstart(nT)+win.tnseq(nT)-3) = ...
                trial_hand(seqstart_hand(nT):seqstart_hand(nT)+win.tnseq(nT)-3);      
            win.trial_trigger(win.seqstart(nT)+win.tnseq(nT)-2:win.seqstart(nT)+win.tnseq(nT)-1) = ...
                cell2mat(randsample({[9,10],[14,13]},1));
        elseif win.block_cross(bb) == 2 % crossed
            win.trial_trigger(win.seqstart(nT):win.seqstart(nT)+win.tnseq(nT)-3) = ...
                trial_hand(seqstart_hand(nT):seqstart_hand(nT)+win.tnseq(nT)-3)+4;
            win.trial_trigger(win.seqstart(nT)+win.tnseq(nT)-2:win.seqstart(nT)+win.tnseq(nT)-1) =....
                cell2mat(randsample({[17,18],[22,21]},1));
        end
       nT = nT+1;
    end
end


% for SOA between tactile stimulation we also used a flat hazard function
% so after the minimum length al latencies
soax                        = round(1000*(win.tact_minlat+...
                                (-1./win.decay .* log(rand([1 win.total_tact])))))/1000;
soax(soax>win.tact_max_lat) = win.tact_max_lat;
% last soa comes from the distribution of possible TOJ SOA
win.trial_tactsoa          = soax;
win.trial_tactsoa(rvals)   = randsample(repmat([20 30 50 90 170 330 650 1290 1800]/1000,1,nBlocks),nTrials); % TODO:this I have to think a little more, specially for the latencies abobe win.tact_minlat
% win.trial_tact_visSOA   = win.tact_visfix + win.tact_visrnd*rand(1,win.exp_trials);

%% remove variables that will not be used again
clear IsConnected wavedata1 wavedata2 exp_path nrchannels bb nT soax tt 
clear OpenError rvals seqstart_hand wrect freq dur                                  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE ACTUAL EXPERIMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fixdots = [win.fixpos(1)-win.fixrad,win.fixpos(2)-win.fixrad,...
        win.fixpos(1)+win.fixrad,win.fixpos(2)+win.fixrad;...
        win.tgtpos(1,1)-win.fixrad,win.tgtpos(1,2)-win.fixrad,...
        win.tgtpos(1,1)+win.fixrad,win.tgtpos(1,2)+win.fixrad;...
        win.tgtpos(2,1)-win.fixrad,win.tgtpos(2,2)-win.fixrad,...
        win.tgtpos(2,1)+win.fixrad,win.tgtpos(2,2)+win.fixrad];
resrect = [win.tgtpos(1,1)-win.fixrad*10,win.tgtpos(1,2)-win.fixrad*10,...
        win.tgtpos(1,1)+win.fixrad*10,win.tgtpos(1,2)+win.fixrad*10;...
        win.tgtpos(2,1)-win.fixrad*10,win.tgtpos(2,2)-win.fixrad*10,...
        win.tgtpos(2,1)+win.fixrad*10,win.tgtpos(2,2)+win.fixrad*10];
b                   = 0;                                                    % block flag
ntt                 = 1;
for nT = 1:nTrials                                                          % loop throught the experiment trials
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BLOCK START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if  win.block_start(nT) == 1                                            % if it is a trial that starts a block   
        PsychPortAudio('Stop', pahandle)                                    % if the white noise was on, we stop it here
        
        if nT ==1 && win.block_cross(nT)==1                                % practice trials and uncrossed
            draw_instructions_and_wait(texts.txt6,win.bkgcolor,win.hndl,win.in_dev,1)
        elseif nT ==1 && win.block_cross(nT)==2                       % practice trials amd crossed, this according to randomization should be unnecesary
            draw_instructions_and_wait(texts.txt7,win.bkgcolor,win.hndl,win.in_dev,1)
    
        elseif win.block_cross(nT)==1 
            texts.txt8    = double(['Block ' num2str(b) '/' num2str(nBlocks+1) ' beendet \n Pause \n  F' 252 'r den n' 228 ... 
            'chsten Block bitte die H' 228 'nde parallel positionieren (parallel). \n Zum Fortfahren die ' texts.txtdev]);
         draw_instructions_and_wait(texts.txt8,win.bkgcolor,win.hndl,win.in_dev,1)
        elseif  win.block_cross(nT)==2         % crossed
            texts.txt8    = double(['Block ' num2str(b) '/' num2str(nBlocks+1) ' beendet \n Pause \n  F' 252 'r den n' 228 ... 
            'chsten Block bitte die H' 228 'nde ' 252 'berkreuzen (crossed). \n Zum Fortfahren die ' texts.txtdev]);
         draw_instructions_and_wait(texts.txt8,win.bkgcolor,win.hndl,win.in_dev,1)
        end
           
        b = b+1;
        if nT>1 %&& ismember(nT, win.t_perblock+win.test_trials+1:win.calib_every*win.t_perblock:nTrials)                              % we calibrate every two small blocks
            EyelinkDoTrackerSetup(win.el);
        end
        
        if win.block_cross(nT)==1 
            DrawFormattedText(win.hndl, texts.txt10,'center','center',255,55);
        else
            DrawFormattedText(win.hndl,texts.txt11,'center','center',255,55);
        end
        Screen('Flip', win.hndl);
        if win.in_dev == 1                                                              
            waitForKB_linux({'space'});                                           
        elseif win.in_dev == 2
            GetClicks(win.hndl,0);                                                      
        end
        Screen('FillRect', win.hndl, win.bkgcolor);                         % remove what was writte or displayed
        Screen('Flip', win.hndl);
        Eyelink('WaitForModeReady', 50);
        PsychPortAudio('Start', pahandle, 0, 0, 1);
        EyelinkDoDriftCorrect2(win.el,win.res(1)/2,win.res(2)/2,1,1)          % drift correction 
       % beginning of each block we have 10 seconds of baseline
        Eyelink('Command','record_status_message ''Block %d Trial %d''',b,nT);
        Eyelink('WaitForModeReady', 25);
        Eyelink('message','TRIALID %d', nT);                                    % message about trial start in the eye-tracker
        Eyelink('WaitForModeReady', 25);
        Screen('FillOval',  win.hndl, 128, fixdots(1,:))
        Screen('Flip', win.hndl); 
        Eyelink('StartRecording');
        Eyelink('WaitForModeReady', 25);
        Eyelink('command', '!*write_ioport 0x378 4')
         WaitSecs(5)
    end
     Eyelink('command', '!*write_ioport 0x378 0')
    % TRIAL START
    if win.block_start(nT) == 0
%         PsychPortAudio('Start', pahandle, 0, 0, 1);               % starts the white noise, third input is set to 0 so it loops until is sopped
        Eyelink('Command','record_status_message ''Block %d Trial %d''',b,nT);
        Eyelink('WaitForModeReady', 25);
        Eyelink('message','TRIALID %d', nT);                                    % message about trial start in the eye-tracker
        Eyelink('WaitForModeReady', 25);
        Eyelink('StartRecording');
    end
      
    Screen('FillOval',  win.hndl, 256, fixdots(1,:))                             % Change in fixation dot indicates the start of the trial sequence
    Screen('Flip', win.hndl);   
    
    while 1                                                                %  check that the gaze is on the fixation point before starting with the sequence of stimulation
        [data,type] = get_ETdata;
            if type ==200 % start fixation
                if abs(data.gx(1)-win.fixpos(1,1))<win.center_thershold && ...
                    abs(data.gy(1)-win.fixpos(1,2))<win.center_thershold 
                   break
                end
            end
    end
    WaitSecs(.1);
     
    last_tact      = GetSecs;
   
    Eyelink('message','METATR block %d',win.block_cross(nT));               % block condition
    Eyelink('WaitForModeReady', 25);
    Eyelink('message','METATR block_start %d',win.block_start(nT));         % if it was the first image in the block
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TACTILE STIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for cs = 1:win.tnseq(nT)
        if cs==1
            Eyelink('message','SYNCTIME'); 
        end
        while GetSecs<last_tact+win.trial_tactsoa(ntt);
            continue
        end
        Eyelink('command', '!*write_ioport 0x378 %d',win.trial_trigger(ntt));                        
        last_tact = WaitSecs(win.stim_dur);
        Eyelink('command', '!*write_ioport 0x378 %d',0);                        % flush the parallel port
        last_trig = win.trial_trigger(ntt);
        ntt = ntt+1; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUBJECT GAZE RESPONSE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WaitSecs(.25)                                                            % 500 ms to present the response display
    
    Screen('FrameRect', win.hndl, [200 200 200], resrect',3)
    Screen('FillOval',  win.hndl, [255 0 0;200 200 200;200 200 200]', fixdots')
    Screen('Flip', win.hndl);
    Eyelink('command', '!*write_ioport 0x378 %d',96);                       % response display
    tvis = GetSecs;
    
    answ = NaN; tAns = [];
    while GetSecs-tvis<5                                                    %  check that the gaze moves to the left or right response sqare
        [data,type] = get_ETdata;
            if type ==6 % start fixation
                if abs(data.genx(1)-win.tgtpos(1,1))<win.center_thershold && ...
                    abs(data.geny(1)-win.tgtpos(1,2))<win.center_thershold 
                    answ = 1;
                    tAns = GetSecs;
                    Eyelink('message','TRIAL_RESULT %d',answ);  
                   break
                elseif abs(data.genx(1)-win.tgtpos(2,1))<win.center_thershold && ...
                    abs(data.geny(1)-win.tgtpos(2,2))<win.center_thershold 
                    answ = 2;
                    tAns = GetSecs;
                    Eyelink('message','TRIAL_RESULT %d',answ);  
                   break
                end
            end
    end
    Eyelink('command', '!*write_ioport 0x378 %d',0);
    Screen('FillOval',  win.hndl, win.fixcolor, fixdots(1,:))
    
%   TOJ Results;
    if (ismember(last_trig,[13,18]) && answ==1) || (ismember(last_trig,[10,21]) && answ==2)
        win.correct(nT) = 1;
    elseif (ismember(last_trig,[13,18]) && answ==2) || (ismember(last_trig,[10,21]) && answ==1)
        win.correct(nT) = 0;
    end
    win.response(nT)= answ;
    if isempty(tAns)
        win.RT(nT) = NaN;
        win.correct(nT) = NaN;
        Eyelink('message','TRIAL_RESULT 0');  
    else
        win.RT(nT)      = tAns-tvis;
    end
    
    Eyelink('WaitForModeReady', 50);
    Eyelink('StopRecording');
    Screen('Flip', win.hndl);
    Eyelink('WaitForModeReady', 50);
    save([pathEDF,win.fnameEDF(1:end-3),'mat'],'win')
end
win.end_time = clock;
PsychPortAudio('Stop', pahandle);                                           % Stop the white noise after the last trial

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finishing EDF file and transfering info (in case experiment is interrupted
% this can be run to save the eye-tracking data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Task Iteration done; save files, restore stuff, DON'T clear vars %%%%
Eyelink('CloseFile');
Eyelink('WaitForModeReady', 500); % make sure mode switching is ok
if ~win.DoDummyMode
    % get EDF->DispPC: file size [bytes] if OK; 0 if cancelled; <0 if error
    rcvStat = Eyelink('ReceiveFile', win.fnameEDF, pathEDF,1);
    if rcvStat > 0 % only sensible if real connect
        fprintf('EDF received to %s (%.1f MiB).\n',pathEDF,rcvStat/1024^2);
    else
        fprintf(2,'EDF file reception error: %d.\n', rcvStat);
    end
end

save([pathEDF,win.fnameEDF(1:end-3),'mat'],'win')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLOSING ALL DEVICES, PORTS, ETC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Most of the commands below are important for being able to restart the
% experiment, in the case the experiment crushs they should be entered
% manually

% A known issue: Eyelink('Shutdown') crashing Matlab in 64-bit Linux
% cf. http://tech.groups.yahoo.com/group/psychtoolbox/message/12732
% not anymore it seems
%if ~IsLinux(true), Eyelink('Shutdown'); end
PsychPortAudio('Close', pahandle);                                          % close the audio device
Eyelink('Shutdown');                                                        % close the link to the eye-tracker

Screen('CloseAll');                                                         % close the PTB screen
Screen('Preference','Verbosity', prevVerbos);                               % restore previous verbosity
Screen('Preference','VisualDebugLevel', prevVisDbg);                        % restore prev vis dbg
ListenChar(1)                                                               % restore MATLAB keyboard listening (on command window)
