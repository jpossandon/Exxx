
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exxx_experiment.m 
% Rapid serial tactile presentation + TOJ task
% This is an EEG/Eye-tracking/Tactile-stimulation experiment
% The experiment consists in the presentation of tactile stimulus to the left
% hand in uncrossed or crossed positions, meanwhile gaze is directed to the 
% center of the screen, at the middle between the position of the hands.
% Tactile stimulation occurs in a sequence of variable length according to a
% flat hazard function. During the sequence stimulus occured randomly to the
% left or the right hand at different SOA, with a minimmum of 750 ms (to
% prevent overlap of ERPs or TFRs). At the end of the sequence, the last
% two stimuli are always to different hands and at asynchronies that permit
% to calibration of TOJ performance curve. Subject respons which side was 
% stimulated by an eye movement to a left or right box displayed in the screen.
%
% Temporal profile
%
%     ts(l/r)--------ts(l/r)---...---ts(l/r)-----ts(!previous)-------response
% ms:       750+rnd(500)   750+rnd(500)     constant
%                                           stimulus
%                                           distribution:(10,30,60,180,360)
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
% time delaz to start audio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all                                                                 % we clear parameters?
% this is for debugging
win.DoDummyMode             = 0;                                            % (1) is for debugging without an eye-tracker, (0) is for running the experiment
win.stim_test               = 1;                                            % (1) for testing the stimulators (always when doing the experiment), (0) to skip
% PsychDebugWindowConfiguration(0,0.7);                                       % this is for debugging with a single screen

% Screen parameters
win.whichScreen             = 0;                                            % (CHANGE?) here we define the screen to use for the experiment, it depend on which computer we are using and how the screens are conected so it might need to be changed if the experiment starts in the wrong screen
win.FontSZ                  = 20;                                           % font size
win.bkgcolor                = 0;                                          % screen background color, 127 gray
win.Vdst                    = 66;                                           % (!CHANGE!) viewer's distance from screen [cm]         
win.res                     = [1920 1080];                                  %  horizontal x vertical resolution [pixels]
win.wdth                    = 51;                                           %  51X28.7 cms is teh size of Samsung Syncmaster P2370 in BPN lab EEG rechts
win.hght                    = 28.7;                                         % 
win.pixxdeg                 = win.res(1)/(2*180/pi*atan(win.wdth/2/win.Vdst));% 
win.center_thershold        = 3*win.pixxdeg;                                % distance from the midline threshold for gaze contingent end of trial
win.fixpos                  = [win.res(1)/2, win.res(2)/2];
win.tgtpos                  = [win.fixpos(1,1)/2+win.fixpos(2,1)/2 win.fixpos(3,1)/2+win.fixpos(1,1)/2];
win.fixcolor                = 127;
win.fixrad                  = 10;
% Tactile stimulation settings (using the box stimulator)
win.tact_freq               = 200;                                          % frequency of stimulation in Hz
win.stim_dur                = .025;                                         % duration of tactile stimulation. The vibrator takes some time to stop its motion so for around 50 ms we use 25 ms of stimulation time (ask Tobias for the exact latencies they have measured)
win.tact_minlat             = .8;
win.tact_rndlat             = .4;

win.halflife                = 30;
win.decay                   = log(2)./halflife;         %lambda
win.min_length              = 1;


% Blocks and trials
win.exp_trials              = 2400;
win.t_perblock              = 30;
win.calib_every             = 4; 

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
    exp_path                = '/Users/jossando/trabajo/Exxx/';              % path in my mac
else
    exp_path                = '/home/th/Experiments/Exxx/';
end

win.s_n                     = input('Subject number: ','s');                % subject id number, this number is used to open the randomization file
win.fnameEDF                = sprintf('s%02d.EDF',str2num(win.s_n));       % EDF name can be only 8 letters long, so we can have numbers only between 01 and 99
pathEDF                     = [exp_path 'data/' sprintf('s%02d/',str2num(win.s_n))];                           % where the EDF files are going to be saved
if exist([pathEDF win.fnameEDF],'file')                                         % checks whether there is a file with the same name
    rp = input(sprintf('!Filename %s already exist, do you want to overwrite it (y/n)?',win.fnameEDF),'s');
    if (strcmpi(rp,'n') || strcmpi(rp,'no'))
        error('filename already exist')
    end
end
if ~isdir(pathEDF)
    mkdir(pathEDF)
end
win.s_age                   = input('Subject age: ','s');
win.s_hand                  = input('Subject handedness for writing (l/r): ','s');
win.s_gender                = input('Subject gender (m/f): ','s');
setStr                      = sprintf('Subject %d\nAge %s\nHandedness %s\nGender %s\n',win.s_n,win.s_age,win.s_hand,win.s_gender); % setting summary
fprintf(setStr); 

AssertOpenGL();                                                             % check if Psychtoolbox is working (with OpenGL) TODO: is this needed?
ClockRandSeed();                                                            % this changes the random seed

[IsConnected, IsDummy] = EyelinkInit(win.DoDummyMode);                      % open the link with the eyetracker
assert(IsConnected==1, 'Failed to initialize EyeLink!')
 
% ListenChar(2)                                                             % disable key listening by MATLAB windows(CTRL+C overridable)

prevVerbos = Screen('Preference','Verbosity', 2);                           % this two lines it to set how much we want the PTB to output in the command and display window 
prevVisDbg = Screen('Preference','VisualDebugLevel',3);                     % verbosity-1 (default 3); vdbg-2 (default 4)
Screen('Preference', 'SkipSyncTests', 2)                                    % for maximum accuracy and reliability

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START PTB SCREEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[win.hndl, win.rect]        = Screen('OpenWindow',win.whichScreen,win.bkgcolor);   % starts PTB screen

[win.cntr(1), win.cntr(2)] = WindowCenter(win.hndl);                        % get where is the display screen center
Screen('BlendFunction',win.hndl, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     % enable alpha blending for smooth drawing
HideCursor(win.hndl);                                                       % this to hide the mouse
Screen('TextSize', win.hndl, win.FontSZ);                                   % sets teh font size of the text to be diplayed
KbName('UnifyKeyNames');                                                    % recommended, called again in EyelinkInitDefaults
win.start_time = clock;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START AUDIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq                        = 48000;                                        % noise frequency in Hz
 dur                         = 5;                                           % noise duration in seconds
 wavedata1                   = rand(1,freq*dur)*2-1;                         % the noise
 wavedata2                   = sin(2.*pi.*win.tact_freq.*[0:1/freq:dur-1/freq]);      % the stimulus
 wavedata                    = [wavedata1 ; wavedata2];     
nrchannels                  = 2;                                            % ditto

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
s = PsychPortAudio('GetStatus', pahandle);                                  % Status of the port, necessary later to defice wheter to start or stop audio


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUCTIONS IN GERMAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if win.in_dev == 1
    txtdev = ['Leertaste dr' 252 'cken (press Space)'];
elseif win.in_dev == 2
    txtdev = 'Maustate klicken';
end

% This is the safest way so umlaut are corretly displayed UTF. 
% code for umlauts 228(a),252(u)
 txt1     = double(['Die Funktionsf' 228 'higkeit der Ger' 228 'te muss '...    
             252 'berpr' 252 'ft werden\n' txtdev ' zum Fortfahren..']);             
txt2     = double(['Stimulatoren werden auf den Handr' 252 'cken besfestigt ...\n' txtdev]);
txt3     = double(['Der linke Stimulator wird getestet (vibriert drei mal), ...\n danach ' txtdev]);
txt4     = double(['Der rechte Stimulator wird getestet (vibriert drei mal), ...\n danach ' txtdev]);
txt6     = double(['Beginn des Experiments \n Block 1 \n Die H' 228 ...
            'nde bitte Links positionieren. \n Zum Beginnen die ' txtdev]);
txt7     = double(['Beginn des Experiments\n Test Block \n  F' 252 'r den n' 228 ...
            'chsten Block bitte Rechts positionieren. \n Zum Fortfahren die ' txtdev]);      
        
txt10    = double(['F' 252 'r den n' 228 'chsten Block die H' 228 ...
        'nde Links positionieren. Zum Fortfahren die ' txtdev]);
txt11    = double(['F' 252 'r den n' 228 'chsten Block die H' 228 ...
        'nde Rechts positionieren. Zum Fortfahren die '  txtdev]);
txt12    = double(['Links oder Rechts']);
%these are for debugging
handstr  = {'Left','Right','','','Left','Right'};
crossstr = {'Uncrossed','Crossed'};
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EYE-TRACKER SETUP, OPEN THE EDF FILE AND ADDS INFO TO ITS HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[win.el, win.elcmds] = setup_eyetrackerxxx(win, 1);                            % Setups the eye-tracker with the before mentioned parameters

DrawFormattedText(win.hndl,txt1,'center','center',255,55);                  % This 'draws' the text, nothing is displayed until Screen flip
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
Eyelink('Command', sprintf('add_file_preamble_text ''%s''', setStr));       % this adds the information about subject to the end of the header  
wrect = Screen('Rect', win.hndl);                                           
Eyelink('Message','DISPLAY_COORDS %d %d %d %d', 0, 0, wrect(1), wrect(2));  % write display resolution to EDF file

ListenChar(2)                                                               % disable MATLAB windows' keyboard listen (no unwanted edits)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIMULATOR TEST
% ALWAYS: left-hand is # 1 and right-hand is # 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if win.stim_test
    DrawFormattedText(win.hndl,txt2,'center','center',255,55);
    Screen('Flip', win.hndl);
    if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
        waitForKB_linux({'space'});                                             % press the space key in the keyboard
    elseif win.in_dev == 2      
        [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
    end

    DrawFormattedText(win.hndl,txt3,'center','center',255,55);                  % TEST LEFT STIMULATOR (three times)
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

    DrawFormattedText(win.hndl,txt4,'center','center',255,55);                  % TEST RIGHT STIMULATOR (three times)
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

[image,map,alpha]   = imread([exp_path 'stimuli/blackongrt.jpg']);          % drift correction dot image
fixIndex            = Screen('MakeTexture', win.hndl, image);               % this is one of the way PTB deals with images, the image matrix is transformed in a texture with a handle that can be user later to draw the image in theb PRB screen


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nBlocks             = win.exp_trials./win.t_perblock;                       % # experimental block without counting the first test one
% win.blockcond       = reshape(repmat(randsample(repmat([1 2,5 6,9 10],1,nBlocks/6),nBlocks),...
%                          win.t_perblock,1),1,win.exp_trials);          % 1 (gc-l), 2 (gc-r),  5 (gr-ll), 6 (gr-l),  9 (gl-r) ,10 (gl-rr)
win.blockcond       = [reshape(randsample(repmat([1;5;9], win.exp_trials/2/3,1),win.exp_trials/2,1),win.t_perblock,nBlocks/2),...
    reshape(randsample(repmat([2;6;10], win.exp_trials/2/3,1),win.exp_trials/2,1),win.t_perblock,nBlocks/2)];
win.blockcond       = win.blockcond(:,randsample(nBlocks,nBlocks));
win.blockcond       = win.blockcond(:)';
nTrials             = win.exp_trials;                                       % Total # of trial
nBlocks             = win.exp_trials./win.t_perblock;                       % # experimental block without counting the first test one
win.block_start     = repmat([1,zeros(1,win.t_perblock-1)],1,nBlocks);     % Trials that are block start
win.visTagtdevx     = round(randn(1,win.exp_trials)*win.stdvistgt*win.pixxdeg);
win.visTagtdevy     = round(randn(1,win.exp_trials)*win.stdvistgt*win.pixxdeg);
win.visTagtSide     =sign(win.visTagtdevx)/2+1.5;                           % 1 - Left, 2 - Right

outy = find(abs(win.visTagtdevy)>win.rect(4)/2);
if ~isempty(outy)
    win.visTagtdevy(outy) = round(ones(1,length(outy))*win.stdvistgt*win.pixxdeg);
end

outx = find(abs(win.visTagtdevx)>win.rect(3)/2);
if ~isempty(outx)
    win.visTagtdevy(outx) = round(ones(1,length(outx))*win.stdvistgt*win.pixxdeg);
end
    win.trial_tactsoa   = win.tact_fixlat + win.tact_rndlat*rand(1,win.exp_trials);
win.trial_tact_visSOA   = win.tact_visfix + win.tact_visrnd*rand(1,win.exp_trials);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE ACTUAL EXPERIMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b                   = 0;                                                    % block flag
for nT = 1:nTrials                                                          % loop throught the experiment trials
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BLOCK START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if  win.block_start(nT) == 1                                                % if it is a trial that starts a block   
        if s.Active                                                         % if the white noise was on, we sopt it here
            PsychPortAudio('Stop', pahandle)
        end
        if nT ==1 && ismember(win.blockcond(nT),[1,5,9])                                % practice trials and uncrossed
            draw_instructions_and_wait(txt6,win.bkgcolor,win.hndl,win.in_dev,1)
        elseif nT ==1 && ismember(win.blockcond(nT),[2,6,10])                     % practice trials amd crossed, this according to randomization should be unnecesary
            draw_instructions_and_wait(txt7,win.bkgcolor,win.hndl,win.in_dev,1)
    
        elseif ismember(win.blockcond(nT),[1,5,9])
            txt2    = double(['Block ' num2str(b) '/' num2str(nBlocks+1) ' beendet \n Pause \n  F' 252 'r den n' 228 ... 
            'chsten Block bitte die H' 228 'nde Links positionieren. \n Zum Fortfahren die ' txtdev]);
         draw_instructions_and_wait(txt2,win.bkgcolor,win.hndl,win.in_dev,1)
        else         % crossed
            txt2    = double(['Block ' num2str(b) '/' num2str(nBlocks+1) ' beendet \n Pause \n  F' 252 'r den n' 228 ... 
            'chsten Block bitte die H' 228 'nde Rechts positionieren. \n Zum Fortfahren die ' txtdev]);
         draw_instructions_and_wait(txt2,win.bkgcolor,win.hndl,win.in_dev,1)
        end
           
       
        b = b+1;
        if nT>1 %&& ismember(nT, win.t_perblock+win.test_trials+1:win.calib_every*win.t_perblock:nTrials)                              % we calibrate every two small blocks
            EyelinkDoTrackerSetup(win.el);
        end
        
        if ismember(win.blockcond(nT),[1,5,9])
            DrawFormattedText(win.hndl, txt10,'center','center',255,55);
        else
            DrawFormattedText(win.hndl,txt11,'center','center',255,55);
        end
        Screen('Flip', win.hndl);
        if win.in_dev == 1                                                              
            waitForKB_linux({'space'});                                           
        elseif win.in_dev == 2
            GetClicks(win.hndl,0);                                                      
        end
        Screen('FillRect', win.hndl, win.bkgcolor);                         % remove what was writte or displayed
        Screen('DrawTexture', win.hndl, fixIndex);                          % drift correction image
        Screen('Flip', win.hndl);
        Eyelink('WaitForModeReady', 500);
        EyelinkDoDriftCorrect2(win.el,win.res(1)/2,win.res(2)/2,0)          % drift correction 
       
    end
    Eyelink('Command','record_status_message ''Block %d Image %d''',b,nT);
    Eyelink('WaitForModeReady', 500);
       
    Eyelink('StartRecording');
    if ismember(win.blockcond(nT),[1,2])   % center gaze
        caux = 1;
    elseif ismember(win.blockcond(nT),[9,10])   % left gaye
        caux = 2;
    else
        caux = 3;
    end
    
    if ismember(win.blockcond(nT),[1,5,9])   %left
        hpos =1;
    else
        hpos=2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TRIAL Start
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     fixrect = [win.fixpos(caux,1)-win.fixrad,...
        win.fixpos(caux,2)-win.fixrad,...
        win.fixpos(caux,1)+win.fixrad,...
        win.fixpos(caux,2)+win.fixrad];
    Screen('FillOval',  win.hndl, 64, fixrect)
    Screen('Flip', win.hndl);   
    
    PsychPortAudio('Start', pahandle, 0, 0, 1);               % starts the white noise, third input is set to 0 so it loops until is sopped
                                                                             %  it goes on onlz when the sound has started aproximatedlz 50 ms in the linux computer EEG nue
   
     while 1                                                             % images change contingent to the end of a fixation and the horixzontal position, we already did the waiting at the end of previous trial
            [data,type] = get_ETdata;
                if type ==200 % start fixation
%                     data
                     if abs(data.gx(1)-win.fixpos(caux,1))<win.center_thershold && ...
                         abs(data.gy(1)-win.fixpos(caux,2))<win.center_thershold 
%                     if abs(data.genx(win.el.eye_used)-win.res(1)./2)<win.center_thershold % (!TODO!) check this
%                        WaitSecs(.04+randsample(.01:.01:.1,1));              % this is a lag+jitter so the change of the image occurs after saccadic supression betwenn .05 and .150 sec
%                         display('Gaze On')
                        break
                     end
                end
     end
%      WaitSecs(1)
     
    Screen('FillOval',  win.hndl, win.fixcolor, fixrect)
    Eyelink('message','TRIALID %d', nT);                                % message about trial start in the eye-tracker
    Screen('Flip', win.hndl);   
    tstart      = GetSecs;
   
  
   
    Eyelink('message','METATR block %d',win.blockcond(nT));                 % block condition
    Eyelink('WaitForModeReady', 50);
    Eyelink('message','METATR block_start %d',win.block_start(nT));         % if it was the first image in the block
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TACTILE STIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    while GetSecs<tstart+win.trial_tactsoa(nT)
        continue
    end
    Eyelink('command', '!*write_ioport 0x378 %d',win.blockcond(nT));                        
    WaitSecs(win.stim_dur)
    Eyelink('command', '!*write_ioport 0x378 %d',0);                        % flush the parallel port
    while GetSecs<tstart+win.trial_tactsoa(nT)+win.trial_tact_visSOA
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUAL STIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fixtgt = [win.fixpos(caux,1)-win.fixrad, win.fixpos(caux,2)-win.fixrad,...
        win.fixpos(caux,1)+win.fixrad,win.fixpos(caux,2)+win.fixrad;
        win.tgtpos(hpos)+win.visTagtdevx(nT)-win.fixrad, win.fixpos(caux,2)+win.visTagtdevy(nT)-win.fixrad,...
        win.tgtpos(hpos)+win.visTagtdevx(nT)+win.fixrad,win.fixpos(caux,2)+win.visTagtdevy(nT)+win.fixrad];
    Screen('FillOval',  win.hndl, win.fixcolor, fixtgt')
    Screen('Flip', win.hndl);
    Eyelink('command', '!*write_ioport 0x378 %d',96);  
    tvis = GetSecs;
    
    Screen('FillOval',  win.hndl, win.fixcolor, fixrect)
    WaitSecs(.05);
     Eyelink('command', '!*write_ioport 0x378 %d',0);         
    Screen('Flip', win.hndl);
     WaitSecs(.05);
     PsychPortAudio('Stop', pahandle);
    click = 0;
    while GetSecs<tvis+3 
        [x,y,buttons] = GetMouse(win.hndl);    
        if sum(buttons)
            Eyelink('command', '!*write_ioport 0x378 %d',92);  
            win.RT(nT) = GetSecs-tvis;
            if length(find(buttons))==1
                win.respose(nT) = find(buttons);
            else
                win.respose(nT)= NaN;
            end
            click =1;
            while sum(buttons)
                [x,y,buttons] = GetMouse(win.hndl);   
            end
            break
        end
    end
    if ~click
        win.RT(nT) = NaN;
        win.respose(nT)= NaN;
    end
    
    Eyelink('WaitForModeReady', 50);
    Eyelink('StopRecording');
    Eyelink('WaitForModeReady', 50);
      Eyelink('command', '!*write_ioport 0x378 %d',0);                  
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
% fclose(obj);                                                                % close the serial port
ListenChar(1)                                                               % restore MATLAB keyboard listening (on command window)
