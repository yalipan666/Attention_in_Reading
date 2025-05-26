%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading Experiment with Rapid Tagging frequency:
% whole dynamics of attention shift in reading(wrd freq)
% Both target and postarget words will be tagged with two different frequency bands, which are switched between participants
% Press Q -quit: to exit exp
% Press S -continue: to skip eyechecker of the start box
% Press C -continue: to skip eyechecker of the end box and continue to
% the next trial
% Press E -eye: to start eyelink setup, calibration or/and validation
% 20230208 Yai Pan
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Read_FullAttention_2wrd
%%======= double check parameters before MEG =======%%
cfg.debugmode = 0; % runnning on local computer for rapidmode testing without triggers
if cfg.debugmode
    Screen('Preference', 'SkipSyncTests', 1); %must be 0 during experiment
    cfg.el.eyelink = 0;              %eyelink on/off
    cfg.el.override = 1;          %No eyelink in actual experiment, use only in case of fault
    cfg.DataPixxOnly = 0;      %for rapidmode testing without triggers with propixx
    tt = 0.1;                 % just to speed up the presentation during debug mode
else
    Screen('Preference', 'SkipSyncTests', 1); %must be 0 during experiment
    cfg.el.eyelink = 1;              %eyelink on/off
    cfg.el.override = 0 ;         %No eyelink in actual experiment, use only in case of fault
    cfg.DataPixxOnly = 1;      %for rapidmode testing without triggers with propixx
    tt = 1;            % just to speed up the presentation during debug mode
end

% Input dialog
prompt = {'Participant ID:(1,2,3...)', ...
    'Test date:', ...
    'Participant code (b123):', ...
    'Screen width (cm): ', ...
    'Screen height (cm): ', ...
    'Screen distance (cm): '};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'46','0901','b5e8','77.5','43.5','118'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
cfg.SubID = str2double(answer{1});
cfg.Date = answer{2};
cfg.SubCode = answer{3};
%Physical screen parameters (cm) (input)
cfg.width = str2double(answer{4});   %projection screen width in cm
cfg.height = str2double(answer{5});  %projection screen height in cm
cfg.dist = str2double(answer{6});    %distance from subject eye to screen in cm
% version for sentence set, tagging frequency [f1 f2], and button response
tmp = mod(cfg.SubID,4);
if tmp == 1
    cfg.SentVersion = 1;
    cfg.FreqVersion = 1; %[60 65]
elseif tmp == 2
    cfg.SentVersion = 1;
    cfg.FreqVersion = 2; %[65 60]
elseif tmp == 3
    cfg.SentVersion = 2;
    cfg.FreqVersion = 1;
else 
    cfg.SentVersion = 2;
    cfg.FreqVersion = 2;
end 
% version for button press is the same as sentence set

%% ========Frequency Tagging Parameters =============%%%%%%%%%%%
% odd sub: sentence_1, tag[60 65]Hz;even sub: sentence_2, tag[65 60]Hz
f1 = 60;
f2 = 65;
if cfg.FreqVersion == 1
    tagfreq = [f1 f2]; %[target postarget]
else
    tagfreq = [f2 f1]; %[target postarget]
end
cfg.TagFreq = tagfreq;


%% initializing logfile, need to change the filename
exp_dir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(exp_dir, 'Data');
% create results dir if it does not exist
if (~exist(resultsDir,'dir'))
    mkdir(resultsDir);
end
expdate_all = clock;
expdate = [num2str(expdate_all(1)) cfg.Date];
datafilename = [resultsDir filesep expdate '_' cfg.SubCode '.mat'];
if exist(datafilename,'file') && ~cfg.debugmode
    error('The data file exists! Please enter a different subject code.');
end

%% %%%%%%%%%%%%%%%=================PTB screen Initialization======================%%%%%%%%%%%%%%%%%%%%%
close all;
Screen('CloseAll');
AssertOpenGL;
PsychDefaultSetup(2);    % call some default settings for setting up Psychtoolbox
%Open screen
screens = Screen('Screens'); % Get the screen numbers
cfg.screenNumber = max(screens); %select screen
%%%  Get the size of the on screen window and set resolution
sc_info = Screen('Resolution', cfg.screenNumber);
resx = sc_info.width;
resy = sc_info.height;
cfg.resx = resx;
cfg.resy = resy;
% initializing keyboard setting
KbName('UnifyKeyNames'); % for easy use of Keyboard keys
cfg.el.eyelinkKey = KbName('E');  %Key used to toggle eyelink validation/calculation on or off during experiment.
cfg.el.SkipStartboxKey = KbName('S'); %skip eye check part in start box 跳过start box的校准，继续实验
cfg.el.SkipEndboxKey = KbName('C'); %skip eye check part in end box 跳过end box的校准，继续实验
escKey = KbName('Q');

%%%%%%%%%%%%%%%%%%%%%%%%================EyeLink settings====================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eyeused = {'LEFT_EYE';'RIGHT_EYE'};
eyeused_id = 1;
cfg.el.Eyeused = eyeused{eyeused_id};
cfg.el.edffile = [cfg.Date cfg.SubCode '.edf']; %EDF filename
%%%%%%=====Eyelink Initializationqq
if cfg.el.eyelink
    %add eyelink script folder (should be in main experiment folder)
    addpath([exp_dir filesep 'Eyelink']);
    %make directory if it doesn't already exist (local computer)
    cfg.el.eyedir = [exp_dir filesep 'Eyelink' filesep ];
    if ~exist(cfg.el.eyedir, 'dir')
        mkdir(cfg.el.eyedir);
    end
    %check whether files already exist for this subject/session
    if exist([exp_dir filesep 'Eyelink' filesep 'Data' filesep  cfg.el.edffile],'file')>0
        cont = input('Warning! Eyelink file will be overwritten, do you want to continue? (y/n) ','s');
        if cont == 'n'
            error('Session aborted')
        end
    end
    cfg.el_rect = [0 0 resx resy]; %% needed for the el_Set_Params function
    % Set parameters, start and calibrate eyelink
else %=% when eyelink is off
    if ~cfg.debugmode %is this is real experiment time, eyelink should be on
        if cfg.el.override
            warning('Eyelink not in use, continuing anyway...')
        else
            error('Eyelink not selected!')
        end
    end
end

%%%%%%%%============= open PTB window
cfg.ScrBgc = [0.5 0.5 0.5]; %middle gray screen
w = PsychImaging('OpenWindow', cfg.screenNumber,cfg.ScrBgc);
cfg.window = w;
%%% for the rapid drawing in the offscreen window
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');

% Query the maximum priority level
topPriorityLevel = MaxPriority(w);
if ~cfg.debugmode
    HideCursor;
    Priority(topPriorityLevel);
end

% enable alpha blending, also required for using offscreen window
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% Flip to clear
Screen('Flip', w);

% Query the frame duration
ifi = Screen('GetFlipInterval', w);
Para.ifi = ifi;

%%%%%%%%%%%%%%%%%%%%%%%%=================Stimuli settings=================%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.WordSpace = 0.43; % unit in visual angle, equal space between each word;
WordSpace = Deg2Pix_RIFT(cfg.WordSpace,cfg);
cfg.WordStart = 3*cfg.WordSpace; %% unit in visual angle, the start point of sentence
WordStart = Deg2Pix_RIFT(cfg.WordStart,cfg);
cfg.TextStyle = 1;                  %0=normal,1=bold,2=italic,4=underline,8=outline,32=condense,64=extend.
cfg.TextFont = 'Courier New';
cfg.TextSize = 20;  %0.43 va degrees
cfg.TextColor = BlackIndex(w);
white = WhiteIndex(w);

%%%Timing
cfg.fix_t = tt*1.2;                 %Duration (s) of the fixation
cfg.fix_jitter = tt*0.4;           %Baseline jitter. fixation duration is fix_t+jitter*[0..1]
cfg.iti = tt*0.5;                    %intertrial interval
cfg.rest = tt*30; %% rest time between blocks
cfg.fedbk = tt*1;

%%%%%%====== fixation cross
% Set size of the arms and linewidth
cfg.fixCross_r = 0.6; %radius of the fixation cross, unit in visual degrees
cfg.fixCross_r = Deg2Pix_RIFT(cfg.fixCross_r, cfg);
lineWidthPix = 3;
% get coordinates
xCoords = [-cfg.fixCross_r cfg.fixCross_r 0 0];
yCoords = [0 0 -cfg.fixCross_r cfg.fixCross_r];
CrossCoords = [xCoords; yCoords];

%%%%%%====== load in Stimulus matrix
load(['SentMat_' num2str(cfg.SentVersion) '.mat']); % SentMat
load(['TarWodLocFreq_' num2str(cfg.SentVersion) '.mat']); % TarWodLocFreq
load Question.mat %same questions for both versions of sentence sets

% parameteres for the reading paradigm
Para.SentMat = SentMat;
FlkWords = TarWodLocFreq(:,1);
if length(tagfreq) == 2
    FlkWords = TarWodLocFreq(:,[1 3]); %% location of the tagging word1 & word2
end
Para.FlkWords = FlkWords;
Para.TarWodLocFreq = TarWodLocFreq;

% condition matrix
Para.CondMat = TarWodLocFreq(:,2); %word freq condition for each sentence [(1--low, 2--high)]
%total number of trials
nTrials = size(TarWodLocFreq,1);
Para.nTrials = nTrials;
n_block = 4;
Para.BreakTrials = ceil(nTrials/n_block); %How many trials between breaks?

%Calculate the stimulus centers
xPos = resx/4;
yPos = resy/4;
pos_1 = round([xPos yPos]); % left upper
pos_2 = round([3*(xPos) yPos]); % right-upper
pos_3 = round([xPos 3*(yPos)]); % left-lower
pos_4 = round([3*(xPos) 3*(yPos)]); % right-lower
qcenters = [pos_1 ; pos_2 ; pos_3 ; pos_4]; % centers for each quandrant

%for every quandrant, calculate rects (useful for texture placement)
q_rects(1,:) = [0 0 qcenters(1,:)*2]; % left upper
q_rects(2,:) = [resx/2 0 resx resy/2]; % right-upper
q_rects(3,:) = [0 resy/2 resx/2 resy];  % left-lower
q_rects(4,:) = [resx/2 resy/2 resx resy]; % right-lower

%%%% define the eyelink monitor windows
%Parameters for fixation control, note conversion to pixels does not
%take into account rapidmode, because eyelink knows only projected screen
cfg.gaze_start = 0.2;  %% duration of gaze within start box before sentence onset
cfg.gaze_end = 0.2; %% duration of gaze within end box
cfg.box_r_el = 1; % actual visual degree for eye-link to monitor %% RADIUS of start and end box of sentences
dot_center_ptb = Deg2Pix_RIFT((cfg.box_r_el+cfg.WordSpace), cfg);
cfg.dot_r_ptb = 0.3*cfg.box_r_el; % small visual degree for ptb to present
dot_r = Deg2Pix_RIFT(cfg.dot_r_ptb, cfg); %% radius, unit in pixel, in small rapid mode
%%%  dot coords for ptb to present
coords_dot = [0 0 2*dot_r 2*dot_r];
DotCoords_start = zeros(4,4);
DotCoords_end = zeros(4,4);
for q = 1:4
    DotCoords_start(q,:) = floor(CenterRectOnPointd(coords_dot, q_rects(q,1)+dot_center_ptb, qcenters(q,2)));
    DotCoords_end(q,:) = floor(CenterRectOnPointd(coords_dot,qcenters(q,1), q_rects(q,4)-2*dot_center_ptb)); % bottom vertical end box
end
%%% box coords for eyelink monitor window, in big normal window not small rapidmode window
box_r_el = Deg2Pix(cfg.box_r_el,cfg); %% convert degree to pix
box_center_ptb = Deg2Pix(cfg.box_r_el+cfg.WordSpace,cfg);
cfg.el.startWindow = [0 0 2*box_r_el 2*box_r_el];
cfg.el.startWindow = floor(CenterRectOnPointd(cfg.el.startWindow,box_center_ptb,resy/2)); %% starting box small window
cfg.el.endWindow = floor(CenterRectOnPointd(cfg.el.startWindow,resx/2,resy-2*box_center_ptb)); %% starting box small window
endbox_color = white./4;

%%%%%%%%%%%%%%================Initalise Labjack / buttonbox settings===============%%%%%%%%%%%%%%%%%%%%%
if ~cfg.debugmode
    %%%% NAta Setting up
    if cfg.SentVersion == 1
        cfg.keyTrue = KbName('7&');  
        cfg.keyFalse = KbName('8*'); 
        fingers = {'True: index finger';'False: middle finger'};
    else
        cfg.keyTrue = KbName('8*');  
        cfg.keyFalse = KbName('7&');
        fingers = {'True: middle finger';'False: index finger'};
    end
    active = [cfg.keyTrue cfg.keyFalse]; %These are the left and right index fingers of the (5-button) NATA boxes
    keylist = zeros(1,256); %Set all keys to zero (ignore)
    keylist(active) = 1; %set active keys to 1 (listen)
    KbQueueCreate(0,keylist);%%Create queue, this is a time consuming operation (relatively), do while non-time critical
    KbQueueStart(); %Start listening
else
    %%% KEY response in debugmode
    cfg.keyTrue = KbName('J'); %% true
    cfg.keyFalse = KbName('K'); %% false
    fingers = {'True: J';'False: K'};
end

%%%%%%%%%%%%%%%%%==================Parallel Port IO & triggers settings===============%%%%%%%%%%%%%%%%%%%%%
% set up triggers
cfg.TriggerFix = 1;   % fixation onset
cfg.TriggerStartBox = 2; % start box onset
cfg.TriggerSentOn = 4;   %  sentence onset
cfg.TriggerSentOff = 8;  % sentence offset
cfg.TriggerITI = 16;   % ITI onset

% set up parallel port
if ~cfg.debugmode
    % initialize the parallel port
    cfg.PortAddress = hex2dec('BFF8');% MEG New PC; hex2dec('BFF8');% MEG Old PC
    cfg.ioObjTrig = io64;
    status = io64(cfg.ioObjTrig); % activation
    io64(cfg.ioObjTrig,cfg.PortAddress,0); %trigger 0 (reset)
end

%%%%%%%%%%%%%%%%%=================Propixx initialization======================%%%%%%%%%%%%%%%%%%%%%
cfg.ProPixxModel = 5;%different models of Propixx, 2 for tagging at 480Hz, 5 for tagging at 1440 Hz, 0 for normal model without tagging
cfg.ProPixxRefresh = 1440;
% Setup Propixx
if  ~cfg.debugmode || cfg.DataPixxOnly
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', cfg.ProPixxModel);
    Datapixx('RegWrRd');
end

%%%%%%%%%%%%%%%%%%=================photoDiode settings=================%%%%%%%%%%%%
cfg.diodeSize = 1.5; %radius of the photodiode in visual angle
diode_size_pix = Deg2Pix_RIFT(cfg.diodeSize, cfg);
%Get the placement of the photoDiode
%calculate diode positions -- left bottom one on each quadrant
diode_pos(:,1) = [0 resy/2-diode_size_pix diode_size_pix resy/2]';
diode_pos(:,2) = [resx/2 resy/2-diode_size_pix resx/2+diode_size_pix resy/2]';
diode_pos(:,3) = [0 resy-diode_size_pix diode_size_pix resy]';
diode_pos(:,4) = [resx/2 resy-diode_size_pix resx/2+diode_size_pix resy]';
if length(tagfreq) == 2 %need to plot the second pd for the 2 tagfreq
    %calculate diode positions -- right bottom one on each quadrant
    diode2_pos(:,1) = [resx/2-diode_size_pix resy/2-diode_size_pix resx/2 resy/2]';
    diode2_pos(:,2) = [resx-diode_size_pix resy/2-diode_size_pix resx resy/2]';
    diode2_pos(:,3) = [resx/2-diode_size_pix resy-diode_size_pix resx/2 resy]';
    diode2_pos(:,4) = [resx-diode_size_pix resy-diode_size_pix resx resy]';
    diode_pos = [diode_pos diode2_pos]; %4row * n columns
end

% create photoDiode
imageSize = [379 379];
ci = [199, 199, 199];     % center and radius of circle ([c_row, c_col, r])
[xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
mask_alpha = uint8((xx.^2 + yy.^2)<ci(3)^2);
mask_alpha(mask_alpha==1) = 255;  % Alpha layer: center is 255, surrounding is 0
diode = uint8(ones(size(mask_alpha))*255);
diode(:,:,2) = mask_alpha;
diode_tex = Screen('MakeTexture', w, diode);

% Get the maximal amount of frames to calculate the timecourse for
max_trialframes = 1000*round((cfg.fix_t+cfg.fix_jitter)/ifi); %max frames per trial---1000s
frame_mult = round(cfg.ProPixxRefresh/(1/ifi)); %every refresh is 12 frames for tagging at 1440hz with a 120 Hz monitor
%Effective presentation frequency. Should be 1440 for propixx
disp(['Effective refresh rate: ' num2str(cfg.ProPixxRefresh,6) 'Hz']);
if sum(cfg.TagFreq>cfg.ProPixxRefresh/2)>0
    warning('Presentation rate too low for the chosen flicker frequencies!')
end

%Frequency timecourse parameters
cfg.patch_amplitude = 0.5;
cfg.f_offset = 0;
cfg.patch_startPhase = [0 0]; % the start phase for both tagging signals
%initialize the table
freqTable = NaN(length(cfg.TagFreq),(max_trialframes*frame_mult));
for f = 1:length(cfg.TagFreq)
    patch_frequency = cfg.TagFreq(f);
    patch_angFreq = 2*pi*patch_frequency;
    start_time = (cfg.fix_t+cfg.fix_jitter)*-1;
    frametime = start_time:ifi/frame_mult:(max_trialframes*frame_mult)*(ifi/frame_mult)+start_time;
    frametime = frametime(1:max_trialframes*frame_mult);
    %sinusoidal tagging signal
    freqTable(f,:) = cfg.patch_amplitude * sin(patch_angFreq * frametime + cfg.patch_startPhase(f)) + cfg.patch_amplitude + cfg.f_offset;
end
cfg.freqTable = freqTable;

%%%%%%%%%%%%%%%%%================= output setting ======================%%%%%%%%%%%%%%%%%%%%%
Result.FixationON = nan(nTrials,1);
Result.StartBoxON = nan(nTrials,1);
Result.SentenceON = nan(nTrials,1);
Result.ITION = nan(nTrials,1);
Result.ProbeON = nan(nTrials,1);
Result.PureTagON = nan(nTrials,1);
Result.RT = nan(nTrials,1);
Result.CORR = nan(nTrials,1);
Result.KeyPress = cell(nTrials,1);
Result.WordLocation = zeros(nTrials,size(SentMat,2),4); %[x-start, x-end, y-start y-end] unit in pixel of the small rapid window

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%===================== EXPERIMENTAL LOOP =======================%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% introduction#
KbReleaseWait; 
if ~cfg.debugmode
   KbQueueFlush();
end
% set up Text parameters
Screen('TextFont', w, cfg.TextFont);
Screen('TextSize', w, cfg.TextSize);
Screen('TextStyle',w, cfg.TextStyle);
%welcome screen at first trial
message = 'Please make sure the back and the top of \n\n your head touch the helmet tightly. \n\n\n\n Press any button to start!';
for q = 1:4
    DrawFormattedText(w,message,'center','center',cfg.TextColor,[],[],[],[],[],q_rects(q,:));
end
Screen('Flip', w);
KbWait; %% waiting for key pressing

%%%%%%%%%%%%%================ Eye link: calibration + validation ===============%%%%%%%%%%%%%%%%%
%%%%% testing eyetracker
if cfg.el.eyelink
    cfg = el_calib_valid(cfg,0); % get the cfg.el.defaults settings
end
if cfg.el.eyelink
    %%Experiment start message to eyelink
    Eyelink('Message', 'Exp start');
end

%%%% no trigger sent yet
cfg.triggerSent=0;

%%%%%%%%%%%%%%================  trial loops  ===============%%%%%%%%%%%%%%%%%
for i = 1:nTrials
    %%%  releasing both key press and button pressq
    KbReleaseWait;
    if ~cfg.debugmode
        KbQueueFlush();
    end
    
    %%%=== define frames for each trial
    Para.FixDuration(i,1) = cfg.fix_t + cfg.fix_jitter*rand;
    fixendframes = round(Para.FixDuration(i,1)/ifi);
    
    %%%set frequency table
    lum = freqTable;
    
    %%% do drift correction for every 3 trails 
    if  cfg.el.eyelink && i ~= 1
        if mod(i,3)==1 || mod(i,Para.BreakTrials)==1
            el_calib_valid(cfg,2); % run EyelinkDoDriftCorrection
        end
    end
    Screen('Flip', w);
    WaitSecs(0.001);
    
    %%%%%%%%%%%%%%%%%%%% ============ 1. fixation ============ %%%%%%%%%%%%%%%%%%%%
    for j = 1:fixendframes
        % draw cross-fixation in all 4 quadrants
        for q = 1:4
            Screen('DrawLines', w, CrossCoords, lineWidthPix, cfg.TextColor, qcenters(q,:), 2);
        end
        % draw photodiode in all 4 quadrants
        DrawPhotodiode(w,diode_tex,diode_pos,lum,j)
        % Flip to the screen
        vbl = Screen('Flip', w);
        %%%=== send triggers
        if j == 1
            Result.FixationON(i,1) = vbl;
            cfg = sendTrigger(cfg,cfg.TriggerFix);
        end
        %we want to reset the trigger 50ms after the last trigger
        if cfg.triggerSent && GetSecs>(cfg.triggerTime+0.05) && ~cfg.debugmode
            io64(cfg.ioObjTrig,cfg.PortAddress,0);
            cfg.triggerSent=0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%% ============ 2. start box ============ %%%%%%%%%%%%%%%%%%%%
    ppp = 1; %%% index to send trigger
    j = 1; %% index to change the frequency table in each frame
    f = 1; %index for the frames
    while f <= round(cfg.gaze_start/ifi)
        %%% draw start box
        Screen('FillRect', w, cfg.TextColor, DotCoords_start');
        % draw photodiode in all 4 quadrants
        DrawPhotodiode(w,diode_tex,diode_pos,lum,j);
        % Flip to the screen
        vbl = Screen('Flip', w);
        
        %%%=== send triggers
        if ppp == 1
            Result.StartBoxON(i,1) = vbl;
            cfg = sendTrigger(cfg,cfg.TriggerStartBox);
            ppp = 0;
        end
        %we want to reset the trigger 50ms after the last trigger
        if cfg.triggerSent && GetSecs>(cfg.triggerTime+0.05) && ~cfg.debugmode
            io64(cfg.ioObjTrig,cfg.PortAddress,0);
            cfg.triggerSent=0;
        end
        
        %check whether fixation is in the start box
        if cfg.el.eyelink
            sample = Eyelink('NewestFloatSample');
            % Get current gaze position from sample
            x = sample.gx(eyeused_id); %first sample should be left eye
            y = sample.gy(eyeused_id);
            %%% compare x and y to start box window
            if sample.pa(eyeused_id)>0
                inbox = (x>cfg.el.startWindow(1) && x<cfg.el.startWindow(3) && y>cfg.el.startWindow(2) && y<cfg.el.startWindow(4));
                if inbox == 0
                    f = 0;
                end
            else
                f = 0;
            end
        else
            f = 0;
        end
        
        %%%% check key press
        [keyIsDown, ~, KeyCode ] = KbCheck;
        if keyIsDown && KeyCode(escKey)
            cleanup(cfg);
            break;
        end
        if keyIsDown && KeyCode(cfg.el.eyelinkKey)
            el_calib_valid(cfg,2)  %% do drift-checking
            f = 0;
        end
        if keyIsDown && KeyCode(cfg.el.SkipStartboxKey)
            break
        end
        j = j + 1;
        f = f + 1;
    end
    
    %%%%%%%%%%%%%%%%%%%% ============ 3. sentence  ============ %%%%%%%%%%%%%%%%%%%%
    %%%%trial start
    if cfg.el.eyelink
        %send trial start trigger to eyelink
        Eyelink('Message', ['Sentence_' int2str(i) ' start' ]);
    end
    %%% get the sentence and its x-coordinate in this trial
    Words = SentMat(i,:); %% words in the present sentences
    for wid = 1:length(Words)
        if strfind(Words{wid},'.')
            break
        end
    end
    NumWord = wid;
    TxtXcoord = zeros(1,NumWord); %% the x-coordinate of each word
    WordWid = zeros(1,NumWord);
    WordHit = zeros(1,NumWord);
    for wid = 1:NumWord
        textmp = Words{wid};
        WordWid(wid) = RectWidth(Screen('TextBounds',w,textmp));
        WordHit(wid) = RectHeight(Screen('TextBounds',w,textmp));
        if wid ~= 1
            TxtXcoord(wid) = TxtXcoord(wid-1)+WordWid(wid-1)+WordSpace;
        end
    end
    TxtXcoord = TxtXcoord + WordStart;
    MostHeight = mode(WordHit);
    % Compute bounding box of textstring:
    TextureBox1 = q_rects(1,:);
    
    %%% recording eye-movement coordinates
    Result.WordLocation(i,1:NumWord,1) = TxtXcoord; %% x-start
    Result.WordLocation(i,1:NumWord,3) = TxtXcoord+WordWid; %% x-end
    ys = 0.5*TextureBox1(4)-1.25.*MostHeight;
    ye = 0.5*TextureBox1(4)-.25.*MostHeight;
    Result.WordLocation(i,1:NumWord,2) = ys.*ones(1,NumWord);
    Result.WordLocation(i,1:NumWord,4) = ye.*ones(1,NumWord);
    
    %%% make an offscreen to draw the sentence
    sents = Screen('OpenOffscreenwindow', w, [cfg.ScrBgc 0],TextureBox1);
    % need to set up the text parameters again here
    Screen('TextFont', sents, cfg.TextFont);
    Screen('TextSize', sents, cfg.TextSize);
    Screen('TextStyle',sents, cfg.TextStyle);
    for wid = 1:NumWord
        Screen('DrawText', sents, Words{wid},TxtXcoord(wid),0.5*TextureBox1(4)+0.25*MostHeight,[1 1 1], [], 1);
    end
    
    %%% flickering words, get [texture_tagwrd, texture_mask, destination_rect4both]
    % flickering the target word1
    [tag_1,mask_1,dest_1] = DrawFlicker(TxtXcoord,FlkWords(i,1),cfg,TextureBox1,WordSpace);
    % flickering the post-target word1
    [tag_2, mask_2,dest_2] = DrawFlicker(TxtXcoord,FlkWords(i,1)+1,cfg,TextureBox1,WordSpace);
    tag_all = [tag_1 tag_2];
    mask_all = [mask_1 mask_2];
    dest_all = [dest_1 dest_2];
    if ~isnan(FlkWords(i,2)) %need to tag two target words and two postarget words (for sentences from Degno JEP 2019)
        % flickering the target word2
        [tag_3, mask_3,dest_3] = DrawFlicker(TxtXcoord,FlkWords(i,2),cfg,TextureBox1,WordSpace);
        % flickering the post-target word2
        [tag_4, mask_4,dest_4] = DrawFlicker(TxtXcoord,FlkWords(i,2)+1,cfg,TextureBox1,WordSpace);
        tag_all = [tag_all tag_3 tag_4];
        mask_all = [mask_all mask_3 mask_4];
        dest_all = [dest_all dest_3 dest_4];
    end
    
    %%%%%========= frames loops in each trial ===========%%%%%%%
    KbReleaseWait;
    if ~cfg.debugmode
        KbQueueFlush();
    end
    ppp = 1; %%% index to send trigger
    j = 1; %% index of the frames to change the frequency table in each frame
    f = 1;
    while f <= round(cfg.gaze_end/ifi) %% word frames during one trial
        % put all the textures together, texture order matters-->layer order of texture
        texts = [tag_all mask_all sents diode_tex diode_tex];
        % multiple offscreen textures to four quadrants
        texts = kron(texts,ones(1,4));
        % get all the destination rect coordinates [4rows*n_column]
        dests = [dest_all dest_all q_rects' diode_pos];
        col_tag1 = [lum(1,(j-1)*12+(1:4)); lum(1,(j-1)*12+(5:8)); lum(1,(j-1)*12+(9:12))];% 3row(RGB)*4column(quadrant)
        col_tag2 = [lum(2,(j-1)*12+(1:4)); lum(2,(j-1)*12+(5:8)); lum(2,(j-1)*12+(9:12))];% 3row(RGB)*4column(quadrant)
        col_wrd = repmat(cfg.TextColor,3,4);
        col_msk = ones(3,4).*white;
        % put color to [tag_1 tag_2 mask_1 mask_2]
        colors = [col_tag1 col_tag2 col_msk col_msk ];
        if ~isnan(FlkWords(i,2)) %need to tag two target words and two postarget words
            % put color to [tag_1 tag_2 tag_3 tag_4 mask_1 mask_2 mask_3 mask_4]
            colors = [col_tag1 col_tag2 col_tag1 col_tag2 col_msk col_msk col_msk col_msk];
        end
        % get the whole color matrix [flickeringwords masks sents diode_1 diode_2]
        colors = [colors col_wrd col_tag1 col_tag2];
        Screen('DrawTextures', w, texts,[], dests,0, [], 1, colors);
        Screen('FillRect', w, endbox_color, DotCoords_end');
        %%% flip the frame
        [vbl] = Screen('Flip', w, vbl + 0.5 * ifi);
        
        %%%=== send triggers
        if ppp == 1
            Result.SentenceON(i,1) = vbl;
            cfg = sendTrigger(cfg,cfg.TriggerSentOn);
            ppp = 0;
        end
        %we want to reset the trigger 50ms after the last trigger
        if cfg.triggerSent && GetSecs>(cfg.triggerTime+0.05) && ~cfg.debugmode
            io64(cfg.ioObjTrig,cfg.PortAddress,0);
            cfg.triggerSent=0;
        end
        
        %%% check end box fixation
        if cfg.el.eyelink
            %Get eyelink sample
            sample = Eyelink('NewestFloatSample');
            % Get current gaze position from sample
            x = sample.gx(eyeused_id);
            y = sample.gy(eyeused_id);
            %%%% compare x and y to fixation_end_window
            if sample.pa(eyeused_id)>0
                inbox = x > cfg.el.endWindow(1) &&  x <  cfg.el.endWindow(3) && y > cfg.el.endWindow(2) && y < cfg.el.endWindow(4);
                if inbox == 0
                    f = 0;
                end
            else
                f = 0;
            end
        else
            f = 0;
        end
        
        %CHECK RESPONSES
        [keyIsDown, ~, KeyCode] = KbCheck;
        if keyIsDown && KeyCode(escKey)
            cleanup(cfg);
            break;
        end
        %%%%check the eye-tracker key-press
        if keyIsDown && KeyCode(cfg.el.eyelinkKey)
            el_calib_valid(cfg,2)  %% do drift-checking
            f = 0;
        end
        if keyIsDown && KeyCode(cfg.el.SkipEndboxKey)
            break
        end
        j = j + 1;
        f = f + 1;
    end
   %%% sentence off trigger
    cfg = sendTrigger(cfg,cfg.TriggerSentOff);
    if cfg.triggerSent && GetSecs>(cfg.triggerTime+0.05) && ~cfg.debugmode
        io64(cfg.ioObjTrig,cfg.PortAddress,0);
        cfg.triggerSent=0;
    end
    %%% close offscreens/textures
    Screen('Close', [tag_all mask_all sents]);
    
    %%%%%%%%%%%%%%%%%%%% ============ 4. probe questions  ============ %%%%%%%%%%%%%%%%%%%%
    KbReleaseWait;
    if ~cfg.debugmode
        KbQueueFlush();
    end
    answer = Question{i,2};
    if  strfind('TF',answer)
        Screen('TextFont', w, cfg.TextFont);
        Screen('TextSize', w, cfg.TextSize);
        Screen('TextStyle',w, cfg.TextStyle);
        message = [Question{i,1} '\n\n\n\n\n ' fingers{1} ' \n\n ' fingers{2}];
        for q = 1:4
            DrawFormattedText(w,message,'center','center',cfg.TextColor,[],[],[],[],[],q_rects(q,:));
        end
        [vbl] = Screen('Flip', w);
        Result.ProbeON(i,1) = vbl;
        
        %%% check the response
        noResponse = 1;
        while noResponse
            if ~cfg.debugmode %% NaTA box
                [pressed, firstpress] = KbQueueCheck(); %check response, return whether pressed, and first press timestamp
                KbQueueFlush(); %only the KbQueues events are deleted.
                %Note that two keys may have been pressed
                KeyCode = find(firstpress);
                if length(KeyCode)>1 %two or more buttons pressed
                    [~,ind] = min(firstpress(KeyCode));
                    KeyCode = KeyCode(ind); %select first response
                end
                t_keypress = firstpress(KeyCode);
            else   %% keyboard
                [pressed, t_keypress, KeyCode] = KbCheck;
            end
            if pressed
                if ~cfg.debugmode  %% NaTA box
                    if (strcmp(Question{i,2},'T') && KeyCode == cfg.keyTrue) || (strcmp(Question{i,2},'F') && KeyCode == cfg.keyFalse)
                        Result.CORR(i) = 1;
                    else
                        Result.CORR(i) = 0;
                    end
                else  %% keyboard
                    if (strcmp(Question{i,2},'T') && KeyCode(cfg.keyTrue)) || (strcmp(Question{i,2},'F') && KeyCode(cfg.keyFalse))
                        Result.CORR(i) = 1;
                    else
                        Result.CORR(i) = 0;
                    end
                    KbReleaseWait;
                end
                Result.RT(i,1) = t_keypress-vbl;
                Result.KeyPress{i,1} = KbName(KeyCode);
                noResponse=0;
            end
            WaitSecs(0.001);
        end
        
        %%%%%%%%%%%%%%%%%%%% ============ 5. feedback to the resposne  ============ %%%%%%%%%%%%%%%%%%%%
        if Result.CORR(i)
            message = 'Correct';
        else
            message = 'Wrong';
        end
        for q = 1:4
            DrawFormattedText(w,message,'center','center',cfg.TextColor,[],[],[],[],[],q_rects(q,:));
        end
        Screen('Flip', w);
        WaitSecs(cfg.fedbk);
    end
    
    %%%%%%%%%%%%%%%%%%%% ============ 6. blank screen  ============ %%%%%%%%%%%%%%%%%%%%
    vbl = Screen('Flip', w);
    Result.ITION(i,1) = vbl;
    if ~cfg.debugmode
        % sent trigger
        cfg = sendTrigger(cfg,cfg.TriggerITI);
        WaitSecs(0.05);
        io64(cfg.ioObjTrig,cfg.PortAddress,0);
        cfg.triggerSent=0;
    end
    WaitSecs(cfg.iti - 0.05);
    
    %%%%%%%%%%%%%%%%%%%% ============ 7. rest  ============ %%%%%%%%%%%%%%%%%%%%
    %ADD extra break after every so many trials
    if ~mod(i,Para.BreakTrials) && i ~= nTrials
        %%% have a rest
        for kkkk = 1:cfg.rest
            message=['You have done ' num2str(i/Para.BreakTrials) ' out of ' num2str(n_block) ' blocks!'...
                '\n\n Please close your eyes and take a break \n\n\n\n ' num2str(cfg.rest-kkkk)];
            for q = 1:4
                DrawFormattedText(w,message,'center','center',cfg.TextColor,[],[],[],[],[],q_rects(q,:));
            end
            Screen('Flip', w);
            WaitSecs(1);
        end
        
        %%% press any button to continue
        message = 'Please make sure the back and the top of \n\n your head touch the helmet tightly. \n\n\n\n Press any button to start!';
        for q = 1:4
            DrawFormattedText(w, message, 'center', 'center',cfg.TextColor,[],[],[],[],[],q_rects(q,:));
        end
        Screen('Flip', w);
        %%% check button press
        if ~cfg.debugmode
            KbQueueFlush();
        end
        KbReleaseWait;
        noResponse = 1;
        while noResponse
            if ~cfg.debugmode
                [pressed] = KbQueueCheck(); %check response, return whether pressed, and first press timestamp
            else
                [pressed] = KbCheck;
            end
            if pressed
                noResponse = 0;
            end
        end
        %add empty screen
        Screen('Flip', w);
    end
    %%% saving data
    save(datafilename,'cfg','Para','Result');
end

%%% end of study
message = 'The End! \n\n Well done! \n\n  THANK YOU !!';
for q = 1: 4
    DrawFormattedText(w, message, 'center', 'center',cfg.TextColor,[],[],[],[],[],q_rects(q,:));
end
Screen('Flip', w);
WaitSecs(1);

%stop eyelink & transfer file
if cfg.el.eyelink
    Eyelink('Message', 'end of block');
    el_Stop(cfg);
end

%return to lower priority
if ~cfg.debugmode
    Priority(0);
end
%ListenChar(0);
ShowCursor;
%set propixx to normal state
if ~cfg.debugmode || cfg.DataPixxOnly
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
end
%Close screen
Screen('CloseAll');
end

%function to send MEG and eyelink triggers
function [cfg] = sendTrigger(cfg,trig)
%send trigger to MEG
if ~cfg.debugmode
    io64(cfg.ioObjTrig,cfg.PortAddress,trig);
    cfg.triggerSent = 1;
end
%send trigger to eyelink
if cfg.el.eyelink
    Eyelink('Message', ['Trigger_' int2str(trig)]);
end
cfg.triggerTime = GetSecs;
end

function [] = cleanup(cfg)
%Return propixx to normal state
if  ~cfg.debugmode
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
end
%lower priority
if ~cfg.debugmode
    Priority(0);
end
%stop eyelinkq
if cfg.el.eyelink
    Eyelink('Message', 'end of block - ABORTED');
    el_Stop(cfg);
end
%close screen
Screen('CloseAll');
%ListenChar(0);
ShowCursor;
%throw warning due to prematurely aborted experiment56
warning('Experiment aborted');
end

function pix = Deg2Pix_RIFT(degree,cfg)
%%% in rapid mode, stimuli are presented at four quarands
pix = tand(degree)*cfg.dist*0.5*cfg.resy/cfg.height;
pix = round(pix);
end

function pix = Deg2Pix(degree,cfg) 
%%% in normal mode, stimuli are presented at the screen centre
pix = tand(degree)*cfg.dist*cfg.resy/cfg.height;
end

function cfg = el_calib_valid(cfg,mode)
%%% set screen back to normal cfg.RT
Datapixx('SetPropixxDlpSequenceProgram', 0);
Datapixx('RegWrRd');
if mode == 0 %% run the eye-tracker setup for the first time
    %%% run el_Start
    cfg = el_Start_SameWindow(cfg);
    cfg.el.defaults.ScrBgc = cfg.ScrBgc;
    cfg.el.defaults.TextColor = cfg.TextColor;
elseif mode == 1 %% run the eye-tracker setup for the non-first time
    %%% re-calibration and re-validation and re-drift correction
    EyelinkDoTrackerSetup(cfg.el.defaults);
elseif mode == 2
    EyelinkDoDriftCorrection(cfg.el.defaults);
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.01);
    % mark zero-plot time in data file
    Eyelink('Message', 'SYNCTIME');
end
%%% set screen to rapid mode
Datapixx('SetPropixxDlpSequenceProgram',cfg.ProPixxModel);
Datapixx('RegWrRd');
%%% set screen to experiment background color
Screen('FillRect', cfg.window, cfg.ScrBgc);
Screen('Flip', cfg.window);
end

function DrawPhotodiode(w,diode_tex,diode_pos,lum,j)
colors = [];
for i = 1: size(lum,1) %loop over number of tagfreq
    col_tag = [lum(i,(j-1)*12+(1:4)); lum(i,(j-1)*12+(5:8)); lum(i,(j-1)*12+(9:12))];
    colors = [colors col_tag]; % 3row*n columns
end
% draw photodiode in all 4 quadrants
texts_pd = kron(diode_tex,ones(1,size(colors,2))); %multiply textures into 4 quadrants
Screen('DrawTextures',w,texts_pd,[],diode_pos,0,[],1,colors);
end

function [woff_rect,woff_mask,dest_rect] = DrawFlicker(TxtXcoord, flk_id,cfg,TextureBox1,WordSpace)
% only tag the location of the word + its left space (no right space)
x_start = round(TxtXcoord(flk_id)-WordSpace);
x_end = round(TxtXcoord(flk_id+1)-WordSpace);
rectwidth = x_end - x_start;
%%% make it a square
rectheight = rectwidth;
y_start = round(0.5*TextureBox1(4)-rectheight/2);
y_end = y_start + rectheight;
woff_rect = Screen('OpenOffscreenwindow', cfg.window, [cfg.ScrBgc 0],[0 0 rectwidth rectheight]);
Screen('FillRect', woff_rect,[1 1 1]);
% We create a Luminance+Alpha matrix for use as transparency mask:
% Layer 1 (Luminance) is filled with luminance value 'gray' of the background.
ms = round(rectwidth)/2;
transLayer = 2;
%%%% square mask
msy = ms;%make it a square
[x,y] = meshgrid(-ms:ms, -msy:msy);
maskblob = uint8(ones(2*msy+1, 2*ms+1, transLayer) * 128);
% Layer 2 (Transparency aka Alpha) is filled with gaussian transparency mask.
xsd = ms/1.2; %bigger than 1.2, more concentrated blob, less smoothing area
ysd = msy/1.2;
maskblob(:,:,transLayer) = uint8(round(255 - exp(-((x/xsd).^2)-((y/ysd).^2))*255));
% Build a single transparency mask texture in offscreen 3
woff_mask = Screen('MakeTexture', cfg.window, maskblob);
wordcenter_x = (x_start + x_end)/2;
wordcenter_y = (y_start + y_end)/2;
dest_rect = [wordcenter_x-ms  wordcenter_y-msy  wordcenter_x+ms  wordcenter_y+msy];
% we need to add the starting points (x1,y1,x2,y2) in each quadrant to the
% dest_rect, which will be used as the destination rect during
% DrawTexture
xx = TextureBox1(3);
yy = TextureBox1(4);
plus_start = [0 0 0 0; xx 0 xx 0; 0 yy 0 yy; xx yy xx yy]; %[left up;right up;left bottom;right bottom]
dest_rect = plus_start + dest_rect;
dest_rect = dest_rect'; %transpose to 4row*n_column
end

