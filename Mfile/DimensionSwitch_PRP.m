function DimensionSwitch_PRP
%==========================================================================
% PURPOSE:
% -PRP task with 4 x 4 bivarient stimuli
% SETTING:
% -4 task sets: Odd-Even, High-Low / Red-Blue, Faint-Bold
% -block-wise rule: any 2 tasks come as a pair
%==========================================================================
clear all;
close all;
sca;
commandwindow;
if 1, Screen('Preference', 'SkipSyncTests', 1); end
%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST TASK_ST FPATH_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST CONFIG TRIGGER_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL PRACEXC BTRIAL BLOCKSTS
%**************************************************************************
% TASK SETTING
TASK_ST.NSTIM=2;
TASK_ST.STIM(1,:)={3,4,6,7};
TASK_ST.STIM(2,:)={'RED','LIGHTRED','LIGHTBLUE','BLUE'};
TASK_ST.STIMPAIRS=allcomb(1:size(TASK_ST.STIM,2), 1:size(TASK_ST.STIM,2));
TASK_ST.TASK(1,:)={'Low-High','Odd-Even'};% SRREF 1, SRREF2
TASK_ST.TASK(2,:)={'Red-Blue','Bold-Faint'};% SRREF 1, SRREF2
TASK_ST.TASKPAIRS=allcomb(1:numel(TASK_ST.TASK), 1:numel(TASK_ST.TASK));
TASK_ST.SRREF(1,:)={[1,1,2,2], [1,2,2,1]};
TASK_ST.SRREF(2,:)={[1,1,2,2], [1,2,2,1]};
TASK_ST.ERRORCUT=1.00;%always go on!

%EXPERIMENT SETTING
EXC_ST.NUMPBLOCK=1;
EXC_ST.NUMEBLOCK=150;
EXC_ST.NUMTBLOCK=EXC_ST.NUMPBLOCK+EXC_ST.NUMEBLOCK;
EXC_ST.BLOCKDUR_P=15.00;% in seconds
EXC_ST.BLOCKDUR_E=20.00;% in seconds
EXC_ST.MAXTRIAL=35;%max # of trials to prepare data(keep it as even number)!
EXC_ST.BREAKBLOCK=EXC_ST.NUMPBLOCK;

%INCENTIVES/POINTS
EXC_ST.ACCCUT=80;%at least n% accuracy
EXC_ST.RESPCUT=6;%at least n responses
EXC_ST.RTCUT=70;%RTs from the current trial are below 70th percentile

%TIME SETTING
TIME_ST.YEILDT=0.002;
TIME_ST.TRIALT_P=inf;
TIME_ST.RSIT=0.300;
TIME_ST.FBT=0.100;
TIME_ST.SOA=[0.10, 0.25, 0.50, 0.75, 1.00];
TIME_ST.SDISP=.200;

%FILE PATH SETTING
FPATH_ST.TXTCODE='dsPRP';
FPATH_ST.EXFILE='DimensionSwitch_PRP.m';
FPATH_ST.PTLPATH=PsychtoolboxRoot;
exFilePath=which(FPATH_ST.EXFILE);
FPATH_ST.EXDIR=exFilePath(1:end-length(FPATH_ST.EXFILE));
FPATH_ST.SPDIR=[FPATH_ST.EXDIR,'z_Support'];
FPATH_ST.DATADIR=[FPATH_ST.EXDIR,'z_Data'];
cd(FPATH_ST.EXDIR);
FPATH_ST.PORTID_S={'BCC8'}; %Oregon (EEG_EC)
FPATH_ST.PORTID_R={'BCD8'}; %Oregon (EEG_EC)
FPATH_ST.IO64ID={hex2dec('D050')}; %Oregon (BP)
FPATH_ST.PORDIR=fileparts(which('write_parallel.m'));
FPATH_ST.IO64DIR=fileparts(which('write_parallel_IO64.m'));
for d={'EXDIR','SPDIR','DATADIR','PORDIR','IO64DIR'},addpath(FPATH_ST.(d{:}));end

%STIMULUS SETTING
STIM_ST.RESOLUTION=[1920,1080];% Oregon = 1024,768, RIKEN = 1920, 1080 
STIM_ST.FONTSIZE=30;
STIM_ST.FONTSTIM=100;
STIM_ST.FONTFC=55;
STIM_ST.PLSIZE=80;
STIM_ST.FRSIZE=7;
STIM_ST.FROFFSET=[-150, 150];
STIM_ST.FRRECT=[0 0 200 200];
STIM_ST.STIMRECT=ceil([0 0 60 60]);
STIM_ST.INSTRUCTION=floor([0,0,1500,900]/2.2);
if ismac,STIM_ST.FXADJ=10;else,STIM_ST.FXADJ=10;end

%COLOR SETTING
COLOR_ST.BLACK = [0,0,0];
COLOR_ST.GRAY = [165,165,165];
COLOR_ST.WHITE = [250,250,250];
COLOR_ST.RED = [255,0,0];
COLOR_ST.LIGHTRED = [250,192,203];
COLOR_ST.BLUE = [0,0,205];
COLOR_ST.LIGHTBLUE = [135,206,250];
COLOR_ST.GREEN = [0,250,0];
COLOR_ST.PURPLE = [155,0,250];
COLOR_ST.YELLOW=[250,250,0];
COLOR_ST.ORANGE = [250,165,0];
COLOR_ST.DARKBLUE = [50 90 170];
COLOR_ST.BACKG = COLOR_ST.GRAY;

%KEY SETTING
KbName('UnifyKeyNames');
% % AT RIKEN
% KEY_ST.LKEY1 = KbName('4');
% KEY_ST.RKEY1 = KbName('5');
% KEY_ST.LKEY2 = KbName('1');
% KEY_ST.RKEY2 = KbName('2');
% Testing on MBP
KEY_ST.LKEY1 = KbName('a');
KEY_ST.RKEY1 = KbName('s');
KEY_ST.LKEY2 = KbName('z');
KEY_ST.RKEY2 = KbName('x');
KEY_ST.REFKEY=[KEY_ST.LKEY1,KEY_ST.RKEY1;KEY_ST.LKEY2,KEY_ST.RKEY2];
KEY_ST.SPACEKEY = KbName('SPACE');
KEY_ST.MAGICKEY = KbName('');
KEY_ST.ESCKEY = KbName('ESCAPE');
keyList = zeros(1, 256);% [~, ~, KEY_ST.KEYLIST] = KbCheck;
keyList(1, [KEY_ST.ESCKEY, reshape(KEY_ST.REFKEY,1,[])]) = 1;
KEY_ST.KEYLIST = keyList;
[KEY_ST.KEYINDEX, ~, ~] = GetKeyboardIndices;

%Get subject ID and open data file for output==============================
EXC_ST.SUBID = input('Enter Subject ID: ','s');
EXC_ST.SESSION = input('Enter Session ID: ','s');
EXC_ST.EXMODE = str2double(input('EEG_BP(2) / EEG_EC(X) / BEH(0): ','s'));
EXC_ST.FNAME = [EXC_ST.SUBID,'_',EXC_ST.SESSION,'_',FPATH_ST.TXTCODE];
RECORD_ST.SUBID=str2double(EXC_ST.SUBID);
RECORD_ST.SESSION=str2double(EXC_ST.SESSION);
RECORD_ST.DFILE = [FPATH_ST.DATADIR,filesep,EXC_ST.FNAME,'_RECORD_ST.mat'];
RECORD_ST.DLOADED = logical(exist(RECORD_ST.DFILE,'file'));
if RECORD_ST.DLOADED,load(RECORD_ST.DFILE);end

%Configure serial port (IOPort or IO64) ===================================
% Configure serial port
if EXC_ST.EXMODE==2,CONFIG = IOPort('OpenSerialPort', 'COM3');end 

% http://apps.usd.edu/coglab/psyc770/IO64.html
% if EXC_ST.EXMODE==2,cd(FPATH_ST.IO64DIR{:});CONFIG=config_io;cd(FPATH_ST.EXDIR);end

%Screen Stuff==============================================================
STIM_ST.SCCODE=max(Screen('Screens'));
% if ~ismac,Screen('Resolution',STIM_ST.SCCODE,STIM_ST.RESOLUTION(1),STIM_ST.RESOLUTION(2));end;% adjust resolution?   
% Screen('Resolution',STIM_ST.SCCODE,STIM_ST.RESOLUTION(1),STIM_ST.RESOLUTION(2));% adjust resolution?   
[W, wRect]=Screen('OpenWindow',STIM_ST.SCCODE, 0,[],32,2);%[0 0 200 200]
Screen('FillRect',W,COLOR_ST.BACKG);
[swidth, sheight]=Screen('WindowSize',W);
STIM_ST.XC=fix(swidth/2);
STIM_ST.YC=fix(sheight/2);
TIME_ST.IFI = Screen('GetFlipInterval', W);
Screen('TextFont',W,'Arial');
Screen('TextSize',W,STIM_ST.FONTSIZE);

%Setting up randomized array for a pracblock===============================
rng(str2double(EXC_ST.SUBID),'twister');
ListenChar(-1);
HideCursor;
locSetUpEx;
locLoadImage;
tic;GetSecs;toc;% doing this once makes later GetSecs faster...?
%YAY!

%Main loop structure=======================================================
for BLOCK=[1:EXC_ST.NUMTBLOCK]
    PRACEXC=BLOCK>EXC_ST.NUMPBLOCK;
    locShowMessege('BINTRO');
    if ~PRACEXC,BLDUR=EXC_ST.BLOCKDUR_P;else, BLDUR=EXC_ST.BLOCKDUR_E;end
        
    % Mark the begining of each block
    locCodeEvent('TEMP',0);%null code! (to prebent mis-trigegr)
    locCodeEvent('TEMP',10);%begining of block @@@@@@@@@@@@@@@@@@@@@@@@@@@@

    % Photosensor trigger for ECoG: flip black mark to white!
    TRIAL=0;
    BLOCKSTS=GetSecs;
    RECORD_ST.EV_BLOCK(BLOCK,:)=BLOCKSTS;%record timing for re-epoching!
    
    %Repeatedly loop trials for the block duration!!
    while GetSecs-BLOCKSTS<BLDUR
        %Within trial computations
        TRIAL=TRIAL+1;
        RECORD_ST.EV_TRIAL(BLOCK,TRIAL)=GetSecs-BLOCKSTS;
        locDoTrial;
        
        % Record each trial?
        %cd(FPATH_ST.DATADIR);
        %save([EXC_ST.FNAME '_RECORD_ST.mat'],'RECORD_ST');
                
        %Make sure that trial does not exceed the limit
        if TRIAL== EXC_ST.MAXTRIAL,break;end;
        
        %Exit out if esc is pressed
        [~,~,keyCode]=KbCheck;
        if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;
    end

    %Block result 
    if BLOCK <= EXC_ST.NUMTBLOCK,locShowMessege('BRESULT');end;
    
    %Save RECORD_ST??
    cd(FPATH_ST.DATADIR);
    save([EXC_ST.FNAME '_RECORD_ST.mat'],'RECORD_ST');
    
    %Breask time Screen
    %if any(ismember(EXC_ST.BREAKBLOCK,BLOCK)),locShowMessege('EBREAK');end;

    %End Screen
    if BLOCK == EXC_ST.NUMTBLOCK,locShowMessege('EXEND');end;
end

%Putting everything back
cd(FPATH_ST.EXDIR);
Screen('CloseAll');
ShowCursor;
ListenChar;
if EXC_ST.EXMODE==2,IOPort('Close', CONFIG);end

end

%##########################################################################

function locSetUpEx
%locSetUpEx

%**************************************************************************
%BLOCK ORGANIZATION/%TIME REGULATION/%STIMULUS/%KEYS/%FONT/%COLOR/%PARAMETER
global EXC_ST TASK_ST FPATH_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST
%**************************************************************************

%PARAMST SETTING===========================================================
PARAM_ST.PARAMETER={
    'CONTEXT';%set up task context parameters
    'CONTENTS';%code parameters coding contents and position of 2-levels
    'TASK';%position of stimulus (1=top_l,2=top_r,3-bot_r,4=bot_l)
    };

%RANDOMIZATION LOOP========================================================
%Parameter Sign
for p=1:length(PARAM_ST.PARAMETER)
    %Parameter Setting!
    fn=PARAM_ST.PARAMETER{p,1};
    numPB=EXC_ST.NUMPBLOCK;
    numEB=EXC_ST.NUMEBLOCK;numET=EXC_ST.MAXTRIAL;
    numTB=EXC_ST.NUMTBLOCK;
    
    switch fn
        case {'PRACEX'}
        case {'CONTEXT'}
            %TASK RULE PAIRS
            taskpairs = RandSample(1:size(TASK_ST.TASKPAIRS, 1), [numTB, 1]);
            PARAM_ST.('TASKPAIR') = repmat(taskpairs,[1, numET]);
            
            % TASK SETS
            tasksets = num2cell(TASK_ST.TASKPAIRS(taskpairs,:),2);
            PARAM_ST.('TASKSET') = repmat(tasksets, [1, numET]);
            
            % RELEVANT FEATURE DIM
            fdim = cellfun(@(x) (x> size(TASK_ST.TASK,1))+1, tasksets, 'Uni', false);
            PARAM_ST.('FDIM') = repmat(fdim, [1, numET]);
           
        case {'TASK'}
            % STIMPOS
            stimpos = repmat({1:TASK_ST.NSTIM},numTB,numET);
            PARAM_ST.('STIMPOS') = cellfun(@Shuffle,stimpos,'Uni',false);
            
            % SOA condition
            PARAM_ST.('SOA') = RandSample(1:length(TIME_ST.SOA), [numTB, numET]);
            
            % Sample stim features randomly
            nfeat = size(TASK_ST.STIM, 2);
            s1_f1 = RandSample(1:nfeat,[numTB,numET]);
            s1_f2 = RandSample(1:nfeat,[numTB,numET]);
            s2_f1 = RandSample(1:nfeat,[numTB,numET]);
            s2_f2 = RandSample(1:nfeat,[numTB,numET]);
            
            % STIM PROPERTY ({1st object, 2nd object})
            bindF = @(x1,x2) cellfun(@(x)squeeze(x)',num2cell(cat(3, x1, x2), 3),'Uni',false);
            divF = @(xc,i) cell2mat(cellfun(@(x)x(i),xc,'Uni',false));
            PARAM_ST.('F_NUMBER') = bindF(s1_f1, s2_f1);% feature1(number): 1st and 2nd item 
            PARAM_ST.('F_COLOR') = bindF(s1_f2, s2_f2);% feature2(color): 1st and 2nd item 
            
            % STIM ID ({NUMBER, COLOR})
            PARAM_ST.('S1') = bindF(s1_f1, s1_f2);% feature1(number)-feature2(color) of 1st item
            PARAM_ST.('S2') = bindF(s2_f1, s2_f2);% feature1(number)-feature2(color) of 2nd item
            PARAM_ST.('S1ID') = (s1_f1 -1) * nfeat + s1_f2;
            PARAM_ST.('S2ID') = (s2_f1 -1) * nfeat + s2_f2;
            
            % SR (depends on taskset, feature dimension, and stim)
            % CORRESP (depends on stimpos and SR): code this locally!
            t = PARAM_ST.TASKSET;
            fd = PARAM_ST.FDIM;
            sr = TASK_ST.SRREF'; %'Low-High','Even-Odd','Red-Blue','Bold-Faint'
            c1 = cellfun(@(t, f, s) sr{t(1)}(s(f(1))), t, fd , PARAM_ST.S1, 'Uni', false);
            c2 = cellfun(@(t, f, s) sr{t(2)}(s(f(2))), t, fd , PARAM_ST.S2, 'Uni', false);
            PARAM_ST.('SR') = bindF(cell2mat(c1) ,cell2mat(c2));   
    end
end

%RECORD SETTING============================================================
cellD=repmat({zeros(1, TASK_ST.NSTIM)},numTB, numET);zeroD=zeros(numTB,numET);
zVars={'SUBID','BLOCK','TRIAL','RSIJIT','EV_BLOCK','EV_TRIAL','MONEY'};
cVars={'TASKTS','RT','RESP','CORRESP','RESPORDER','ACC','EV_STIM','EV_RESP'};

if ~RECORD_ST.DLOADED    
    for p=zVars(2:end),if ~isfield(RECORD_ST,p),RECORD_ST.(p{:})=zeroD;end;end
    for p=cVars,if ~isfield(RECORD_ST,p),RECORD_ST.(p{:})=cellD;end;end
    RECORD_ST.PARAMORDER=[zVars,cVars];
    RECORD_ST.SUBID=str2double(EXC_ST.SUBID);
end

%RECORD PARAMST JUST IN CASE!!=============================================
cd(FPATH_ST.DATADIR);
save([EXC_ST.FNAME '_PARAM_ST.mat'],'PARAM_ST');
save([EXC_ST.FNAME '_RECORD_ST.mat'],'RECORD_ST');
disp('YAY! Completed!!');

% sca
% keyboard

end

%##########################################################################

function locDoTrial
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPATH_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BLOCKSTS
%**************************************************************************

%PHASE0:Trial Preparation--------------------------------------------------
[~,~, keyCode]=KbCheck;
if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;

%PHASE0:Neuralize screen---------------------------------------------------
locDrawFC('NEUTRAL');
Screen('Flip',W);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE1:Fixed RSI----------------------------------------------------------
locCodeEvent('TEMP',9);%begining of trial @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
locCodeEvent('CHAR',[]);
if TRIAL ==1,ads=1.5;else, ads=1;end
WaitSecs(TIME_ST.RSIT * ads);
RECORD_ST.RSIJIT(BLOCK,TRIAL)=TIME_ST.RSIT*1000;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE1:STIMULUS DISPLAY & RECORDING---------------------------------------
%NO LOCAL FUNCTIONS ALLOWED TO OPTIMIZE TIMING CONTROL!!
spos = PARAM_ST.STIMPOS{BLOCK, TRIAL};%first-second
sr = PARAM_ST.SR{BLOCK, TRIAL};% first-second
for s = 1:numel(spos), coresp(s) = (KEY_ST.REFKEY(spos(s), sr(s))); end

% Initialize
phase = 1;% each phase corresponds to ON1-ON2-OFF1-OFF2
dispT = [0, TIME_ST.SOA(PARAM_ST.SOA(BLOCK, TRIAL))];
evT = [dispT, dispT+TIME_ST.SDISP];[~, evI] = sort(evT);%ON1-ON2-OFF1-OFF2
stateV = [false, false]; % [stim1, stim2] true = on, false = off
schedule = [evT(1),NaN;NaN,evT(2); evT(3),NaN;NaN, evT(4)];%ON1-ON2-OFF1-OFF2
schedule = schedule([evI],:);% sequence of phases: (1,:) = S1, (2,:) = S2
stateR = [0, 0];% [keyset ID, count]
refRT = [0, 0]; % reference time for RTs 

% Start loop
startT = GetSecs;
KbQueueCreate(KEY_ST.KEYINDEX, KEY_ST.KEYLIST);
KbQueueStart(KEY_ST.KEYINDEX);

% Detect 2 responses
while stateR(2) < TASK_ST.NSTIM || phase < size(schedule,1)
    % Always measure time and detect responses!
    curTime = (GetSecs - startT);
    %[keyP, secs, keyC] = KbCheck;
    [keyP, keyC, firstR, lastP, lastR] = KbQueueCheck(KEY_ST.KEYINDEX);
    if phase <= size(schedule,1),curEV = schedule(phase,:);end
    if keyC(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;break;end;

    % Stim presentation
    if curTime > curEV(~isnan(curEV)) && phase <= size(schedule,1)
        % Change status of stimuli (off ->on, on ->off)
        fIDX = find(~isnan(curEV));
        update = 1 - stateV(fIDX);
        stateV(fIDX) = logical(update);
        phase = phase + 1; % move on to next phase

        % Update screen
        locDrawStim(stateV);
        locDrawFC('NEUTRAL');
        Screen('Flip',W);
        if update,refRT(fIDX)=GetSecs;locCodeEvent('TEMP',fIDX);end
        if update,RECORD_ST.('EV_STIM'){BLOCK,TRIAL}(fIDX)=GetSecs-BLOCKSTS;end
    end
    
    % Always detect responses
    if keyP && stateR(2)<TASK_ST.NSTIM 
        % Check response
        resp=find(keyC);
        if length(resp)~=1,resp=99;end
        secs=keyC(resp);
        [referedKey, ~] = find(KEY_ST.REFKEY==resp);
        
        % When specified key sets are pressed
        if ~isempty(referedKey)
            % Response states 
            stateR(1) = referedKey;
            stateR(2) = stateR(2)+1;
            sIDX = spos == stateR(1);% index of referred SR
            
            sca
            keyboard
            
            % Record responses (first-second item)
            % 2nd column could be faster than first
            RECORD_ST.BLOCK(BLOCK,TRIAL)=BLOCK;
            RECORD_ST.TRIAL(BLOCK,TRIAL)=TRIAL;
            RECORD_ST.RT{BLOCK,TRIAL}(sIDX)=(secs-refRT(sIDX))*1000;
            RECORD_ST.TASKTS{BLOCK,TRIAL}(sIDX)=(secs-BLOCKSTS)*1000;
            RECORD_ST.RESP{BLOCK,TRIAL}(sIDX)=resp;
            RECORD_ST.CORRESP{BLOCK,TRIAL}(sIDX)=coresp(sIDX);
            RECORD_ST.RESPORDER{BLOCK,TRIAL}(stateR(2))=stateR(1); % refered order!
            if coresp(sIDX)==resp,acc=1; else, acc=0;end
            RECORD_ST.ACC{BLOCK,TRIAL}(sIDX)=acc;
            RECORD_ST.('EV_RESP'){BLOCK,TRIAL}(sIDX)=GetSecs-BLOCKSTS;
            KbReleaseWait;% Wait for keyboard release
        end
    end    
end

%PHASE1:FEEDBACK-----------------------------------------------------------
% Debugging
% 4  22
% 29 27
check = RECORD_ST;
check = rmfield(check, {'SUBID','SESSION','DFILE','DLOADED','PARAMORDER'});
disp(structfun(@(x) x(BLOCK,TRIAL),check, 'Uni', false))

% Display feedback
fb = {'FB_BAD', 'FB_GOOD'};
acc = RECORD_ST.ACC{BLOCK,TRIAL};
locDrawFC(fb{all(acc)+1});
Screen('Flip', W);
WaitSecs(TIME_ST.FBT);

% Clean up
KbQueueFlush;
KbEventFlush;

end

%##########################################################################

function locCodeEvent(codeT,code)
%==========================================================================
% locCodeEvent creates three different kind of event codes
% temporal event codes:codes for individual events regardress of contents
% characteristic codes:codes for the contents of one trial(all conditions)
%==========================================================================
%STIMULUS PROPERTY / TIME REGULATION / COLOR
global BLOCK TRIAL RECORD_ST CONFIG TRIGGER_ST W
global EXC_ST TASK_ST PARAM_ST FPATH_ST STIM_ST TIME_ST COLOR_ST KEY_ST
%**************************************************************************

switch codeT
    case {'TEMP'}%Temporal codes
        portCode=FPATH_ST.PORTID_S{:}; % for Oregon ECAP
        eventCode={code};
        
    case {'CHAR'}%Characteristic codes
        portCode=FPATH_ST.PORTID_S{:}; % for Oregon ECAP
                
        %-Block (50~): (block+50)
        eventCode(1,:)={BLOCK+50};
        
        %-Trial (20~): (trial+20)
        eventCode(2,:)={TRIAL+20};
        
    case {'RESP'}%Response codes (only for EEG_EC)
        portCode=FPATH_ST.PORTID_R{:};
        
        %-Accuracy:11 (correct),10 (wrong)
        eventCode(1,:)={RECORD_ST.ACC(BLOCK,TRIAL)+10};
end

%Actually write event code here!
for i=1:size(eventCode,1)
    switch [EXC_ST.EXMODE]
        case {2}, IOPort('Write', CONFIG, uint8(eventCode{i,:}), 0);% BP EEG (RIKEN)
        case {1},write_parallel(portCode,eventCode{i,:});% ECAP EEG    
        otherwise,disp(strcat(portCode,'_',num2str(eventCode{i})));
    end
    
    % IOPor requires a short pause and neutral code (unit8(0))
    pause(TIME_ST.YEILDT); % for 1000Hz wait for 1ms 
    if EXC_ST.EXMODE==2,IOPort('Write', CONFIG, uint8(0), 0);end
end

% sca
% keyboard

end

%##########################################################################

function locDrawStim(stateV)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPATH_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL TRIGGER_ST 
%**************************************************************************

% Initialize
Screen('TextSize',W,STIM_ST.FONTSTIM);
adj = floor(STIM_ST.FONTSTIM/10);
pos = PARAM_ST.STIMPOS{BLOCK,TRIAL};
feat_n = [TASK_ST.STIM(1, PARAM_ST.F_NUMBER{BLOCK, TRIAL})];
feat_c = [TASK_ST.STIM(2, PARAM_ST.F_COLOR{BLOCK, TRIAL})];

% Override with neutral color (after display time)
for i = 1:TASK_ST.NSTIM
    num = num2str(feat_n{i});
    color = COLOR_ST.(feat_c{i});
    if ~stateV(i),color = COLOR_ST.BACKG;end
    posY = STIM_ST.FROFFSET(pos(i));
    CenterTextOnPoint(W,num,STIM_ST.XC-4,STIM_ST.YC+posY+adj,color);
end

end

%##########################################################################

function locDrawFC(workT)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPATH_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL TRIGGER_ST 
%**************************************************************************
% Place holder
for offset = STIM_ST.FROFFSET
    frame=CenterRectOnPoint(STIM_ST.FRRECT,STIM_ST.XC,STIM_ST.YC+offset);
    Screen('FrameRect', W, COLOR_ST.WHITE, frame, STIM_ST.FRSIZE);
end

% Fixation cross
switch workT
    case {'NEUTRAL'},color='WHITE';
    case {'FB_GOOD'},color='GREEN';
    case {'FB_BAD'},color='RED';
end
fixation='+';Screen('TextSize',W,STIM_ST.FONTFC);
DrawFormattedText(W,fixation,'center',STIM_ST.YC+STIM_ST.FXADJ,COLOR_ST.(color));

end


%##########################################################################

function locShowMessege(msgtime)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPATH_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL PRACEXC IMAGE_ST
%**************************************************************************

%Consistant settings
Screen('FillRect',W,COLOR_ST.BACKG);
Screen('TextSize',W,STIM_ST.FONTSIZE);

switch msgtime
    case {'BINTRO'}% Begining Block BREAK screen
        [ConText]=locCondMessage(msgtime);
        TEXT(1)={ConText{1}};
        TEXT(2)={ConText{2}};
        TEXT(3:8)={''};
        TEXT(9)={'Press any key to proceed'};
        TEXTC = repmat({'WHITE'},1,length(TEXT));
        TEXTC(2:3) = {'GREEN'};
        
        % Draw instruction image
        s=CenterRectOnPoint(STIM_ST.INSTRUCTION,STIM_ST.XC+100,STIM_ST.YC+50);
        Screen('DrawTexture',W,IMAGE_ST.IMARRAY{1,end},[],s,0,[],[],[]);
        
    case {'BRESULT'}% Ending Block BREAK screen (no money) 
        [ConText]=locCondMessage(msgtime);
        TEXT(1)={sprintf('You  have completed block %d of %d',BLOCK,EXC_ST.NUMTBLOCK)};
        TEXT(2)={ConText{1}};
        TEXT(3)={ConText{2}};
        TEXT(4:8)={''};
        TEXT(9)={'Press any key to proceed'};
        TEXTC = repmat({'WHITE'},1,length(TEXT));
        
    case {'EBREAK'}% Break time
        TEXT(1)={'If you have questions about the task, please ask the experimenter'};
        TEXT(2)={'If not, please continue whenever you are ready'};
        TEXT(3:8)={''};
        TEXT(9)={'Press any key to proceed'};
        TEXTC = repmat({'WHITE'},1,length(TEXT));
        
    case {'EXEND'}% Ending Experiment screen
        TEXT(1)={'Thank you'};
        TEXT(2)={'You have completed the experiment'};
        TEXT(3:8)={''};
        TEXT(9)={'Press any key to proceed'};
        TEXTC = repmat({'WHITE'},1,length(TEXT));
end

%Drawing messages
for i=1:size(TEXT,2)
    text=TEXT{i};xp=STIM_ST.XC-375;yp=STIM_ST.YC-300+(80*(i-1));
    Screen('DrawText',W,text,xp,yp,COLOR_ST.(TEXTC{i}));
end

% Photosensor trigger for ECoG: flip black mark to white!
Screen('Flip',W);

%Prevents jumping due to key pressed ahead
%and waits for keyboard input
while KbCheck; end;KbWait;WaitSecs(1.5);

%Exit out if esc is pressed
[~,~,keyCode]=KbCheck;
if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;

end

%##########################################################################

function [ConText]=locCondMessage(msgtime)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPATH_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL PRACEXC
%**************************************************************************

switch msgtime
    case {'BINTRO'}% Begining Block BREAK screen
         if ~PRACEXC, cond='a PRACTICE'; else, cond='an EXPERIMENT';end
         ConText(1)={sprintf('ATTENTION!! This is %s!',cond)};
         tasksets = reshape(TASK_ST.TASK',1,[]);%'Low-High','Odd-Even','Red-Blue','Bold-Faint'
         tasksetsC = tasksets(PARAM_ST.TASKSET{BLOCK, 1});
         ConText(2)={sprintf('1:%s 2:%s',tasksetsC{1},tasksetsC{2})};
        
    case {'BRESULT'}% Begining Block BREAK screen  
        %Basic information
        prevRT=RECORD_ST.RT(EXC_ST.NUMPBLOCK+1:max(BLOCK-1,1),:);
        prevRT = [prevRT{:}];% unpack from cell
        useIDX=prevRT>0;prevRT=prevRT(useIDX);
        cutRT=prctile(prevRT,EXC_ST.RTCUT);
        acc=cellfun(@all,RECORD_ST.ACC(BLOCK,1:TRIAL));
        fastRTs=cellfun(@(x) x<cutRT,RECORD_ST.RT(BLOCK,1:TRIAL),'Uni',false);
        fastRTs = [fastRTs{:}];
        
        %Critria for incentives/points
        goodF(1)=(length(find(acc))/TRIAL)*100>EXC_ST.ACCCUT;
        goodF(2)=true;
        goodF(3)=TRIAL>EXC_ST.RESPCUT;
        if all(goodF)&&PRACEXC,moneyE=length(fastRTs);else moneyE=0;end
        
        %Number of response check
        moneyE = ceil(moneyE/2);
        RECORD_ST.MONEY(BLOCK,:)=moneyE;
        moneyT=sum(RECORD_ST.MONEY(1:BLOCK,1));
        
        %If both accuracy and RT are good, then add money!
        ConText(1)={sprintf('You missed %d out of %d trials',length(find(~acc)),length(acc))};
        ConText(2)={sprintf('Please go on to the next block!')};
        if PRACEXC,ConText(2)={sprintf('You earned %d points, and the total is %d points',moneyE,moneyT)};end
        if PRACEXC,disp(sprintf('You earned %d points, and the total is %d points',moneyE,moneyT));end
end

end

%##########################################################################

function locLoadImage
%==========================================================================
% locLoadStim loads up all of visual and auditory data. Right now
% compatible for jpg and wav data. Provide the directory path for the file
% location and exact file names, and it will create a summary cell array.
%==========================================================================
global FPATH_ST IMAGE_ST W 
%**************************************************************************

%Image File Names
%1=filename
%2=imagedata
IMAGE_ST.IMFPATH=FPATH_ST.SPDIR;
IMAGE_ST.IMARRAY=cell(1,3);
IMAGE_ST.IMARRAY(1,1)={'DimensionSwitch_PRP_instruction_J.jpg'};

% Visual stimuli preparation
for i=1:size(IMAGE_ST.IMARRAY,1)
    %reading jpg image data
    imaCode = char([IMAGE_ST.IMFPATH filesep IMAGE_ST.IMARRAY{i,1}]);
    IMAGE_ST.IMARRAY{i,2}=imread(strtrim(imaCode));
    IMAGE_ST.IMARRAY{i,3}=Screen('MakeTexture',W,IMAGE_ST.IMARRAY{i,2});
end

end