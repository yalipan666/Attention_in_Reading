function cfg = el_Start_SameWindow(cfg)
% Used in FG experiment
% Open screen for calibration, calibrate and start recording


%try
    % STEP 1
    % Open a graphics window on the main screen
    % using the PsychToolbox's Screen function.    
    % use the shrunk version of the window
    %window=Screen('OpenWindow', cfg.screenNumber, [] ,cfg.el_rect);
%     cfg.window=Screen('OpenWindow', cfg.screenNumber);
    
    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    % Psychtoolbox defaults function
    cfg.el.defaults = EyelinkInitDefaults(cfg.window);
    
    % Disable key output to Matlab window:
    ListenChar(2);
    
    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    else
        disp('Eyelink initizalized')
    end
    
    % open file to record data to
    disp('Opening EDF file');
    status = Eyelink('Openfile', cfg.el.edffile);
    
    if ~status
        disp('EDF file opened on Eyelink computer')
    else
        error(['Could not open EDF file on Eyelink computer, error: ' int2str(status)])
    end
    
    % set custom parameters
    disp('Setting parameters')
    cfg = el_Set_Params(cfg);
  
        % Calibrate the eye tracker
        disp('Starting calibration')
        EyelinkDoTrackerSetup(cfg.el.defaults);
        % do a final check of calibration using driftcorrection
%         EyelinkDoDriftCorrection(cfg.el.defaults);

    
    % STEP 5
    % start recording eye position
    disp('Start recording')
    %sca
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.1);
    % mark zero-plot time in data file
    disp('Sending message')
    Eyelink('Message', 'SYNCTIME');
    
    %sca
    ListenChar(0);    
% catch
%     %this "catch" section executes in case of an error in the "try" section
%     %above.  Importantly, it closes the onscreen window if its open.
%     cleanup;
%     psychrethrow(psychlasterror);
% end %try..catch.


% Cleanup routine:
function cleanup
% Shutdown Eyelink:
Eyelink('Shutdown');

% Close window:
%sca;

% Restore keyboard output to Matlab:
ListenChar(0);

