function [cfg] = el_Set_Params(cfg)

%el_Set_Params 
%Custom parameter file for eyelink

el=cfg.el.defaults;

el.eye_used                = cfg.el.Eyeused;
el.calibrationtargetsize   = 1;
el.calibrationtargetwidth  = 0.5;
el.targetbeep              = 0;
el.feedbackbeep            = 0;
el.displayCalResults       = 1;
el.eyeimagesize            = 50;  % percentage of screen

%%% color settings
el.backgroundcolour = cfg.ScrBgc;
el.foregroundcolour = cfg.TextColor;
el.msgfontcolour    = cfg.TextColor;
el.imgtitlecolour   = cfg.TextColor;
el.calibrationtargetcolour = cfg.TextColor;

disp('Updating Parameters')
EyelinkUpdateDefaults(el);

cfg.el.defaults=el;
% make sure we're still connected.
if Eyelink('IsConnected')~=1
    warning('eyelink is not connected! restart the tracker');
    cleanup;
    return;
end

% make sure that we get gaze data from the Eyelink
Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT'); 

% This Command is crucial to map the gaze positions from the tracker to
% screen pixel positions to determine fixation
Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld',  cfg.el_rect(1)-cfg.el_rect(1), cfg.el_rect(2)-cfg.el_rect(2), cfg.el_rect(3)-cfg.el_rect(1), cfg.el_rect(4)-cfg.el_rect(2));
Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld',         cfg.el_rect(1)-cfg.el_rect(1), cfg.el_rect(2)-cfg.el_rect(2), cfg.el_rect(3)-cfg.el_rect(1), cfg.el_rect(4)-cfg.el_rect(2));

% Use Psychophysical setting
Eyelink('Command', 'recording_parse_type = GAZE');
Eyelink('Command', 'saccade_motion_threshold = 0.1');
Eyelink('Command', 'saccade_pursuit_fixup = 60');
Eyelink('Command', 'fixation_update_interval = 0');
if strcmp(cfg.el.Eyeused,'BOTH_EYE')
   Eyelink('Command', 'binocular_enabled = YES'); % binocular
else
   Eyelink('Command', 'binocular_enabled = NO'); %monocular
end
%%% YPan: changed accoording to the user mannual to be better for reading study
Eyelink('Command', 'saccade_velocity_threshold = 30');
Eyelink('Command', 'saccade_acceleration_threshold = 8000');

% Other tracker configurations

%% these might crash:
Eyelink('Command', 'heuristic_filter = 0');
Eyelink('Command', 'pupil_size_diameter = YES');

%%
%use 9 point calibration (Default)
% Eyelink('Command', 'calibration_type = HV5');
%Eyelink('Command', 'calibration_type = HV3');
Eyelink('Command', 'calibration_type = HV9');
Eyelink('Command', 'generate_default_targets = YES');
Eyelink('Command', 'enable_automatic_calibration = YES');
Eyelink('Command', 'automatic_calibration_pacing = 1000');
Eyelink('Command', 'use_ellipse_fitter = NO');
Eyelink('Command', 'sample_rate = 1000');
%Eyelink('Command', 'elcl_tt_power = %d', 2); % illumination, 1 = 100%, 2 = 75%, 3 = 50%

%%% enable eye drift checking
Eyelink('Command', 'driftcorrect_cr_disable = NO');

switch cfg.el.Eyeused
    case 'RIGHT_EYE'
        Eyelink('Command', 'file_event_filter = RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,INPUT');
        Eyelink('Command', 'link_event_filter = RIGHT,FIXATION,FIXUPDATE,SACCADE,BLINK,MESSAGE,INPUT');
    case  'LEFT_EYE'
        Eyelink('Command', 'file_event_filter = LEFT,FIXATION,SACCADE,BLINK,MESSAGE,INPUT');
        Eyelink('Command', 'link_event_filter = LEFT,FIXATION,FIXUPDATE,SACCADE,BLINK,MESSAGE,INPUT');
end

Eyelink('Command', 'file_sample_data  = GAZE,GAZERES,HREF,PUPIL,AREA,STATUS,INPUT');

% Cleanup routine:
function cleanup
% Shutdown Eyelink:
    Eyelink('Shutdown');
    el.online = 0;
end
end

