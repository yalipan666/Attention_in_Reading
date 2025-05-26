function el_Stop(cfg)
el=cfg.el;

%try
    % STEP 7
    % finish up: stop recording eye-movements,
    % close graphics window, close data file and shut down tracker
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    % download data file
    %try
        fprintf('Receiving data file ''%s''\n', cfg.el.edffile);
        %status=Eyelink('ReceiveFile');
        status=Eyelink('ReceiveFile',cfg.el.edffile,cfg.el.eyedir,1); %transfer file to experiment directory
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(cfg.el.edffile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', cfg.el.edffile, cfg.el.eyedir );
        end
%     catch rdf
%         fprintf('Problem receiving data file ''%s''\n', cfg.edffile );
%         rdf;
%     end
    
    cleanup;
    
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
sca;

% Restore keyboard output to Matlab:
ListenChar(0);