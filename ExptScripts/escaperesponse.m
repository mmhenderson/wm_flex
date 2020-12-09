%DEFINE ESCAPE RESPONSE FUNCTION------------------------------------------
function escaperesponse(OriginalCLUT)
    %fix the color table, otherwise it keeps the corrected one
    if ~isempty(OriginalCLUT) && ~isempty(OriginalCLUT)
        if exist('ScreenNr','var')
            Screen('LoadCLUT', ScreenNr, OriginalCLUT);
        else
            Screen('LoadCLUT', 0, OriginalCLUT);
        end
    end
    Screen('CloseAll');                
    ShowCursor;
    if IsWin
        ShowHideWinTaskbarMex;     
    end
    ListenChar(1)
    error('User exited program.');
end