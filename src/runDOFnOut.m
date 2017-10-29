function [state,options,optchanged] = runDOFnOut(options,state,flag)
persistent pophist scrhist
optchanged = false;
switch flag
    case 'init'
        pophist(:,:,1) = state.Population;
        scrhist(:,:,1) = state.Score;
        assignin('base','pophistory',pophist);
        assignin('base','scrhistory',scrhist);
    case 'iter'
        % Update the pophist every 1 generations.
        if rem(state.Generation,1) == 0
            ss = size(pophist,3);
            pophist(:,:,ss+1) = state.Population;
            assignin('base','pophistory',pophist);
            ss = size(scrhist,3);
            scrhist(:,:,ss+1) = state.Score;
            assignin('base','scrhistory',scrhist);
        end
        % Update the plot.
    case 'done'
        % Include the final population in the pophist.
        ss = size(pophist,3);
        pophist(:,:,ss+1) = state.Population;
        assignin('base','pophistory',pophist);
        ss = size(scrhist,3);
        scrhist(:,:,ss+1) = state.Score;
        assignin('base','scrhistory',scrhist);
end
