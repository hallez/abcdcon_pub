function [b] = run_exceptions_ABCDCon(b) 

if strcmp(b.curSubj,'s001')
    fprintf('\nReplacing values for %s.\n',b.curSubj)
    b.beta_matrix_size = [144,144,52];
end

if strcmp(b.curSubj,'s003')
    fprintf('\nReplacing values for %s.\n',b.curSubj)
    b.runs = {'run1' 'run4'};
end

if strcmp(b.curSubj,'s015')
    fprintf('\nReplacing values for %s.\n',b.curSubj)
    b.runs = {'run1' 'run2' 'run4'};
end

end

