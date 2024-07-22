function print_indented_message(message_string, hyperlink)
%PRINT_INDENTED_MESSAGE Print a message string to the console, the string is indented appropiately based on the function call stack
% Depending of the value silent_flag() returns, the message is or is not logged
%
% INPUTS:
% message_string: string or char to print
% optional bold: logical, if true the message will be bold and have a hyperlink to it (default is false)

arguments
    message_string = ""
    hyperlink = false
end

silent = silent_flag();
if silent == 2 || silent == 1 && hyperlink == false
    return % Do not print
end

% Get function call stack
stack = dbstack();
% Check if top of call stack is a script
topIsScript = 0;
try
    nargin(stack(end).name);
catch 
    % error: file is script
    topIsScript = 1;
end
stack_depth = length(dbstack())-1-topIsScript-hyperlink; 
% indent depth = stack depth - 1 (from current function: print_indented_message) -
% topIsScript (if script, then -1) - hyperlink (hyperlink is one less indent)

fprintf(repmat('\t',1,stack_depth))
if hyperlink && numel(stack) > 1
    message_string = sprintf('<a href="matlab:open %s" style="font-weight:bold" rel="nofollow">%s</a>', stack(2).name, message_string);
    % make message bold and hyperlink
    % https://undocumentedmatlab.com/articles/bold-color-text-in-the-command-window/#comment-103035
end
fprintf(message_string + "\n")

end