function flag = silent_flag(set_flag)
%SILENT_FLAG Set or get the silent flag (a persistent variable). The silent flag is used in print_indented_message to log messages to the console.
% The silent flag can be:
%   0 (default): print_indented_message prints all messages that it receives to the console
%   1: print_indented_message only prints the bold hyperlinked function names to the console
%   2: print_indented_message is completely silent and prints nothing to the console
% INPUTS:
%   optional set_flag: double (0, 1 or 2), if specified set the silent flag to set_flag
% OUTPUTS:
%   flag: get the current silent flag (if set_flag is not specified) or the previous silent flag (if set_flag is specified)
% EXAMPLES:
%   silent_flag(2) sets the silent flag to 2
%   flag = silent_flag() gets the silent flag

persistent silent_flag;
if isempty(silent_flag)
    silent_flag = 0;
end

flag = silent_flag;

if nargin > 0
    if isnumeric(set_flag) && isscalar(set_flag) && any(set_flag == [0,1,2])
        silent_flag = set_flag;
    else
        error('set_flag argument should be either 0, 1 or 2')
    end
end

end