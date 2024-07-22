function [f,g] = get_fg(L, initial_condition,kwargs)
%GET_FG Returns default initial conditions
% INPUTS:
%   L: function_handle, boundary length of the domain
%   initial_condition: str, either 'gauss', 'sin', 'wavepacket' are supported
%   optional s: double, std deviation in gaussian in 'gauss' and 'wavepacket'
%   optional x0: double, shift in gaussian in 'gauss' and 'wavepacket'
%   optional T: double, period in 'sin' and 'wavepacket'
% OUTPUTS:
%   f: function_handle, initial position
%   g: function_hanlde, initial velocity
% EXAMPLES:
%   [f,g] = get_fg(L, "wavepacket", s = 0.3) 
    arguments
        L
        initial_condition
        kwargs.s
        kwargs.x0
        kwargs.T
    end
    if ~isempty(initial_condition)
        switch initial_condition
            case "gauss"
                if isfield(kwargs,'s'); s = kwargs.s; else; s = L(0)/16; end
                if isfield(kwargs,'x0'); x0 = kwargs.x0; else; x0 = L(0)/2; end
                f = @(x) 2*exp(-0.5*((x-x0)/s).^2);
                g = @(x) 0; % no initial speed
            case "sin"
                if isfield(kwargs,'T'); T = kwargs.T; else; T = 2; end
                f = @(x) 2*sin(T*2*pi*x/L(0));
                df = derivative(f); dL = derivative(L);
                g = @(x) -x/L(0) * dL(0) .* df(x);
            case "wavepacket"
                if isfield(kwargs,'s'); s = kwargs.s; else; s = L(0)/16; end
                if isfield(kwargs,'x0'); x0 = kwargs.x0; else; x0 = L(0)/2; end
                if isfield(kwargs,'T'); T = kwargs.T; else; T = 20; end
                f = @(x) 2*exp(-0.5*((x-x0)/s).^2).*sin(T*2*pi*x/L(0));
                df = derivative(f); 
                g = @(x) -df(x); % one pulse  
        end
    end
end