function [f,g] = get_fg(L,initial_condition,kwargs)
%GET_FG Returns default initial conditions
% INPUTS:
%   L: function_handle, boundary length of the domain
%   initial_condition: str, either 'gauss', 'sin', 'wavepacket' are supported
%   optional s: double, std deviation in gaussian in 'gauss' and 'wavepacket'
%   optional x0: double, shift in gaussian in 'gauss' and 'wavepacket'
%   optional n: double, standing wave mode in 'sin' and 'wavepacket'
%   optional A: double, amplitude in 'sin' and 'wavepacket'
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
        kwargs.n
        kwargs.A
    end
    if ~isempty(initial_condition)
        switch initial_condition
            case "gauss"
                L0 = L(0);
                if isfield(kwargs,'s'); s = kwargs.s; else; s = L0/16; end
                if isfield(kwargs,'x0'); x0 = kwargs.x0; else; x0 = L0/2; end
                if isfield(kwargs,'A'); A = kwargs.A; else; A = 1; end
                f = @(x) A*exp(-0.5*((x-x0)/s).^2);
                g = @(x) 0; % no initial speed
            case "sin"
                L0 = L(0);
                dL = derivative(L); dL0 = dL(0);
                if isfield(kwargs,'n'); n = kwargs.n; else; n = 4; end
                if isfield(kwargs,'A'); A = kwargs.A; else; A = 1; end
                f = @(x) A*sin(n*pi*x/L0);
                g = @(x) -A*n*pi*dL0/L0^2 * x .* cos(n*pi*x/L0);
            case "wavepacket"
                L0 = L(0);
                if isfield(kwargs,'s'); s = kwargs.s; else; s = L0/16; end
                if isfield(kwargs,'x0'); x0 = kwargs.x0; else; x0 = L0/2; end
                if isfield(kwargs,'n'); n = kwargs.n; else; n = 40; end
                if isfield(kwargs,'A'); A = kwargs.A; else; A = 1; end
                f = @(x) A*exp(-0.5*((x-x0)/s).^2) .* sin(n*pi*x/L0);
                g = @(x) - A*(-(x-x0)/s^2 .* sin(n*pi*x/L0) + n*pi/L0*cos(n*pi*x/L0)) .* exp(-0.5*((x-x0)/s).^2);
                % = -df/dx, one pulse
            case "shiftedsin"
                L0 = L(0);
                if isfield(kwargs,'n'); n = kwargs.n; else; n = 4; end
                if isfield(kwargs,'A'); A = kwargs.A; else; A = pi*n/L0; end
                f = @(x) 0;
                g = @(x) A*sin(n*pi*x/L0);
        end
    end
end