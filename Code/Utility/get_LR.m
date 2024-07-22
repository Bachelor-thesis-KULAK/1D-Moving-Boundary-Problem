function [L,R] = get_LR(boundary_condition,kwargs)
%GET_LR Returns default boundary conditions
% INPUTS:
%   boundary_condition: str, 
%       linearL : L(t) = L0+v*t
%       sinhR   : R(xi) = A*sinh(k*(xi-xi0))
%       expR    : R(xi) = exp(k*(xi-xi0))
%       expL    : L(t) = L0*exp(-k*t)
%       sinL    : L(t) = L0 + A*sin(k*t)
%   optional L0: double or sym
%   optional v: double or sym
%   optional k: double or sym
%   optional xi0: double or sym
%   optional A: double or sym
% OUTPUTS:
%   L: function_handle
%   R: function_handle
% EXAMPLES:
%   [L,R] = get_LR("linearL", v = 0.5)
%   [L,R] = get_LR("sinhR", k = 0.5, A = 2, xi0 = 0)
    arguments
        boundary_condition
        kwargs.L0  (1,1) {mustBeA(kwargs.L0,  {'double', 'sym'})}
        kwargs.v   (1,1) {mustBeA(kwargs.v,   {'double', 'sym'})}
        kwargs.k   (1,1) {mustBeA(kwargs.k,   {'double', 'sym'})}
        kwargs.xi0 (1,1) {mustBeA(kwargs.xi0, {'double', 'sym'})}
        kwargs.A   (1,1) {mustBeA(kwargs.A,   {'double', 'sym'})}
    end

    R = [];
    if ~isempty(boundary_condition)
        switch boundary_condition
            case "linearL"
                if isfield(kwargs,'L0'); L0 = kwargs.L0; else; L0 = .5; end
                if isfield(kwargs,'v'); v = kwargs.v; else; v = .3; end
                L = @(t) L0+v*t;
                if v ~= 0; R = @(xi) 2*log(1+v/L0*xi)/(log((1+v)/(1-v))); else; R = @(xi) xi/L0; end
            case "sinhR"
                if isfield(kwargs,'k'); k = kwargs.k; else; k = 1; end
                if isfield(kwargs,'xi0'); xi0 = kwargs.xi0; else; xi0 = 1; end
                if isfield(kwargs,'A'); A = kwargs.A; else; A = 1; end
                L = @(t) 1/k * asinh(1./(A*cosh(k*(t-xi0))));
                R = @(xi) A*sinh(k*(xi-xi0));
            case "expR"
                if isfield(kwargs,'k'); k = kwargs.k; else; k = .8; end
                if isfield(kwargs,'xi0'); xi0 = kwargs.xi0; else; xi0 = -1; end
                L = @(t) 1/k * asinh(1./exp(k*(t-xi0)));
                R = @(xi) exp(k*(xi-xi0));
            case "expL"
                if isfield(kwargs,'L0'); L0 = kwargs.L0; else; L0 = 1; end
                if isfield(kwargs,'k'); k=kwargs.k; else; k=0.5; end
                L = @(t) L0*exp(-k*t);
            case "sinL"
                if isfield(kwargs,'L0'); L0 = kwargs.L0; else; L0 = 1; end
                if isfield(kwargs,'A'); A = kwargs.A; else; A = 0.01; end
                if isfield(kwargs,'k'); k = kwargs.k; else; k = 50; end
                L = @(t) L0 + A*sin(k*t);
        end
    end
end