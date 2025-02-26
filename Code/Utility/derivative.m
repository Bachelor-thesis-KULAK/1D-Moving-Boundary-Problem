function [dfun] = derivative(fun, dorder, kwargs)
% DERIVATIVE Calculates the n'th derivative of a function handle (or interpolant)
% INPUTS:
%   fun: function_handle or 'interpolant' object (See interpolant.m) 
%   optional dorder: positive integer, order of derivative (default 1)
%   optional replacenan: boolean, replaces NaN and Inf in output function (default false)
% OUTPUTS:
%   dfun: function_handle 
% EXAMPLE:
%   dfun = derivative(f) calculates the first derivative of f
%   dfun = derivative(f,dorder=4) calculates the fourth derivative of f
    arguments
        fun
        dorder (1,1) double {mustBeInteger} = 1
        kwargs.replacenan logical = false
    end
    if class(fun) == "function_handle" % if L(t) is a function
        syms f(z) df(z)
        f(z) = fun(z);
        df(z) = diff(f, dorder);
        if ~isempty(findSymType(f,'floor')) 
            % define derivative of floor = 0
            df(z) = mapSymType(df(z),'D(floor)',@(Dfloor) 0);
            replaceDiffFloor = @(dif) ( string(feval(symengine,'op',feval(symengine,'op',dif,1),0)) ~= "floor" ) * dif;
            df(z) = mapSymType(df(z),'diff',replaceDiffFloor);
        end
        dfun = matlabFunction(df);
        if kwargs.replacenan
            dfun = @(z) replaceNaN(dfun,z);
        end
    elseif class(fun) == "interpolant"
        dfun = fun.derivative(dorder);
    else
        error(['Function type is not supported.' ...
            ' Specify a function_handle or interpolant (see interpolant.m)'])
    end
end

function y = replaceNaN(f,x)
    y = f(x);
    problempoints = isnan(y) | isinf(y);
    if any(problempoints)
        warning('Inf or NaN in the derivative detected (see derivative.m)')
        deltax = max( eps(x(problempoints)), 10^(-100) );
        yproblempoints = f(x(problempoints)+deltax)/2 + f(x(problempoints)-deltax)/2;
        y(problempoints) = yproblempoints;
        if any( isnan(yproblempoints) | isinf(problempoints) )
            error('Replacing NaN failed.')
        end
    end
end