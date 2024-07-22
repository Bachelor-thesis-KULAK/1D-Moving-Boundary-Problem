function [dfun] = derivative(fun, dorder, replacenan)
% DERIVATIVE Calculates the n'th derivative of a function handle (or interpolant)
% INPUTS:
%   fun: function_handle or 'interpolant' object (See interpolant.m) 
%   dorder: positive integer, order of derivative
%   replacenan: boolean, replaces NaN and Inf with 0 in output function
% OUTPUTS:
%   dfun: function_handle 
% EXAMPLE:
%   dfun = derivative(f) calculates the first derivative of f
%   dfun = derivative(f,4) calculates the fourth derivative of f
    arguments
        fun
        dorder (1,1) double {mustBeInteger} = 1
        replacenan logical = false
    end
    if class(fun) == "function_handle" % if L(t) is a function
        syms f(z) df(z)
        f(z) = fun(z);
        df(z) = diff(f, dorder);
        dfun = matlabFunction(df);
        if replacenan
            dfun = @(z) replaceNaN(dfun(z));
        end
    elseif class(fun) == "interpolant"
        dfun = fun.derivative(dorder);
    else
        error(['Function type is not supported.' ...
            ' Specify a function_handle or interpolant (See interpolant.m)'])
    end
end

function O = replaceNaN(I)
    problempoints = isnan(I) | isinf(I);
    O = I;
    if any(problempoints)
        warning('Inf or NaN in the derivative detected (See derivative.m)')
        O(problempoints) = 0;
    end
end