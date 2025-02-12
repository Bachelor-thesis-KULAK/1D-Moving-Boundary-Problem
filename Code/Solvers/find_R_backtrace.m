function R = find_R_backtrace(L,kwargs)
%FIND_R_BACKTRACE Finds the transformation R by backtracing
%Based on R(t+L(t)) = 2+R(t-L(t))
% INPUTS:
%   L: function_handle or interpolant object, length of the domain
%   optional initialR: 'quadratic' or 'cubic' (= default) (or 'linear' or 'sigmoid-like'), the choice of R on [-L0,L0]
%   optional derivL0: specify the first (and second) derivative of L at t = 0 instead of calculating it
% OUTPUTS:
%   R: function handle, returns the transformation R
% EXAMPLE:
%   [L,R] = get_LR("linearL")
%   R = find_R_interpolation(L)

arguments
    L
    kwargs.initialR = "cubic"
    kwargs.derivL0 = []
end

print_indented_message("Find R Backtrace",true)
timer = tic;

L0 = L(0);
if isempty(kwargs.derivL0)
    Ldot = derivative(L); Ldot0 = Ldot(0);
    if kwargs.initialR == "cubic"
        Lddot = derivative(Ldot); Lddot0 = Lddot(0);
    end
else
    Ldot0 = kwargs.derivL0(1);
    if numel(kwargs.derivL0) > 1
        Lddot0 = kwargs.derivL0(2);
    end
end

if kwargs.initialR == "linear"
    R_0_fun = @(xi_0) (xi_0+L0)/(L0);
elseif kwargs.initialR == "quadratic"
    deriv_left = (1+Ldot0)/L0; deriv_right = (1-Ldot0)/L0; % make sure derivative is continuous
    R_0_fun = @(xi_0) (xi_0 + L0) .* (deriv_left + 1/(4*L0) * (xi_0 + L0) * (deriv_right-deriv_left) ); % first block, used for interpolation for next blocks
elseif kwargs.initialR == "cubic"
    b = 3/L0 * (1+Ldot0+Ldot0^2+Ldot0^3) / (3+Ldot0^2+L0*Lddot0);
    c = -3/(2*L0^2) * (Ldot0+2*Ldot0^2+Ldot0^3-L0*Lddot0) / (3+Ldot0^2+L0*Lddot0);
    d = 1/(2*L0^3) * (2*Ldot0^2-L0*Lddot0) / (3+Ldot0^2+L0*Lddot0);
    R_0_fun = @(xi_0) b * (xi_0 + L0) + c * (xi_0 + L0).^2 + d * (xi_0 + L0).^3;
elseif kwargs.initialR == "sigmoid-like"
    R_0_fun = @(xi_0) 2/pi * atan(sinh(5*tan(pi/2*xi_0/L0))) + 1;
end
R = @(xi) backtrace(xi,L,R_0_fun);

print_indented_message("Finished in " + string(toc(timer)))
end

function y = backtrace(x,L,R_0_fun)
    y = x;
    plus2 = zeros(size(x));
    L0 = L(0);
    for i = 1:numel(y)
        j = 0;
        while y(i) > L0
            t = fzero(@(t) t+L(t)-y(i), [0,y(i)]); % look between t=0 and t=y(i)
            y(i) = t-L(t);
            j = j+1;
        end
        plus2(i) = j;
    end
    y = R_0_fun(y) + 2*plus2;
end