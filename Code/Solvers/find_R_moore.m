function [R] = find_R_moore(L,kwargs)
%FIND_R_MOORE Finds the function R on [-L(0), tmax+L(tmax)] using Moore's method and symbolic derivatives
%By default the optimal number of gamma_l's is sought and selected
% INPUTS:
%   L: function_handle or interpolant object, length of the domain
%   optional t_max: double, maximum time to calculate R for
%   optional density: double, density of xi-points in the interpolant object
%   optional Kmax: maximum number of gamma_l's to be calculated (more than 4 is not feasible except for simple examples)
%   optional margin: if specified, the function R will be calculated from -L(0)-margin to tmax+L(tmax)+margin (a little outside the domain specified)
%   optional K_array: if specified, ignore the minimization process of the error and instead
%           return a cell array with an interpolant R for each value Kmax in the K_array
%   optional param_array: it is possible to have a symbolic parameter/parameters in the function L,
%           in that case K_array and param_array should be specified. Param_array is a cell array of size
%           1xN, with N the number of parameters in L. Each element in param_array is a double containing
%           each value at which the given parameter should be evaluated (see examples), e.g. [0.3;0.5;0.7].
%           If only 1 parameter is present then param_array can also be a double containing each value 
%           instead of a cell array with 1 entry. If there is more than 1 parameter, specify them alphabetically.
% OUTPUTS:
%   R: interpolant object, returns the function R
% EXAMPLES:
%   [L,R] = get_LR("linearL")
%   R = find_R_moore(L,Kmax=4) 
%       finds the optimal number of terms (between 1 and 5) and returns an R with the optimal number
%   R = find_R_moore(L,K_array=[1,2,3,4])
%       returns a 1x4 cell array with an R for Kmax = 1, 2, 3 and 4
%   [L,R] = get_LR("linearL",v=sym('v'))
%   R = find_R_moore(L,K_array=[1,2,3,4],param_array=[0.3;0.5;0.7])
%       returns a 3x4 cell array with an R for Kmax=1,2,3,and 4 (axis 2) and for v=0.3,0.5,and 0.7 (axis 1)
%   R = find_R_moore(L,K_array=[1,2,3,4],param_array={[0.3;0.5;0.7]})
%       same as above
%   [L,R] = get_LR("linearL",v=sym('v'),L0=sym('L0'))
%   R = find_R_moore(L,K_array=[1,2,3,4],param_array={[1;2], [0.3;0.5]})
%       returns a 2x4 cell array with an R for Kmax=1,2,3,4 (axis 2) and for L0=1,v=0.3 and L0=2,v=0.5 (axis 1)
%   R = find_R_moore(L,K_array=[1,2,3,4],param_array={[1;2], reshape([0.3,0.5,0.7],[1,1,3])})
%       returns a 2x4x3 cell array with an R for Kmax=1,2,3,4 (axis 2) and for L0=1,2 (axis 1) and for v=0.3,0.5 (axis 3)
arguments
    L
    kwargs.t_max = 5
    kwargs.density = 1000
    kwargs.Kmax = 100
    kwargs.margin = 0
    kwargs.K_array = []
    kwargs.param_array = {}
end

print_indented_message("Find R Moore",true)
timer = tic;

if class(kwargs.param_array) == "double"
    % kwargs.param_array should be a cell with arrays for each parameter to
    % fill in, in alphabetical order
    kwargs.param_array = {kwargs.param_array};
end
for j = 1:numel(kwargs.param_array)
    if isempty(kwargs.param_array{j})
        % If any of the parameters are not specified, set params_array empty
        kwargs.param_array = {};
        break
    end
end

if ~isempty(kwargs.param_array)
    try
        % Expand the dimensions of K_array, and of each element in param_array
        for j = 1:numel(kwargs.param_array)
            kwargs.K_array = kwargs.K_array + zeros(size(kwargs.param_array{j})); % widen the dimensions
        end
        for j = 1:numel(kwargs.param_array)
            kwargs.param_array{j} = kwargs.param_array{j} + zeros(size(kwargs.K_array)); % widen the dimensions
        end
    catch
        error('K_array and the arrays in param_array should have the same dimension or broadcasting should be possible')
    end
    param_array_reordened = cell(size(kwargs.K_array));
    % Order param_array differently: each element now has a set of parameters for one R
    for j = 1:numel(kwargs.K_array)
        param_array_reordened{j} = cellfun(@(x) x(j), kwargs.param_array, 'UniformOutput', false);
    end
end

syms dR(z)
% derivative of R, namely R'
tmax = kwargs.t_max;

tt = linspace(0,tmax,17281); % test array for minimization

if ~isempty(kwargs.K_array)
    Kmax = max(kwargs.K_array,[],"all");
else
    Kmax = kwargs.Kmax; % make variable
    min_xi = -kwargs.margin-L(0-kwargs.margin);
    max_xi = tmax+kwargs.margin+L(tmax+kwargs.margin);
    xi = linspace(min_xi,max_xi,kwargs.density*(max_xi-min_xi));
end
gamma = zeros(Kmax+1,1,'sym');
derivs = zeros(Kmax+1,1,'sym');
gamma(1) = 1/L(z);
dR(z) = 0;

if ~isempty(kwargs.K_array)  % store the R'
    dR_arr = cell(1,Kmax+1);
end

k = 0;
while k <= Kmax
    print_indented_message("Iteration " + string(k)) % print which iteration
    % Calculate next gamma_k
    for i = 1:k
        derivs(k+1-i) = diff(derivs(k+1-i),2);
        gamma(k+1) = gamma(k+1)-(L(z)^(2*i))/ factorial(2*i+1) * derivs(k+1-i);
    end
    derivs(k+1) = gamma(k+1);
    
    dR(z) = dR + gamma(k+1); % make sure dR depends on a parameter z
    dR_handle = fastMatlabFunction(vpa(dR));

    if isempty(kwargs.K_array)
        % minimize if no K_array is given
        R = interpolant(xi,dR_handle(xi)).integral;
        new_error = rmse(R(tt+L(tt))-R(tt-L(tt)),2);
    
        % check for convergence
        if (k > 0 && new_error >= prev_error) || k == Kmax
            if k > 0 && new_error >= prev_error
                R = prev_R;
                best_k = k-1;
            else
                best_k = k;
            end
            print_indented_message("Best k: " + best_k)
            print_indented_message("Finished in " + string(toc(timer)))
            return
        end
        prev_error = new_error;
        prev_R = R;
    else
        % Add R' handle based on the first k terms to dR_arr (cell array)
        dR_arr{k+1} = dR_handle;
    end

    k = k+1;
end

R = cell(size(kwargs.K_array));
if ~isempty(kwargs.param_array)
    syms L_sym(z); 
    L_sym(z) = L(z);
    L_handle = fastMatlabFunction(L_sym); % make a handle with the symbolic parameter(s) included as (an) argument(s)
end
for j = 1:numel(kwargs.K_array)
    print_indented_message("Calculating R for k = " + kwargs.K_array(j))
    % Numerical integral
    dR_handle = dR_arr{kwargs.K_array(j)+1};
    if ~isempty(kwargs.param_array)
        dR_handle = @(z) dR_handle(z,param_array_reordened{j}{:}); % fill in the parameters before integrating
        L = @(z) L_handle(z,param_array_reordened{j}{:}); % fill in the parameters before integrating
    end
    min_xi = -kwargs.margin-L(0-kwargs.margin);
    max_xi = tmax+kwargs.margin+L(tmax+kwargs.margin);
    xi = linspace(min_xi,max_xi,kwargs.density*(max_xi-min_xi));

    R{j} = interpolant(xi,dR_handle(xi)).integral;
end

print_indented_message("Finished in " + string(toc(timer)))
end

function f = fastMatlabFunction(f_sym)
    vars = num2cell(symvar(f_sym));
    vars_str = cellfun(@(x) [char(x) ','], vars, UniformOutput=false);
    vars_str = [vars_str{:}];
    vars_str(end) = []; % no comma at end
    f = str2func("@(" + vars_str + ")" + vectorize(f_sym));
end