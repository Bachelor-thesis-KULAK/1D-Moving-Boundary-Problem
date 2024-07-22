classdef interpolant
    % Result of an interpolation. Can be evaluated, differentiated and integrated
    % CONSTRUCTOR:
    %   t: array of values OR pp object (if nargin = 1)
    %   y: array of function values,
    %   optional method: 'linear', 'cubic' or 'spline' (= default), the method of interpolation 
    % EXAMPLE:
    %   I = interpolant(0:0.01:1, sin(0:0.01:1)) 
    %   pp = interp1(0:0.01:1, sin(0:0.01:1), 'linear', 'pp')
    %   I = interpolant(pp)
    % PROPERTIES:
    %   PP: piecewise polynomial object (pp)
    % METHODS:
    %   subsref: evaluates in a double of values. Extrapolation is not allowed
    %   derivative: differentiates interpolation objects or integration (by specifying a negative order) 
    %   integral: integrates interpolant while specifying a value at the left endpoint   
    properties
        PP  % Store the piecewise polynomial object
    end
    
    methods
        % Constructor: initialize using the same arguments as interp1
        function obj = interpolant(t,y,method)
            if nargin == 1 && isa(t, 'struct') % Check if the first argument is a pp object
                obj.PP = t; % Initialize using a pp object
            elseif nargin < 2
                error('Not enough input arguments, or wrong PP.')
            elseif nargin > 3
                error('Too many input arguments.')
            else
                if nargin == 2
                    method = 'spline';
                    % methods: linear (C^0), cubic (C^1) and spline (C^2)
                end
                obj.PP = interp1(t,y,method,'pp');
            end
        end
        
        % Subscript method for interpolation
        function values = subsref(obj, S)
            %SUBSREF Evaluates the interpolant in a double of values. Extrapolation is not allowed
            % Check for extrapolation
            if isscalar(S) && S.type == "()" 
                % Important: subsref is of type "()" and must contain only one element. 
                % Otherwise obj.derivative() will trigger this case
                t_eval = S.subs{:};
                pol_left = obj.PP.breaks(1); eval_left = min(t_eval,[],'all');
                pol_right = obj.PP.breaks(end); eval_right = max(t_eval,[],'all');
                out_left = pol_left - eval_left > 10^(-12);
                out_right = eval_right - pol_right > 10^(-12);
                if out_left || out_right
                    if out_right
                        helpful_string = sprintf('Tried to evaluate in %d but right boundary of interpolant is %d', eval_right, pol_right);
                    else
                        helpful_string = sprintf('Tried to evaluate in %d but left boundary of interpolant is %d', eval_left, pol_left);
                    end
                    error('Interpolant is being evaluated out of bounds (extrapolation error):\n%s', helpful_string)
                end
                % Evaluate the values
                values = ppval(obj.PP, S.subs{:});
            else
                % Call builtin behaviour in other cases
                values = builtin('subsref', obj, S);
            end
        end
        
        % Derivative method
        function der = derivative(obj, order)
            %DERIVATIVE Differentiates interpolation objects or integration (by specifying a negative order) 
            if nargin == 1
                order = 1;
            end
            pp_der = fnder(obj.PP, order);
            der = interpolant(pp_der); % Create a CustomInterpolant from the derivative
        end

        % Integral method
        function int = integral(obj, yvalue)
            %INTEGRAL Integrates interpolant while specifying a yvalue at the left endpoint 
            if nargin == 1
                yvalue = 0; % default
            end
            pp_int = fnint(obj.PP, yvalue);
            int = interpolant(pp_int);
        end
    end
end
