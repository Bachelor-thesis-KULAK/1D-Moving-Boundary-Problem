classdef solutionObject
    % Object representing a solution of the wave equation (either numerical or analytical). Contains u, x, t.
    % CONSTRUCTOR:
    %   u: MxN double, containing the wave amplitudes corresponding to x and t   
    %   x: MxN (or 1xN) double, space vector
    %   t: 1xN double, time vector
    % EXAMPLE:
    %   solution = solutionObject(u,x,t)   
    % PROPERTIES:
    %   u: MxN double, containing the wave amplitudes corresponding to x and t   
    %   x: MxN (or 1xN) double, space vector
    %   t: 1xN double, time vector
    % METHODS:
    %   checkdims: checks if the dimensions of u, x and t agree
    %   reduceToTimes: makes a new solutionObject where the time is reduced to tnew. Interpolation is applied if needed

    properties
        u
        x
        t
    end

    methods
        % Constructor
        function obj = solutionObject(u,x,t)
            arguments
                u double
                x double
                t double
            end
            obj.u = u;
            obj.x = x;
            obj.t = t;
            % check if the dimensions of x, t and u agree
            obj.checkdims() 
        end

        function checkdims(obj)
            %CHECKDIMS Checks if the dimensions of x, t and u agree
            nx = size(obj.u,1);
            nt = size(obj.u,2);
            if size(obj.x,1) ~= nx || ( size(obj.x,2) ~= nt && size(obj.x,2) ~= 1 )
                size_x = num2str(size(obj.x));
                size_u = num2str(size(obj.u));
                error("Dimensions of x (["+size_x+"]) and u (["+size_u+"]) do not agree")
            elseif size(obj.t,2) ~= nt
                size_t = num2str(size(obj.t));
                size_u = num2str(size(obj.u));
                error("Dimensions of x (["+size_t+"]) and u (["+size_u+"]) do not agree")
            end
        end

        function obj_new = reduceToTimes(obj,tnew)
            %REDUCETOTIMES Makes a new solutionObject where the time is reduced to tnew. Interpolation is applied if needed
            if size(obj.t,1) ~= 1
                error('In reduceToTimes: tnew must be a 1xN array')
            end
            if size(obj.t,2) == size(tnew,2) && all(abs(obj.t-tnew) < 1e-12)
                obj_new = obj;
                return
            end
            timer = tic;

            [found, indices] = ismembertol(tnew,obj.t,1e-12);
            if all(found,'all') % If no interpolation needed, just reduce the existing arrays
                if size(obj.x,2) ~= 1
                    xnew = obj.x(:,indices);
                end
                unew = obj.u(:,indices);
            else
                [~,~,indices] = histcounts(tnew,obj.t); 
                if any(indices == 0)
                    error('The t array that was given has an element not in [0, tmax]')
                end
                % Warn if interpolation is going to be executed
                warning('Interpolation is being used, consider using a different discretization to avoid interpolation')
                left_neighbor = obj.t(indices); right_neighbor = obj.t(indices+1);
                modulo = (tnew-left_neighbor)./(right_neighbor-left_neighbor);
                % modulo can be 1 at the end of the time array (otherwise < 1)
                
                % use linear interpolation between two values
                if size(obj.x,2) ~= 1
                    xnew = (1-modulo) .* obj.x(:,indices) + modulo .* obj.x(:,indices+1);
                end
                unew = (1-modulo) .* obj.u(:,indices) + modulo .* obj.u(:,indices+1);
            end
            obj_new = solutionObject(unew,xnew,tnew); % save as new solutionObject
            print_indented_message("Reducing to desired times in "+string(toc(timer)))
        end
    end
end