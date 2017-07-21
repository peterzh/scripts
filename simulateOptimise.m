function varargout = simulateOptimise(what,varargin)

switch(what)
    case 'loglik'
        d = varargin{1};
        v = varargin{2};
        p = simulateOptimise('pvec2struct',v);
        
        lik = nan(size(d.response));
        for t = 1:length(d.response)
            p_hat = simulateMechanisticFcn(d.stimulus(t,:),p,d.areaIdx(t),d.bilateral(t));
            lik(t) = p_hat(d.response(t));
        end
        
        varargout{1} = sum(log2(lik));
                
    case 'pvec2struct'
        v = varargin{1};
        
        p = struct;
        p.weights = reshape(v(1:end-2),2,6);
        p.baselineActivity = v(end-1);
        p.InactivationScalingFactor = v(end);
        varargout{1} = p;
        
    case 'pstruct2vec'
        p = varargin{1};
        
        varargout{1} = [p.weights(:); p.baselineActivity; p.InactivationScalingFactor];
end
        
end