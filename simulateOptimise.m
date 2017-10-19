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
        p.V1i = v(1);
        p.V1c = v(2);
        p.S1i = v(3);
        p.S1c = v(4);
        p.M2i = v(5);
        p.M2c = v(6);
        p.baseline = v(7);
        p.inact  = v(8);
        
        varargout{1} = p;
        
    case 'pstruct2vec'
        p = varargin{1};
        varargout{1} = [p.V1i p.V1c p.S1i p.S1c p.M2i p.M2c p.baseline p.inact];
         
end
        
end