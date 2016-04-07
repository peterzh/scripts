function InteractiveGLM_callbackPlot(which,g,es,step)
global GLM_PARAMS;
switch(which)
    case 'bL'
        GLM_PARAMS(1) = GLM_PARAMS(1)+step;
%         es.String = num2str(GLM_PARAMS(1));
    case 'bR'
        GLM_PARAMS(3) = GLM_PARAMS(3)+step;
%         es.String = num2str(GLM_PARAMS(3));
    case 'sL'
        GLM_PARAMS(2) = GLM_PARAMS(2)+step;
%         es.String = num2str(GLM_PARAMS(2));
    case 'sR'
        GLM_PARAMS(4) = GLM_PARAMS(4)+step;
%         es.String = num2str(GLM_PARAMS(4));
end

evalCL = linspace(0,max(g.data.contrast_cond(:,1)),100);
evalCR = linspace(0,max(g.data.contrast_cond(:,1)),100);
prop=nan(length(evalCL),length(evalCR),3);
titles = {'pred P( left | contrast)','pred P( right | contrast)','pred P( nogo | contrast)'};
for cl = 1:length(evalCL)
    for cr = 1:length(evalCR)
        p = g.calculatePhat(GLM_PARAMS,[evalCL(cl) evalCR(cr)]);
        for i=1:3
            prop(cl,cr,i) = p(i);
        end
    end
end

for i=1:3
    subplot(1,3,i);
    imagesc(evalCR,evalCL,prop(:,:,i),[0 1]);
    set(gca,'YDir','normal');
    xlabel('C Right');
    ylabel('C Left');
    axis square;
    title(titles{i});
end

disp(GLM_PARAMS);
% text(-0.4,-0.1,num2str(step));
end