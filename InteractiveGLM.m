%Visualise the effect of parameters on GLM fits
g = GLM('2016-03-16_1_Spemann').setModel('C50-subset').fit;
which = 'sL';

global GLM_PARAMS;
GLM_PARAMS = g.parameterFits;

evalCL = linspace(0,max(g.data.contrast_cond(:,1)),100);
evalCR = linspace(0,max(g.data.contrast_cond(:,1)),100);
prop=nan(length(evalCL),length(evalCR),3);

for cl = 1:length(evalCL)
    for cr = 1:length(evalCR)
        p = g.calculatePhat(g.parameterFits,[evalCL(cl) evalCR(cr)]);
        for i=1:3
            prop(cl,cr,i) = p(i);
        end
    end
end

fig=figure;
titles = {'pred P( left | contrast)','pred P( right | contrast)','pred P( nogo | contrast)'};
for i=1:3
    subplot(1,3,i);
    imagesc(evalCR,evalCL,prop(:,:,i),[0 1]);
    set(gca,'YDir','normal');
    xlabel('C Right');
    ylabel('C Left');
    title(titles{i});
    axis square;
end
b1 = uicontrol(fig,'Style','pushbutton','String','+ bL',...
                'Position',[50 60 60 40]);
b1.Callback = @(es,ed) InteractiveGLM_callbackPlot('bL',g,es,+0.5);

b2 = uicontrol(fig,'Style','pushbutton','String','- bL',...
                'Position',[50 20 60 40]);
b2.Callback = @(es,ed) InteractiveGLM_callbackPlot('bL',g,es,-0.5);

b3 = uicontrol(fig,'Style','pushbutton','String','+ sL',...
                'Position',[100 60 60 40]);
b3.Callback = @(es,ed) InteractiveGLM_callbackPlot('sL',g,es,+0.5);

b4 = uicontrol(fig,'Style','pushbutton','String','- sL',...
                'Position',[100 20 60 40]);
b4.Callback = @(es,ed) InteractiveGLM_callbackPlot('sL',g,es,-0.5);

b5 = uicontrol(fig,'Style','pushbutton','String','+ bR',...
                'Position',[250 60 60 40]);
b5.Callback = @(es,ed) InteractiveGLM_callbackPlot('bR',g,es,+0.5);

b6 = uicontrol(fig,'Style','pushbutton','String','- bR',...
                'Position',[250 20 60 40]);
b6.Callback = @(es,ed) InteractiveGLM_callbackPlot('bR',g,es,-0.5);

b7 = uicontrol(fig,'Style','pushbutton','String','+ sR',...
                'Position',[300 60 60 40]);
b7.Callback = @(es,ed) InteractiveGLM_callbackPlot('sR',g,es,+0.5);

b8 = uicontrol(fig,'Style','pushbutton','String','- sR',...
                'Position',[300 20 60 40]);
b8.Callback = @(es,ed) InteractiveGLM_callbackPlot('sR',g,es,-0.5);
            
% b = uicontrol('Parent',f,'Style','slider','Position',[81,130,419,23],...
%     'value',0, 'min',-10, 'max',10);
% c = uicontrol('Parent',f,'Style','slider','Position',[81,100,419,23],...
%     'value',0, 'min',-10, 'max',10);
% d = uicontrol('Parent',f,'Style','slider','Position',[81,70,419,23],...
%     'value',0, 'min',-10, 'max',10);
% e = uicontrol('Parent',f,'Style','slider','Position',[81,40,419,23],...
%     'value',0, 'min',-10, 'max',10);
% b.Callback = @(es,ed) plotterClbk('bL',g,es.Value);
% c.Callback = @(es,ed) plotterClbk('bR',g,es.Value);
% d.Callback = @(es,ed) plotterClbk('sL',g,es.Value);
% e.Callback = @(es,ed) plotterClbk('sR',g,es.Value);
