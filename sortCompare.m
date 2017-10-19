
% Script to compare sorting 

%% load mine

d = '\\zserver.cortexlab.net\Data\Subjects\Cori\2016-12-18\ephys_AM\sorting';
cluN = readNPY(fullfile(d, 'spike_clusters.npy'));
st = readNPY(fullfile(d, 'spike_templates.npy'));
[cidN, cgN] = readClusterGroupsCSV(fullfile(d, 'cluster_group.tsv'));

%% load peter's 

d = '\\zserver.cortexlab.net\Data\Subjects\Cori\2016-12-18\ephys_AM\sorting\peter';
cluP = readNPY(fullfile(d, 'spike_clusters.npy'));
% stP = readNPY(fullfile(d, 'spike_templates.npy')); % these are the same
[cidP, cgP] = readClusterGroupsCSV(fullfile(d, 'cluster_group.tsv'));

%% find out how spikes from each template were labeled

outputF = 1;

uST = unique(st);

fprintf(outputF, 'For template X, Nick eventually labeled N, Peter eventually labeled P\n')
fprintf(outputF, '[ ]=same label; [x]=different\n')
totals = zeros(numel(uST),3);
for q = 1:numel(uST)
    [v,i] = countUnique(cluN(st==uST(q)));
    cluForN = v(find(i==max(i),1));
    cgForN = cgN(cidN==cluForN);
    if isempty(cgForN); cgForN=3; end %happens when no label was given
    [v,i] = countUnique(cluP(st==uST(q)));
    cluForP = v(find(i==max(i),1));
    cgForP = cgP(cidP==cluForP);
    if isempty(cgForP); cgForP=3; end %happens when no label was given

    if cgForN==cgForP
        fprintf(outputF, ' ');
    else
        fprintf(outputF, 'x');
    end
    fprintf(outputF, ' %d: %d %d\n', uST(q), cgForN, cgForP);
    totals(q,:) = [uST(q), cgForN, cgForP];
end

totalsDisagree = totals(totals(:,2)~=totals(:,3),:); 
[x1,x2,x3]=unique(totalsDisagree(:,2:3),'rows');
tab = tabulate(x3);
x1 = round([x1 tab(:,3)])


%% compare merges

outputF = 1;

uST = unique(st);

fprintf(outputF, 'For template X (groupN groupP), Nick merged with [A B C], Peter merged with [M N O]\n');
fprintf(outputF, '*=same merges; x=different\n')
for q = 1:numel(uST)
    [v,i] = countUnique(cluN(st==uST(q)));
    cluForN = v(find(i==max(i),1));
    mergedN = unique(st(cluN==cluForN));
    cgForN = cgN(cidN==cluForN);
    if isempty(cgForN); cgForN=3; end %happens when no label was given
    
    [v,i] = countUnique(cluP(st==uST(q)));
    cluForP = v(find(i==max(i),1));
    mergedP = unique(st(cluP==cluForP));
    cgForP = cgP(cidP==cluForP);
    if isempty(cgForP); cgForP=3; end %happens when no label was given

    % take just the other templates that were merged
    mergedN = mergedN(mergedN~=uST(q));
    mergedP = mergedP(mergedP~=uST(q));
    
    if numel(mergedN)>0 || numel(mergedP)>0
        if all(ismember(mergedN, mergedP)) && all(ismember(mergedP, mergedN))
            fprintf(outputF, ' ');
        else
            fprintf(outputF, 'x');
        end
        fprintf(outputF, ' %d (%d %d): [%s] [%s]\n', uST(q), cgForN, cgForP,...
            num2str(mergedN(:)'), num2str(mergedP(:)'));
    end
end