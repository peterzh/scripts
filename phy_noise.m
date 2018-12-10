function phy_noise(ksDir)
    s = loadKSdir(ksDir,struct('excludeNoise',false));
    
    %Immediate set all clusters where numSpikes<=100 to noise
    tab=tabulate(s.clu);
    nIdx = tab(:,2)<=100;
    writePhyTSV(ksDir, 'group', tab(nIdx), repmat({'noise'},sum(nIdx),1));
    
%     %Get ACG for each cluster
%     cluIDs = unique(s.clu);
%     
%     SK = nan(length(cluIDs),1);
%     ACGRISE = nan(length(cluIDs),1);
%     for clu = 1:length(cluIDs)
%         spikeTimes = s.st(s.clu==cluIDs(clu));
%         amps = s.tempScalingAmps(s.clu==cluIDs(clu));
%         [xLin, nLin,xLog, nLog] = myACG(spikeTimes,[],[]);
%         nLin=nLin/max(nLin);
%         try
%         ACGRISE(clu)=xLin(find(nLin>0.2,1,'first'));
%         catch
%         end
%         SK(clu)=skewness(amps);
%     end
%     writePhyTSV(ksDir, 'amp_skew', cluIDs, SK );
%     writePhyTSV(ksDir, 'acg_rise', cluIDs, ACGRISE );
%     
    
    templateIDs = unique(s.spikeTemplates);
    for clu = 1:length(templateIDs)
        cluID = templateIDs(clu);
        
        template = double(squeeze(s.temps(cluID + 1,:,:))); %+1 because cluID is 0-indexed 
        marginal = mean(abs(template),1);
        marginal = marginal/max(marginal);
        
        maxMarginal = max(abs(template),[],1);
        maxMarginal = maxMarginal/max(maxMarginal);

        probabilityNoise(clu,1) = trapz(marginal)/length(marginal); %Proportion of AUC      
%         tempScalingAmp(clu,1) = median(s.tempScalingAmps(s.spikeTemplates==cluID));
    end

    writePhyTSV(ksDir, 'pNoise', templateIDs, probabilityNoise);
%     writePhyTSV(ksDir, 'Amps', templateIDs, tempScalingAmp);

end