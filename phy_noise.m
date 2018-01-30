function phy_noise(ksDir)
    s = loadKSdir(ksDir);
    templateIDs = unique(s.spikeTemplates);

    f = figure('color','w');
    tempAx = subplot(8,1,1:7);
    marginalAx = subplot(8,1,8); xlim(marginalAx, [1 size(s.temps,3)]);
    
    for clu = 1:length(templateIDs)
        cluID = templateIDs(clu);
        
        template = double(squeeze(s.temps(cluID + 1,:,:))); %+1 because cluID is 0-indexed 
        marginal = mean(abs(template),1);
        marginal = marginal/max(marginal);
        
        maxMarginal = max(abs(template),[],1);
        maxMarginal = maxMarginal/max(maxMarginal);
        
%         imagesc(tempAx,template); 
%         plot(marginalAx,maxMarginal);
        
        probabilityNoise(clu,1) = trapz(marginal)/length(marginal); %Proportion of AUC      
        tempScalingAmp(clu,1) = median(s.tempScalingAmps(s.spikeTemplates==cluID));
    end

    writePhyTSV(ksDir, 'pNoise', templateIDs, probabilityNoise);
    writePhyTSV(ksDir, 'Amps', templateIDs, tempScalingAmp);

end