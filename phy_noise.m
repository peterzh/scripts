function phy_noise(ksDir)
    s = loadKSdir(ksDir);
    templateIDs = unique(s.spikeTemplates);

    for clu = 1:length(templateIDs)
        cluID = templateIDs(clu);
        
        template = double(squeeze(s.temps(cluID + 1,:,:))); %+1 because cluID is 0-indexed (stupid!)
        marginal = mean(abs(template),1);
        marginal = marginal/max(marginal);
        
        probabilityNoise(clu,1) = trapz(marginal)/length(marginal); %Proportion of AUC      
        
%         spikeTimes = s.st(s.clu==cluID);
%         [xLin, nLin] = myACG(spikeTimes,[],[]);
%         nLin = nLin/max(nLin);
%         zeroACGTime(clu,1) = xLin(find(nLin>0.3,1,'first'));
    end

    writePhyTSV(ksDir, 'pNoise', templateIDs, probabilityNoise);
%     writePhyTSV(ksDir, 'zACGT', templateIDs, zeroACGTime);

end