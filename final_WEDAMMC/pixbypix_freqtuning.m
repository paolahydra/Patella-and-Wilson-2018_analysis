function peakPrefFreq = pixbypix_freqtuning(pxk, pp, tn)

windowResp = 3:7;

%% intermediate output will be a matrix (nSignPixels x 19)
matrix = nan(sum(pxk), 19);

% tone 4Hz: simply calculate percent dff of max reponse
resp = tn.responses(pxk,windowResp, 1);
basl = tn.baselines(pxk,:,1);
matrix(:,1) =  100 * ( ( max(resp,[],2) -  nanmean(basl,2) ) ./ nanmean(basl,2) );

% all pips: simply calculate percent dff of max reponse, and concatenate
resp = pp.responses(pxk,windowResp, :);
basl = pp.baselines(pxk,:,:);
matrix(:,2:end) =  100 * squeeze( ( max(resp,[],2) -  nanmean(basl,2) ) ./ nanmean(basl,2) );

%% for each ipxl, find peak index
[peakPrefFreq.dff , peakPrefFreq.idx] = max(matrix, [], 2);

end