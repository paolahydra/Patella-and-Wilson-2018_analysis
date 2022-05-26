% make ts_dec_all timestamps for a given subset of stimulusTypes (and
% stimulusClasses)
c = 1;
sortedStimNs_thisclass = sortedStimNs{c};
effectiveFirst  = T(sortedClassNs(c)).ts_dec( find(T(sortedClassNs(c)).ts_dec >= cropbaseline, 1, 'first') );
effectiveLast   = T(sortedClassNs(c)).ts_dec( find(T(sortedClassNs(c)).ts_dec <= T(sortedClassNs(c)).stimDur + 10, 1, 'last') );
iKeepResponse = T(sortedClassNs(c)).ts_dec >= effectiveFirst   & T(sortedClassNs(c)).ts_dec <= effectiveLast;
ts_dec = T(sortedClassNs(c)).ts_dec(iKeepResponse) ;
ts =  T(sortedClassNs(c)).ts; 
ts(length(T(sortedClassNs(c)).sensor) +1 : end) = []; %fix for courtshp songs messed up
deltaEndTss = ts(end) - T(sortedClassNs(c)).ts_dec(end);
deltaStartTss = -(cropbaseline - ts_dec(1));

ts_dec_all = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
zero_tr_idxs = find(ts_dec>=0,1);
stN_tr = sortedStimNs_thisclass(1) *ones(length(ts_dec), 1);

for i = 2 : length(sortedStimNs_thisclass) %assuming they ally belong to the same or analogous class
    zero_tr_idxs = cat(1, zero_tr_idxs, find(ts_dec>=0,1) + length(ts_dec_all));
    ts_dec_all = cat(1, ts_dec_all, ts_dec - ts_dec(1) + ts_dec_all(end) + deltaStartTss + deltaEndTss );
    stN_tr = cat(1, stN_tr, sortedStimNs_thisclass(i)*ones(length(ts_dec), 1) );
end
    
for c = 2:length(sortedClassNs)
    sortedStimNs_thisclass = sortedStimNs{c};
    effectiveFirst  = T(sortedClassNs(c)).ts_dec( find(T(sortedClassNs(c)).ts_dec >= cropbaseline, 1, 'first') );
    effectiveLast   = T(sortedClassNs(c)).ts_dec( find(T(sortedClassNs(c)).ts_dec <= T(sortedClassNs(c)).stimDur + 10, 1, 'last') );
    iKeepResponse = T(sortedClassNs(c)).ts_dec >= effectiveFirst   & T(sortedClassNs(c)).ts_dec <= effectiveLast;
    ts_dec = T(sortedClassNs(c)).ts_dec(iKeepResponse) ;
    ts =  T(sortedClassNs(c)).ts;
    ts(length(T(sortedClassNs(c)).sensor) +1 : end) = [];
    deltaEndTss = ts(end) - T(sortedClassNs(c)).ts_dec(end);
    deltaStartTss = -(cropbaseline - ts_dec(1));
    
    for i = 1 : length(sortedStimNs_thisclass) %assuming they ally belong to the same or analogous class
        zero_tr_idxs = cat(1, zero_tr_idxs, find(ts_dec>=0,1) + length(ts_dec_all));
        ts_dec_all = cat(1, ts_dec_all, ts_dec - ts_dec(1) + ts_dec_all(end) + deltaStartTss + deltaEndTss );
        stN_tr = cat(1, stN_tr, sortedStimNs_thisclass(i)*ones(length(ts_dec), 1) );
    end
end