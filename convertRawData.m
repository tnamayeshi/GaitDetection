function convertRawData(ddir, mass, subName, domLeg, ptype, pleg, ponset, ptime, intcty, tnum)

    mot = convertMot2Struct(fullfile(ddir, 'Trial.mot'));
    trc = convertTrc2Struct(fullfile(ddir, 'Trial.trc'));
    
    data = struct();
    data.mass = mass;
    data.subject = subName;
    data.dominantLeg = domLeg;
    data.mot = mot;
    data.trc = trc;
    data.dataType = 'trial';
    data.perturbationType = ptype;
    data.perturbationTime = ptime;
    data.perturbationOnset = ponset;
    data.perturbedLeg = pleg;
    data.intensity = intcty;
    data.trial = tnum;
    
    save(fullfile(ddir, 'rawdata.mat'), '-struct', 'data');
end