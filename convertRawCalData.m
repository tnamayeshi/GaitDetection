function convertRawCalData(ddir, mass, subName, domLeg, staticDuration)

    mot = convertMot2Struct(fullfile(ddir, 'Calibration.mot'));
    trc = convertTrc2Struct(fullfile(ddir, 'Calibration.trc'));
    
    data = struct();
    data.mass = mass;
    data.subject = subName;
    data.dominantLeg = domLeg;
    data.mot = mot;
    data.trc = trc;
	data.noMoveTime = staticDuration;
    data.dataType = 'cal';  
	
    save(fullfile(ddir, 'caldata.mat'), '-struct', 'data');
end