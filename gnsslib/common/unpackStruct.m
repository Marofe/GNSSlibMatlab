function unpackStruct (structure)
    fn = fieldnames(structure);
    
    for i = 1:numel(fn)
        fni = string(fn(i));
        field = structure.(fni);
        if (isstruct(field))
            unpackStruct(field);
            continue;
        end
        assignin('caller', fni, field);
    end
end