function [ period ] = period_zero_cross( data )
    %PERIOD_ZERO_CROSS Calculate period estimation using zero-crossing method.
    %
    %   Copyright (c) 2014 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0 2014.04.12: initial version
    %
    
    [indmin, indmax, indzer] = extr(data');
    period = 4 * (size(data, 1) / (size(indmin, 2) + size(indmax, 2) + size(indzer, 2)));
end

