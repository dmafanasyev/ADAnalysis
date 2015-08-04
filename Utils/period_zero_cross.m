function [ period, indmin, indmax, indzer ] = period_zero_cross( data )
    %PERIOD_ZERO_CROSS Calculate period estimation using zero-crossing method.
    %
    %   Copyright (c) 2014 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0 2014.04.12: initial version
    %       1.1 2015.08.04: added output params for indexes of min, max & zero
    %
    
    [indmin, indmax, indzer] = extr(data');
    period = 4 * (size(data, 1) / (size(indmin, 2) + size(indmax, 2) + size(indzer, 2)));
end

