function [ tdc, tlc, tic, instant ] = time_dependent_intrinsic_corr( imfs1, imfs2, bootstrap, plots, prints, subranges, periodn, alpha )
    %TIME_DEPENDENT_INTRINSIC_CORR Calculate time-dependent intrinsic correlation
    %   
    %   This is proxy-function for short-named routine tdic
    %   
    %   Copyright (c) 2014-2016 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0  2014.04.05: initial version
    %       1.1  2015.04.13: added time-dependent local correlation, refactor internal routines
    %       1.11 2015.12.25: added compatibility with Octave environment
    %       1.2 2016.12.11: the code moved to short-titled tdic routine
    %   
    
    [tdc, tlc, tic, instant] = tdic(imfs1, imfs2, bootstrap, plots, prints, subranges, periodn, alpha);
    
end