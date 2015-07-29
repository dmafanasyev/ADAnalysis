function [ h, amplitude, phase, omega, period ] = hilbert_transform( data, smooth, plots )
    %HILBERT_TRANSFORM Perform Hilbert transform with or without smoothing of phase.
    %   Zero-phase FIR (Finite impulse response) filter with averaging (equals denominator coefficients) used 2 times for phase smoothing.
    %   Source data can be reconstructed through data = real((sum(amp.*exp(i*phase), 2)));
    %
    %   Input:
    %       data - source data matrix
    %       smooth - boolean flag (0 or 1) for phase smoothing
    %       plots - boolean flag (0 or 1) for graphics plot
    %
    %   Output:
    %       h - analytic signal
    %       amplitude - instantaneous amplitude
    %       phase - instantaneous phase
    %       omega - instantaneous frequency
    %       period - instantaneous period
    %
    %   Copyright (c) 2014 by Dmitriy O. Afanasyev
    %   Versions:
    %       1.0 2014.04.05: initial version
    %       
    
    if (nargin < 2)
        smooth = 0;
    end
    if (nargin < 3)
        plots = 0;
    end
    
    h = hilbert(data);
    amplitude = abs(h);
    phase = unwrap(angle(h));
    
    if(smooth > 0)
        % smooth phase by applying zero-phase FIR (Finite impulse response) filter with averaging (equals denominator coefficients) 2 times
        denominator = [1/3 1/3 1/3];
        numerator = 1;
        phase = filtfilt(denominator, numerator, phase);
        phase = filtfilt(denominator, numerator, phase);
    end
    
    omega = (1/(2*pi)) * abs(diff(phase));
    omega = [omega(1,:); omega];
    
    period = 1./omega;
    
    if(plots > 0)
        figure;
        subplot(2, 2, 1);
        plot(amplitude);
        
        subplot(2, 2, 2);
        plot(omega);
        
        subplot(2, 2, 3);
        plot(phase);
        
        subplot(2, 2, 4);
        plot(period);
    end
end
