%--------------------------------------------------------------------------
%
% The external electric field is assumed to be a sequence of (overlapping 
% or non-overlapping) pulses, thereby mimicking modern time-resolved 
% spectrosopic experiments employing (ultra-)short laser pulses.  
%
% Each pulse has a constant carrier frequency modulated by an envelope of 
% different shape. Furthermore it is characterized by the field amplitude,
% polarization, time delay (center of pulse), duration(FWHM), and phase
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function init
global time

%% Set the number of pulses to default values and start output
if isfield(time, 'pulse')

    if ~isfield(time,'efield') % If not already existing
        time.efield = tmp.efield; % Construct object
    end
    
    log.disp ('***************************************************************')
    log.disp ('External electric field as a sequence of pulses')
    log.disp ('***************************************************************')
    log.disp ( ' ' )
    log.disp ( [ 'Using (Floquet) dressed states  : ' int2str(time.efield.dressed) ] )
    log.disp ( [ 'Using complex-valued fields     : ' int2str(time.efield.complex) ] )
    log.disp ( [ 'Number of pulses                : ' int2str(length(time.pulse) ) ] )
    log.disp ( ' ' )
        
    %% Set up individual pulses
    for p=1:length(time.pulse)
        
        log.disp ('***************************************************************')
        log.disp ( [ 'Parameters of electric field pulse : ' int2str(p) ] )
        
        if ~isempty (time.pulse{p})
            time.pulse{p}.ind = p; % Tell each pulse its index
            init_efi ( time.pulse{p} ); % Initialize pulse
            disp     ( time.pulse{p} ); log.disp ( ' ' ) % Display pulse
        else
            log.disp ( 'Not available' )
            log.disp ('***************************************************************')
            log.disp ( ' ' )
            
        end
    end
    
    % Initialize fields
    init (time.efield);
    
    % Optionally switch to Floquet dressed
    if time.efield.dressed
        floquet (time.efield);
    end
    
    
else
    
    log.disp ('***************************************************************')
    log.disp ('No electric field available                                    ')
    log.disp ('***************************************************************')
    log.disp (' ')
end

end




