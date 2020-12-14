function Win1D = Window1D(Nk,rW,tukey_param,radius_choice)

%% Win2D.m
%
% Create a 1D window to be applied in k-space
% Will have a radius of rW (1 = full width)
% tukey_param specifies cutoff at edge
%
% Inputs: Nk - Length of whole window
%         rW - radius of window in relative terms (>2.5 = no filter, 0 = full filter) 
%         tukey_param - slope at edge of window
%         radius_choice - 'FWHM' rW defines Full Width at Half Maximum (over 2)
%                         'Full' rW defines Full radius of window
%
% Outputs: Win2D - a 2D windowing mask
%
% Created by Harry Mason, University of Oxford

%% Initialisation
if nargin<4
    radius_choice = 'FWHM';
end

if rW <= 0
    Win1D = zeros(Nk,1);
    return
end

%% Definining Relative Width of window
if strcmpi(radius_choice,'FWHM')
    kRad = rW/(1-tukey_param/2); %The new window radius, as the input defines the mid-point of the slope ( rW*(1-tukey_param/2) is midpoint if we don't redefine ) 
    
elseif strcmpi(radius_choice,'Full')||strcmpi(radius_choice,'Edge')
    kRad = rW; %Here, our radius is also the input
else
    disp('Please choose a valid option of ''Edge'' or ''FWHM'' for the final input, to indicate whether the window radius (input 2) is the full radius or radius to the FWHM')

end

%% Defining actual discrete width
if mod(Nk,2)                                                            %Nx is odd, so we need an odd shaped window
    kWid = 2 * (round(Nk*kRad/2 + 0.5) - 0.5 );                            %Diameter of the window
else                                                                    %Nx is even, so we need an even shaped window
    kWid = 2 *  round(Nk*kRad/2);                                          %Diameter of the window
end

%% Create window
Win1D = tukeywin(kWid,tukey_param);

%% zeros-pad/cut window to right size
if Nk>kWid                                      % zero-pads windows which are too small
    Win1D = padarray(Win1D,(Nk-kWid)/2); 
else                                            % It will fail if we want a bigger radius than 1, this will solve that (cuts off the edge of the window).
    Win1D = Win1D(1+(kWid-Nk)/2:end-(kWid-Nk)/2);
end
