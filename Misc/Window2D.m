function Win2D = Window2D(Nx,rW,tukey_param,radius_choice)

%% Win2D.m
%
% Create a 2D window to be applied in Cartesian k-space
% Will have a circular radius of rW (1 = full width)
% tukey_param specifies cutoff at edge
%
% Inputs: Nx - Size x of x*x square image
%         rW - radius of window (>2.5 = no filter, 0 = full filter) 
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

if rW <= 0;
    Win2D = zeros(Nx,Nx);
    return
end

if strcmpi(radius_choice,'FWHM')
    xRad = rW/(1-tukey_param/2); %The new window radius, as the input defines the mid-point of the slope ( rW*(1-tukey_param/2) is midpoint if we don't redefine ) 
    
elseif strcmpi(radius_choice,'Full')||strcmpi(radius_choice,'Edge')
    xRad = rW; %Here, our radius is also the input
else
    disp('Please choose a valid option of ''Edge'' or ''FWHM'' for the final input, to indicate whether the window radius (input 2) is the full radius or radius to the FWHM')

end



%% Create 1D Window

if mod(Nx,2)                                                            %Nx is odd, so we need an odd shaped window
    xWid = 2 * (round(Nx*xRad/2 + 0.5) - 0.5 );                            %Diameter of the window
else                                                                    %Nx is even, so we need an even shaped window
    xWid = 2 *  round(Nx*xRad/2);                                          %Diameter of the window
end

xWin = tukeywin(xWid,tukey_param);

%% Interpolate to 2D window

try %The proper way
    Win2D = zeros(Nx);
    iiWid = linspace(-1,1,Nx);
    jjWid = linspace(-1,1,Nx); % Slightly redundant, but keeps the two dimensions clear
    
    for ii = 1:Nx
    for jj = 1:Nx
        r     = sqrt( iiWid(ii)^2 + jjWid(jj)^2 );                                 % Radius
        Win2D(ii,jj) = interp1(linspace(-xRad,xRad,xWid), xWin, r);    % Equivalent value in mask
    end
    end
    
    
catch % Create bottom-right corner of the mask, then rotate it (Proper way above fills the whole thing properly, I just like this solution and it makes the interpolation more intuitive, so I'm keeping it as reference).
    Win2D = zeros(Nx/2); 

    for ii = 1:Nx/2
    for jj = 1:Nx/2
        r     = sqrt( ii^2 + jj^2 );                                 % Radius
        Win2D(ii,jj) = interp1((1:xWid/2), xWin(end/2+1:end), r);    % Equivalent value in mask
    end
    end
    
    Win2D = [rot90(Win2D,2) rot90(Win2D,1); rot90(Win2D,3) Win2D]; %Rotate to fill space

end

%% Zeroing
Win2D(isnan(Win2D)) = 0;