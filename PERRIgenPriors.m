function [Up, Vp] = PERRIgenPriors(d, E, N, rW, S, Em, svFlag,f)
%% A function for generating priors from a pre-existing dataset
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The function creates a truncated version of the data in k-space. An LSVD
% transform then produces U_prior and V_prior.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%       d        - The k-space data from which the priors will be generated (Nk, Nr)
%       E        - The structure of the NUFFT transform into k-space        (xfm_nufft structure)
%       N.r      - Number of components to include in U_prior and V_prior   (scalar)
%       N.x, N.t - Size of the final image-space                            (scalars)
%       rW       - Radius of the tukey window (0 = no data allowed, 1 = to edge, 2.5 = no windowing) (scalar)
%       S.p      - The prior load string                                    (string)
%       S.d      - The load string of the original dataset                  (string)
%       Em       - The structure of the temporal mean NUFFT transform       (xfm_nufft structure)
%       svFlag   - Chose whether to save the output                         (boolean)
%       f        - A set of flag strutures, only used to pass through to the optimisation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs:
%       Up          - The spatial prior    (Nx , Nr)
%       Vp          - The temporal prior   (Nt , Nr)           
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by Harry Mason, M.Eng, University of Oxford
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial Parameters
if nargin<8; f.rmet  = 'rbyr';     % assume rank-by-rank-update generation. 'irpf' and 'full' are the other options. 
             f.print = false;      % whether we print the output
end  
             f.upsf  = false;      % regardless of what we're doing in the main thread, we can't do a psf construction here!
             f.vpsf  = false;
if nargin<7; svFlag  = 1;   end    % By default, save the output
if nargin<6; Em = E.mean;   end
if nargin<5; S.p     = 'g';        % Save string of prior
             S.d     = '4';        % Load string of the dataset (only used for truth-based priors)
end     
if nargin<4; rW      = 0.5; end     % Using half the distance as the radius is a good default
try
addpath MC_utils/misc Misc MC_utils/optim/transforms
end
%% Creating 1D Mask (to apply in radial k-space)
tukey_param = 0.4;
Nk    = E.Kd(1);          % Ensure k has twice as many points as x
kWin  = Window1D(Nk,rW,tukey_param);

%% Creating 2D mask ( to apply in cartesian k-space (non-lossy) )
Win2D = Window2D(Nk,rW,tukey_param);
Win2D = interp2(Win2D,linspace(1,Nk,sqrt(N.x)),linspace(1,Nk,sqrt(N.x))'); %Only really an issue at very small window sizes

%% Create orthogonal decomposition of dataset
if strcmpi(S.p(1),'g')
    disp('Generating k-t FASTER priors')
    
    %% Save String
    Sp = strcat(    regexprep(S.p,'^g(enerated)?(Priors)?','generatedPriors'),   '.mat');

    %% Applying 1D mask in radial k-space
    dP = d.*repmat(kWin,N.l,N.t,N.c); %Replicating the k-space window across every line, frame and coil
    DP = E'*dP;

    %% Initialise
    if exist(strcat(Sp(1:end-4),'Conv.mat'),'file')  
      
      fprintf('Converged version of prior exists, stopping recon now\n')
       [Up,Vp] = load_prior_n_data(Sp,'p');
    else
        try
            [Up, Vp] = init_UV(Sp,N.r,dP,E,Em); %If there's an existing prior, just plug that in
            fprintf('Starting from previous version of Prior\n')
        catch
            [Up, Vp] = init_UV('tr',N.r,dP,E,Em); %Can't build from previous version, so start anew
            fprintf('Can''t find previous prior, starting anew\n')
        end
    
        if f.print; fprintf('PERRIgenPriors:: Nl: %1.2f; Iters: %i; Prior ID: %s; Win_Rad: %2.2e; Tukey: %2.2e, Init Method: %s; Rank Method: %s\n', N.l, N.outer*N.r,S.p, rW, tukey_param, f.init, f.rmet); end

        %% Perform decomposition
        S.uv = Sp;
        mask = load_prior_n_data(S.d,'m');  % mask = repmat(mask,1,N.r); %Uncomment for backwards compatability
        Up   = Up.*mask;
        
        if contains(Sp,'Tik')
            lam = regexprep(Sp,'[a-zA-Z0-9]*Tik','');
            lam = regexp(lam,'[0-9]e?[+-]?[0-9]*','match');
            lam = str2double(lam);
        else
            lam = 0;
        end
        
        [Up, Vp] = PERRIoptm(dP, DP, E, Up,Vp, lam,Up*0,Vp*0, N, f,S,Win2D);    
    end
    DP       = Up*Vp';
    
elseif strcmp(S.p(1),'d') %use a gridded density-weighted recon    
    disp('Doing Density-Weighted Gridded recon')
    
    %% Save String
    Sp = strcat(    regexprep(S.p,'^d(ensity)?(Priors)?','densityPriors'),   '.mat');
       
    %% Applying 1D mask in radial k-space
    dP = d.*repmat(kWin,N.l,N.t,N.c); %Replicating the k-space window across every line, frame and coil
    DP = E'*dP;

elseif strcmp(S.p(1),'l') %use an iterative least squares densityweighted recon  
    disp('Doing LSQR gridded recon')
    
    %% Save String
    Sp = strcat(    regexprep(S.p,'^l(sqr)?(Priors)?','lsqrPriors'),   '.mat');
    
    %% Applying 1D mask in radial k-space
    dP = d.*repmat(kWin,N.l,N.t,N.c); %Replicating the k-space window across every line, frame and coil
    DP = E.iter(dP, @pcg, 1e-6,1e2); %1e-6 tol, 1e2 it
    
elseif strcmp(S.p(1),'0')||strcmp(S.p(1),'z') %use a prior of zeros (regularisation)
    disp('Using zero-filled Priors')
    
    %% Save String
    Sp = strcat(    regexprep(S.p,'^[0z]','zeroPrior'),   '.mat');
    
    %% Applying 1D mask in radial k-space
    Up = zeros(N.x,N.r);
    Vp = zeros(N.t,N.r); %So quick, no need to save
    return

elseif strcmpi(S.p(1),'r') %use a randomised prior
    disp('Making randomised priors')
    
    %% Save String
    Sp = strcat(    regexprep(S.p,'^r(andomised)?(Priors)?','randomisedPriors'),   '.mat');

    %% Load Data
    DP      = orth(randn(N.x,N.t) + randn(N.x,N.t)*1i); 
    
elseif strcmpi(S.p(1),'c') % using k-t PERRI for Time, Tik regularization for Space
    disp('Making chronos (time) priors')
    
    %% Save String
    Spg = strcat(    regexprep(S.p,'^c(hronos)?(Priors)?','generatedPriors'),   '.mat');
    Sp  = strcat(    regexprep(S.p,'^c(hronos)?(Priors)?','chronosPriors'),   '.mat');

    %% Load Data
    if isfile(Spg); [Up,Vp] = load_prior_n_data(Spg,'p'); 
    else;           S.p = Spg;
                    [Up, Vp] = PERRIgenPriors(d, E, N, rW, S, Em, svFlag,f);
    end
    
    Up = Up*0;
    
    return
    
elseif strcmpi(S.p(1),'s') % % using k-t PERRI for Space, Tik regularization for Time
    disp('Making space priors')
    
    %% Save String
    Spg = strcat(    regexprep(S.p,'^s(pace)?(Priors)?','generatedPriors'),   '.mat');
    Sp  = strcat(    regexprep(S.p,'^s(pace)?(Priors)?','spacePriors'),   '.mat');

    %% Load Data
    if isfile(Spg); [Up,Vp] = load_prior_n_data(Spg,'p'); 
    else;           S.p = Spg;
                    [Up, Vp] = PERRIgenPriors(d, E, N, rW, S, Em, svFlag,f);
    end
    
    Vp = Vp*0;
    
    return

elseif strcmpi(S.p(1),'t') %use a windowed version of the truth  
    disp('Making truth-based priors')
    %% Save String
    Sp = strcat(    regexprep(S.p,'t(ruth)?(Priors)?','truthPriors'),   '.mat');
    
    %% Load Data
    if f.k
        dP = d.*repmat(kWin,N.l,N.t,N.c); %Replicating the k-space window across every line, frame and coil
        DP = E.iter(dP, @pcg, 1e-6,1e2); %may be no underlying truth for prospective (may rewrite later for oversampled versions)
    else
    [~,DP]      = load_prior_n_data(S.d,'d'); 
    end 
else
    try
    [Up, Vp] = load_prior_n_data(S.p,'p');
    DP = Up*Vp';

    disp('Loading pre-existing priors and windowing')
    S.p   = regexprep(S.p,'rW(\d{4})?',sprintf('rW%04d',round(rW*1e3)));
    Sp = S.p;
    catch
    disp('No valid prior selected')
    return
    end
    svFlag = 0; %Don't really want to save over it in this case.
end

%% Post creation windowing (to ensure empty k-space outside window)
DP       = Window2Ddata(DP,Win2D); % function to auto apply the window in k-space and then convert back
DP       = Window2Ddata(DP); %Also make sure there's nothing outside the outer ring of k-space (some radii allow that)
[Up, Vp] = lsvdn( DP, N.r,'S');

%% Final save

if svFlag
    fprintf('Saving Prior as %s\n',Sp);
    save(Sp,'Up','Vp')
end
end

