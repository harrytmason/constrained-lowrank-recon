function [E, Em, effAcc, fsRad, maxR] = PERRIgenE(Sd, Nl, fk,start_frame)
%% PERRIgenE - a combined sampling & Fourer transform algorithm
%
% This function mainly generates E, an xfm_nufft structure which converts a
% dataset between k-space an image space (K = E*X, X=E'*K)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs :  Sd:    The string of the desired dataset. This input could be a (Nk * Nt) or (Nx * Nt) matrix as well, with the size dependent on fk
%           Nl:    Number of sampling lines per frame. Sample (2*SL)*(Nx*Nt) of k-t space. (scalar)
%           fk:    if true, X is in k-t space. if false, X should be in x-y-t space.       (boolean)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs : E:     the NUFFT structure, used to convert between k-space.
%           Em:    the NUFFT structure of the temporal mean.
%           effAcc:the effective acceleration, given the size and number of sample lines
%           fsRad: the radius at which the k-space is still fully sampled (potentially use as radius to geenrate priors)           
%           maxR:  the maximum rank to possibly consider
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Notes : Requires M. Chiew's optim/transforms and related irt toolbox
%         Only needs regenerating if dataset or Nl changes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Harry Mason, University of Oxford, 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialise General Variables
% sets values for X, Nl, fk
%
if nargin < 4; start_frame = 0; end
if nargin < 3; fk = false; end % Are we working with k-space data?
if nargin < 2; Nl = 1;     end % Assume trivial transform

if ischar(Sd) || isscalar(Sd)  % Intended usage - pass through data identifier (use with load_prior_n_data.m)
[X] = load_prior_n_data(Sd,'d');
[~, fk] = load_prior_n_data(Sd,'dID'); %lets you know if it's a k-t space input
else
   X = Sd; %User has plugged a dataset straight into the function, rather than an identifier
end

%% Initialising Sizes
% sets values for Nx, Ny, Nt, Nk, ppl (points per line)
%
if fk
    ppl = size(X,1);   
    Nx = ppl/2; Ny = ppl/2;
    Nt = floor(size(X,2)/Nl); 
    Nc = size(X,3); %Number of coils
else                                
    [Nx, Ny, Nt] = size(X);           
    ppl = Nx*2;
    Nc = 1;
end

% load coils
try     c = load_prior_n_data(Sd,'s');
catch;  c = [];
end



%% Creating k indices
try
    [k,   fkind] = load_prior_n_data(Sd,'k');
    [ang, fang ] = load_prior_n_data(Sd,'a');
    fcart = islogical(k); %Nl is irrelevant here, Cartesian patterns must be pre-defined
catch
    fkind = false;
    fang = false;
    fcart = false;
end
if fkind
    if ~fcart
    k = k(:,1:Nl*Nt,:);
    end
    ppl = size(k,1);
else
    gAng = 360/(1+sqrt(5));   % The amount each angle will increment in degrees (golden angle)
    if ~fang
        ang = mod(0:gAng:(gAng*(Nl*Nt-1)),360); %don't have to use mod(angle,360), I'm doing it for ease in understanding when you read the matrix
    else
        fprintf('Using predetermined angles\n')
        if max(abs(ang))<7 %angles likely given in radians, and I want degrees.
            ang = ang/pi*180;
        end
            ang = ang(1:Nl*Nt); % allow different Nl per frame modifier here.
    end
    k = zeros(ppl,Nl*Nt,2);
    k(:,:,1) = bsxfun(@times,linspace(-pi,pi,ppl)',cosd(ang));
    k(:,:,2) = bsxfun(@times,linspace(-pi,pi,ppl)',sind(ang)); 
end

Nk = ppl*Nl; %Number of points in k per frame

%%%%%%%%%%%%%%%
%k=fftshift(k,2); %temporary line, REMOVE
%%%%%%%%%%%%%%
%fs=1;

%% Print Out
effAcc = (Ny/Nl) * (pi/2);                % Effect acceleration
fsRad  = Nl/(Ny*pi/2);                    % Radius at which k-space is still full sampled
maxR = ceil(((Nx*Ny+Nt) - sqrt((Nx*Ny+Nt)^2-4*(Nk*Nt)))/2);   % Maximum rank (nothing useful beyond here)
fprintf('Using %d sample lines per frame. Effective acceleration = %1.3f. Fully Sampled Radius = %1.3f. Theoretical maximum rank = %d.\n',Nl,effAcc,fsRad,maxR);

%% Trimming frames
if start_frame > 0 & ~fcart
    
    [~,temp] = load_prior_n_data(Sd,'t');
    Nt_true = size(temp,1); %Actual number of time frames we want/Number of time frames in the true reconstruction
    
    if Nt_true >= Nt %Warning to not allow reduced DoF where it isn't possible to do so.
          fprintf('Truth has %d DoF, greater/equal to the proposed recon DOF: %d, so can''t use reduced temporal DoF recon\n',Nt_true, Nt);
          return
    end
    
    if Nt_true>1
        Nl_true     = Nt*Nl/Nt_true;              % Spacing between the start of frames, as determined by the number of lines used in the true reconstruction
        
        if start_frame>(Nl_true-Nl)
            fprintf('Start frames available are 1 to %d. User selected %d, which is outside of the allowed range\n',Nl_true-Nl,start_frame);
            return
        end
        
        start_index = start_frame:Nl_true:Nt*Nl;   % Index of the start line of each frame
        index       = vec( start_index + (1:Nl)' ); % All the lines needed for the recon
        k           = k(:,index,:);               % Resizing k
        Nt          = size(k,2)/Nl;                    % New Number of frames
        fprintf('PERRIgenE:: Doing reduced temporal DoF recon, Nl: %d, Nt: %d, start_frame: %d\n',Nl,Nt,start_frame)
    else
        fprintf('PERRIgenE:: No underlying Truth, can''t use reduced temporal DoF recon\n')
    end
end

%% Create structures
if ~fcart

if ~fk;         E  = xfm_NUFFT([Nx Ny 1 Nt],c,[],reshape(k,[],Nt,2)); %second parameter is coil sensitivity. dsize = [Nk, Nt, Ncoils]. c = [] if there are no coils
else;           E  = xfm_NUFFT([Nx Ny 1 Nt],c,[],reshape(k,[],Nt,2),'wi',1); %'wi' 1 applies density weighting initially, to recognise that the original data is in k-space
end

else
    E  = xfm_FFT([Nx Ny 1 Nt],c,[],'mask',reshape(k,[Nx Ny 1 Nt]));
   % E.Kd=E.Nd;
end


if nargout>1 %A lot of computing power is saved by not computing the temporal mean structure if it isn't asked for
    if ~fcart
    if ~fk;     Em = xfm_NUFFT([Nx Ny 1  1],c,[],reshape(k,[], 1,2),'table',true); % Using table interpolation to make the storage size a bit smaller, otherwise the files get so big they can make things run slow.
    else;       Em = xfm_NUFFT([Nx Ny 1  1],c,[],reshape(k,[], 1,2),'table',true,'wi',1); % Using table interpolation to make the storage size a bit smaller, otherwise the files get so big they can make things run slow.
    end
    else; Em  = xfm_FFT([1 1 1 1],c,[],'mask',reshape(k(1,1,1,1),[1 1 1 1]));
       %   Em.Kd=Em.Nd;
    end
end



