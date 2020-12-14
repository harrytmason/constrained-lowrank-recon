function [E, Em, effAcc, fsRad] = CLRgenE(Sd, Nl, fk)
%% CLRgenE - a combined sampling & Fourer transform algorithm
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
%
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
else                                
    [Nx, Ny, Nt] = size(X);           
    ppl = Nx*2;
end

% load coils
try     c = load_prior_n_data(Sd,'s');
catch;  c = [];
end


%% Creating k indices
try
    [k,   fkind] = load_prior_n_data(Sd,'k');
    [ang, fang ] = load_prior_n_data(Sd,'a');
catch
    fkind = false;
    fang = false;
end
if fkind
    k = k(:,1:Nl*Nt,:);
else
    gAng = 360/(1+sqrt(5));   % The amount each angle will increment in degrees (golden angle)
    if ~fang
        ang = mod(0:gAng:(gAng*(Nl*Nt-1)),360); %don't have to use mod(angle,360), I'm doing it for ease in understanding when you read the matrix
    else
        fprintf('Using predetermined angles\n')
        if max(abs(ang))<7 && max(diff(ang))>0 %angles likely given in radians, and I want degrees.
            ang = ang/pi*180;
        end
            ang = ang(1:Nl*Nt); % allow different Nl per frame modifier here.
    end
    k = zeros(ppl,Nl*Nt,2);
    k(:,:,1) = bsxfun(@times,linspace(-pi,pi,ppl)',cosd(ang));
    k(:,:,2) = bsxfun(@times,linspace(-pi,pi,ppl)',sind(ang)); 
end

%% Print Out
effAcc = (Ny/Nl) * (pi/2);                % Effect acceleration
fsRad  = Nl/(Ny*pi/2);                    % Radius at which k-space is still full sampled
fprintf('Using %d sample lines per frame. Effective acceleration = %1.3f. Fully Sampled Radius = %1.3f.\n',Nl,effAcc,fsRad);

%% Create structures
fprintf('Generating NUFFT structures\n')
if ~fk;         E  = xfm_NUFFT([Nx Ny 1 Nt],c,[],reshape(k,[],Nt,2)); %second parameter is coil sensitivity. dsize = [Nk, Nt, Ncoils]. c = [] if there are no coils
else;           E  = xfm_NUFFT([Nx Ny 1 Nt],c,[],reshape(k,[],Nt,2),'wi',1); %'wi' 1 applies density weighting initially, to recognise that the original data is in k-space
end

if nargout>1 %A lot of computing power is saved by not computing the temporal mean structure if it isn't asked for
    if ~fk;     Em = xfm_NUFFT([Nx Ny 1  1],c,[],reshape(k,[], 1,2),'table',true); % Using table interpolation to make the storage size a bit smaller, otherwise the files get so big they can make things run slow.
    else;       Em = xfm_NUFFT([Nx Ny 1  1],c,[],reshape(k,[], 1,2),'table',true,'wi',1); % Using table interpolation to make the storage size a bit smaller, otherwise the files get so big they can make things run slow.
    end
end
fprintf('NUFFT structures generated\n')


