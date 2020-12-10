function PERRIzstat(varargin)

%% The core centre of the k-t PERRI algorithm. Input the basic arguments, and this algorithm will take care of the rest.
%
% In: S.d, N.l, lam, S.p, rW, f.init, N.inner, N.outer, N.r, f.ld.p, f.ld.m, f.rmet, S.uv, S.num
% Saves nii.gz file with z'rest_of_name', similar to UV'rest_of_name that is typically used
% It then calls fslScript and performs a FEAT analysis.

%% Sorting Inputs
if nargin < 1; S.d     = 2.0;   else; S.d    = varargin{1}; end %ID of data to load
if nargin < 2; N.l     =  15;   else; N.l    = varargin{2}; end %number of sampling lines
if nargin < 3; S.uv    ='UV';   else; S.uv   = varargin{3}; end %Save String
if nargin < 4; S.num   = [] ;   else; S.num  = varargin{4}; end %Save Number modifier
if nargin < 5; f.w     = true; else; f.w    = varargin{5}; end

addpath MC_utils/misc Misc 
if isfolder('/usr/local/fsl/etc/matlab')
    addpath /usr/local/fsl/etc/matlab
else
    addpath /opt/fmrib/fsl/etc/matlab
end
%% Save locations
if isnumeric(S.num); S.num = num2str(S.num); end
S.uv   = regexprep(S.uv,'.mat',strcat(S.num,'.mat')); %adds qualifying number
S.z = fliplr(regexprep(fliplr(S.uv),'(VU)?/','z/','once')); %The flip protects any folder name, but ensures it starts with a z
S.z = regexprep(S.z,'.mat','');
S.z = strcat('FEAT/',fliplr(regexprep(fliplr(S.z),'\/.*','','once'))); %removes all other folders, places in FEAT folder
S.uv = regexprep(S.uv,'.*(?<!.mat)$',strcat(S.uv,'.mat')); %Make sure it ends as a .mat file


%% Load file
load(S.uv,'U','V');
if ~exist('U','var')
    load(S.uv,'Up','Vp');
    U=Up;
V=Vp;
end

%% Load Mask
mask = load_prior_n_data(S.d,'m');

%% Initialise values
N.x = size(U,1);
N.t = size(V,1);

[~,f.k] = load_prior_n_data(S.d,'dID'); %Throws up a flag if the data is retrospective (and not if it's prospective)
[spf] = load_prior_n_data(S.d,'TR');
if f.k
    spf = spf*N.l; %seconds per frame - each line took roughly 0.2 seconds
end

D = U*V';

if f.w

D =  Window2Ddata(U*V').*mask;

S.z = strcat(S.z,'w');
end

D = abs(D)./max(abs(D(:)));

%% Create nii.gz file
i1 = 1;

% We found that save_avw would sometimes save a file that MATLAB couldn't
% open. This loop basically allows multiple attempts. Given we found an
% approximately 10% chance of failure, allowing 10 attempts pretty much
% guarantees success, if the odd saving is the thing which is wrong.
while i1 < 10
    %save_avw(reshape(abs(Data),#x,#y,#z,#t]),?filename','f?,[xcm ycm zcm secondsperframe]) maybe xmm ymm zmm spf
    save_avw( reshape(   abs(D),[sqrt(N.x) sqrt(N.x) 1 N.t]  ) ,  S.z,  'f',  [2 2 2 spf])
    
    try
    % Check a legitimate file has been loaded
    read_avw(S.z);
    i1 = 11; %If it succeeds, we're good
    catch
    i1 = i1+1 %Otherwise try again
    end
end

%save_avw(reshape(abs(Data),#x,#y,#z,#t]),?filename','f?,[xcm ycm zcm secondsperframe]) maybe xmm ymm zmm spf
save_avw( reshape(   abs(D),[sqrt(N.x) sqrt(N.x) 1 N.t]  ) ,  S.z,  'f',  [2 2 2 spf])

fprintf('Saved as %s\n',S.z)

%% Carry out analysis
%system(sprintf('./zScript %s',S.z)) %fine for now - will neeed a way to get a new design.fsf file if we start analysing retrospective data/the phantom dataset.

%final mv *.mix ./FEAT? Considering it may have moved everything to the dataHoldingPen, this could be tricky


