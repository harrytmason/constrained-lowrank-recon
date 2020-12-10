function [U,V,C] = PERRIcore(varargin)

%% The core centre of the k-t PERRI algorithm. Input the basic arguments, and this algorithm will take care of the rest.
%
% In: S.d, N.l, lam, S.p, rW, f.init, N.inner, N.outer, N.r, f.ld.p, f.ld.m, f.rmet, S.uv, S.num
% Out: U (the spatial components), V (the temporal components), C (the costs)

%% Sorting Inputs
if nargin < 1; S.d     = 2.0; else; S.d    = varargin{1}; end %ID of data to load
if nargin < 2; N.l     =  15; else; N.l    = varargin{2}; end %number of sampling lines
if nargin < 3; lam     = 0.0; else; lam    = varargin{3}; end %lambda values [lambda_u lambda_v]
    lamu = lam(1);
    if numel(lam)==1
        lamv = lam(end); %if only one value used, assume lamu==lamv
        lams = 0;
    elseif numel(lam)>1
        lamv = lam(2);
        if numel(lam)>2
            lams = lam(3);
        else
            lams = 0;    % if less than 2 values speicfied, assume no smoothing
        end
    end


if nargin < 4; S.p     = 'g'; else; S.p    = varargin{4}; end %ID of prior to load
if nargin < 5; rWin    = 0.5; else; rWin   = varargin{5}; end %Size of Tukey Window
if nargin < 6; f.init  = 't'; else; f.init = varargin{6}; end %Method of Initialisation of Uo amd Vo
    if (lamu == inf && lamv == inf) || strcmpi(f.init,'uvpsf')
        lamu = 0; lamv=0; f.vpsf = true;  f.upsf = true; %f.init = 'uvpsf';
    elseif lamu == inf || strcmpi(f.init,'upsf')
        lamu = 0;         f.vpsf = false; f.upsf = true; %f.init = 'upsf';
    elseif lamv == inf || strcmpi(f.init,'vpsf')
        lamv = 0;         f.vpsf = true;  f.upsf = false;%f.init = 'vpsf';
    else;                 f.vpsf = false; f.upsf = false;
    end
if nargin < 7; N.inner =  10; else; N.inner= varargin{7}; end %#iterations per Rank per rank
if nargin < 8; N.outer =  10; else; N.outer= varargin{8}; end %#iterations per sub-problem
if nargin < 9; N.r     = 8.0; else; N.r    = varargin{9}; end %End Rank
if nargin <10; f.ld.p  =false;else; f.ld.p = varargin{10};end %Load previous priors
    if numel(f.ld.p)==2; f.ld.m = f.ld.p(2); f.ld.p = f.ld.p(1);  %allow load mask as a array in load prior, if so desired.
    else;                f.ld.m = false;
    end
if nargin <11; N.noise = 0.0; else; N.noise= varargin{11};end %Load previous structures
if nargin <12; f.rmet ='full';else; f.rmet = varargin{12};end %Save Number modifier
if nargin <13; S.uv    ='UV'; else; S.uv   = varargin{13};end %Save String
    S.uv = regexprep(S.uv,'.*(?<!.mat)$',strcat(S.uv,'.mat')); %Make sure it ends as a .mat file
if nargin <14; S.num   = [] ; else; S.num  = varargin{14};end %Save Number modifier
    if isnumeric(S.num); S.num = num2str(S.num); end

%% Basics for whole function
f.print = true;
rng('shuffle'); %ensure that wherever this is run, we have a unique randomisation seed
try %addpath causes problem during runtime instantiations of function, when all required functions are autocompiled anyway.
addpath MC_utils/optim/transforms MC_utils/misc Misc irt
catch
end
setup; clear tmp irtdir;

%% Save locations
S.uv   = regexprep(S.uv,'.mat',strcat(S.num,'.mat')); %adds qualifying number
S.temp = regexprep(S.uv,'(DataHoldingPen/)?UV','UselessDataHoldingPen/temp'); %Save tempfiles in the UselessDataHoldingPen (create this file)
fprintf('File save: %s\n', S.uv)
[~,start_frame] = load_prior_n_data(S.d,'ID'); %At some point, make this a separate input

%% Gen flag
[~,f.k] = load_prior_n_data(S.d,'dID');   %is the loaded data in k-t space?

%% Gen transform operator
if f.ld.m; load(S.mask,'E'); %loads dUS, kt_ind, and kt_st from previous cycle
       Em  = xfm_NUFFT([E.Nd(1) E.Nd(1) 1  1],[],[],reshape(E.k,[], 1,2),'table',true);
else       
   [E, Em] = PERRIgenE(S.d, N.l, f.k, start_frame);
end
N.t        = E.Nt;
N.x        = prod(E.Nd);
N.c        = E.Nc;
N.k        = E.Kd(1);

%% Gen k-space data
d          = load_prior_n_data(S.d,'d');

if f.k; if start_frame > 0
            d_index     = vec( (start_frame:size(d,2)/N.t:size(d,2)) + (1:N.l)' );   % % All the lines needed for the recon
            d           = d(:,d_index,:);               % Resizing k
        end
        d = reshape(d(:,1:N.l*N.t,:),N.k*N.l,N.t,N.c); % Remove any leftover data that didn't fit into a whole sampling line, and reshape

else;   d = E * reshape(d,N.x,N.t);                          % Else, create k-space data for retrospectives by nufft transform
end


%% Adding noise?
if N.noise>0
    rngsum = sum(S.p) + strcmp(S.p(1),'c')*4 -  strcmp(S.p(1),'s')*12 - sum(S.p(regexp(S.p,'Tik'):regexp(S.p,'Tik')+7)); %ugly, but gives low res generated priors the same rng as the hybrid prior
    rng(rngsum);
    if strcmp(S.d,'32')&&~strcmp(S.p(1),'g')&&~strcmp(S.p(1),'c')&&~strcmp(S.p(1),'s')&&~strcmp(S.p(1),'t')&&~strcmp(S.p(1),'l')
        rng('shuffle');
    end

    d = d+randn(N.k*N.l,N.t,N.c) / sqrt(N.k*N.l*N.t*N.c*2) * norm(d(:)) * N.noise * 1 + randn(N.k*N.l,N.t,N.c) / sqrt(N.k*N.l*N.t*N.c*2) * norm(d(:)) * N.noise * 1i ;

rng('shuffle')
end

%% Variable shortcut
D          = E'*d;

%% Set recon parameters
N.tol     = 1E-15;
f.vpsf     = contains(f.init,'vpsf')||contains(f.init,'uvpsf')||contains(f.init,'-p')||f.vpsf;  % Whether we are doing a Partially separable functions recon in time
f.upsf     = contains(f.init,'upsf')||contains(f.init,'uvpsf')||f.upsf;                         % Whether we are doing a Partially separable functions recon in space

%% Gen priors
if f.ld.p; 
    try [Up, Vp] = load_prior_n_data(S.p,'p');
    catch
        if strcmp(S.p(1),'z')
           [Up, Vp] = PERRIgenPriors(d, E, N, rWin, S, Em, 1,f); %if it's a zero prior, just create it
        end
    end
else;      [Up, Vp] = PERRIgenPriors(d, E, N, rWin, S, Em, 1,f); % If priors are to be generated this cycle (and not loaded from a previous version) 
end
Up      = Up(:,1:N.r); 
Vp      = Vp(:,1:N.r);

if strcmp(S.p,'g')
try         [Up,Vp] = lsvdn(Up*Vp',N.r,'S'); %makes Up = Up*Sp
    if sum(isnan(Up(:)))>0||sum(isnan(Vp(:)))>0; Up(:) = 0;Vp(:) = 0; end %being ultra safe
catch;       Up(:)=0; Vp(:) = 0; %just make them zero if there's an error (likely due to them being zeros initially anyway)
end
end

%% Init U and V
[U,V,C] = init_UV(f.init,N.r,d,E,Em,Up,Vp,S);
V=fliplr(V);
if size(U,1)~=size(Up,1) || size(V,1)~=size(Vp,1) %Sanity check
    warning('Prior size isn''t the same size as initialised basis set')
    fprintf('Size Space: %d, Size Space Prior:%d\nSize Time: %d, Size Time Prior:%d\n',size(U,1),size(Up,1),size(V,1),size(Vp,1))
end

if f.upsf; U=Up; end
if f.vpsf; V=Vp; end
    
%% Masking 
mask = load_prior_n_data(S.d,'m');  % mask = repmat(mask,1,N.r); %Uncomment for backwards compatability
U    = U.*mask;
Up   = Up.*mask;

%% Full rank recon
if f.print; fprintf('PERRIcore:: Nl: %1.2f; Iters: %i; Data ID: %s; Prior ID: %s; lambda_U: %2.2e; lambda_V: %2.2e; lambda_s: %2.2e; Init Method: %s; Rank Method: %s\n', N.l, N.outer,S.d,S.p, lamu, lamv, lams, f.init, f.rmet); end
fprintf('flag v psf: %d, flag u psf: %d\n',f.vpsf,f.upsf)
[U, V, Ctemp] = PERRIoptm(d,D,E,U,V,[lamu lamv lams],Up,Vp,N,f,S);
%U = Window2Ddata(U);
C = [C;Ctemp];

%% Save output
save(S.uv,'U','V','C')
fprintf('Final Matrix saved\n')

end

%% Glossary
%
% C: Costs: [by column] cost function, u residual, v residual, internal u loops, internal v loops, true error
% d: nufft(Truth), aka the undersampled k-space data
% D: nufft_adj(d), or underlying truth in retrospective data
% E: Nufft structure
% f: Flags: .k, whether data is in k-space; .psf, whether we're doing a psf recon; rmet: , type of recon we're doing; .ld.p, if we're loading a previous prior; .ld.m; if we're loading a previous mask
% N: Structure containing all sizes .l:lines, .x:pixels, .t: frames, .r: rank, .inner: sub_problem iterations, .outer: switches between u and v, .tol = tolerence , .noise = noise level
% S: Strings to save/load: .uv: U, V, C; .mask: E; .temp: nothing; .p: Priors; .d: data; .init: initialisaiton type
