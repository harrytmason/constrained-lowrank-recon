function testCLR(Sd)

%%
% testCLR1.m is the first of three test functions, demonstrating the algorithm on
% a small dataset.
%
% This code was written for clarity, rather than for optimised behaviour or
% even to produce optimised results (which would take a much longer time to
% run).
%
%% Basics for whole function

rng('shuffle'); %ensure that wherever this is run, we have a unique randomisation seed
try 
addpath MC_utils/optim/transforms MC_utils/misc Misc irt Datasets
catch
    fprintf('Are you currently in the correct folder? You should be in the folder containing testCLR, or edit the line above this warning in the code.')
    return
end
setup; clear tmp irtdir;

%% Initialisation

% S is a structure storing all strings
%S.d     = '1';        % Input dataset 
S.d=Sd;
                      % '1' = test dataset (recommended) 
                      % '2' = retrospective datset A,
                      % '3' = 1st slice of prospective dataset, 
S.xt    = strcat('XT',S.d,'.mat');  % Output save string
fprintf('File save location: %s\n', S.xt)

% N is a structure storing integers. Default values are for speed (i.e. to just let the algorithm run, since full optimisations may run for over an hour)
N.l       = 5;      % Number of radial sampling lines (5 is default)
N.outer   = 100;     % Max number of alternations between subproblems (5 is default for speed purposes, we set to 1000 and let it converge)
N.inner   = 5;      % Number of iterations within each subproblem (2 is default for speed purposes, we used 10)
N.r       = 16;     % Desired rank (16 is default)
N.stol    = 1E-15;  % Subproblem tolerance (1e-15 is default)
N.tol     = 1E-2;   % Convergence criterion (1e-2 is default, we used 1e-3 for prospective data and 1e-5 for retrospective data)

% f is a structure storing flags
[~,f.k] = load_prior_n_data(S.d,'dID');   % Generate retrospective/prospective flag (could remove considering I'm pre-selecting the dataset)
f.psf   = 0; %Initialise as not doing psf reconstruction

%% Generate NUFFT transform operator E

[E, Em] = CLRgenE(S.d, N.l, f.k); % Em is the temporal mean nufft operator

 N.t        = E.Nt;        % Number of frames
 N.x        = prod(E.Nd);  % Nummber of voxels
 N.c        = E.Nc;        % Number of coils
 N.k        = E.Kd(1);     % Number of k-space points

 %% Generate k-space data (could remove considering I'm pre-selecting the dataset)

 d          = load_prior_n_data(S.d,'d');

if f.k; d = reshape(d(:,1:N.l*N.t,:),N.k*N.l,N.t,N.c);       % Remove any leftover data that didn't fit into a whole sampling frame, and reshape
else;   d = E * reshape(d,N.x,N.t);                          % Else, create k-space data for retrospectives by nufft transform
end

D       = E'*d; %This is the x-t transform of the undersampled data, with d now the k-t form of the data

%% Initialise X and T

[X, T] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,N.t),N.r); %Initialise the first component as the temporal mean
 X(:,2:N.r) = zeros(N.x,N.r-1);
 T(:,2:N.r) = orth(randn(N.t,N.r-1));
 
%% Masking 
mask = load_prior_n_data(S.d,'m');  % masking out the values which are zero in sensitivity maps
X    = X.*mask;


%% Reconstructions
%%% 
%%%% 
%%%%%

%% k-t FASTER reconstruction
Xp = X.*0;
Tp = T.*0;
lam = [0 0 0];      % [lamx, lamt, lams] are all zero for k-t FASTER

fprintf('k-t FASTER:: Nl: %1.2f; Iters: %i; Data ID: %s; lambda_X: %2.2e; lambda_T: %2.2e; lambda_s: %2.2e\n', N.l, N.outer,S.d, lam(1), lam(2), lam(3));

[Xktf, Tktf, Cktf] = CLRoptm(d,D,E,X,T,lam,Xp,Tp,N,f,S); 

%% Tikhonov-constrained reconstruction
Xp = X.*0;
Tp = T.*0;
lam = [1e-7 1e-7 0];      % [~,~,lams] is zero for Tikhonov, [lamx, lamt, ~] are non-zero

fprintf('Tikhonov:: Nl: %1.2f; Iters: %i; Data ID: %s; lambda_X: %2.2e; lambda_T: %2.2e; lambda_s: %2.2e\n', N.l, N.outer,S.d, lam(1), lam(2), lam(3));

[Xtik, Ttik, Ctik] = CLRoptm(d,D,E,X,T,lam,Xp,Tp,N,f,S);

%% Smoothness-constrained reconstruction
Xp = X.*0;
Tp = T.*0;
lam = [0 0 1e-6];      % [lamx, lamt, ~] are zero for Temporal Smoothness, [~,~,lams] is non-zero

fprintf('Smoothness:: Nl: %1.2f; Iters: %i; Data ID: %s; lambda_X: %2.2e; lambda_T: %2.2e; lambda_s: %2.2e\n', N.l, N.outer,S.d, lam(1), lam(2), lam(3));

[Xsmo, Tsmo, Csmo] = CLRoptm(d,D,E,X,T,lam,Xp,Tp,N,f,S);

%% LRP-constrained reconstruction
%generate priors first, then generate LRP-constrained recon

%initialise for priors
lam = [0 0 0];
N.rad  = 0.3; % Window of low res prior radius. We used a radius based on acceleration factor.
tukey_parameter = 0.4;
kWin  = Window1D(N.k,N.rad,tukey_parameter);
dWin = d.*repmat(kWin,N.l,N.t,N.c); %Replicating the k-space window across every line, frame and coil, then multiplying by data
DWin = E'*dWin;

%generate priors
fprintf('Priors:: Nl: %1.2f; Iters: %i; Data ID: %s; lambda_X: %2.2e; lambda_T: %2.2e; lambda_s: %2.2e\n', N.l, N.outer,S.d, lam(1), lam(2), lam(3));
[Xp, Tp] = CLRoptm(dWin,DWin,E,X,T,lam,Xp*0,Tp*0,N,f,S);  
[Xp,Tp] = lsvdn(Xp*Tp',N.r,'S');
 Xp = Xp.*mask; %ensure output is masked for safety

%then generate LRP-constrained output
lam = [1e-7 1e-7 0];      % [lamx, lamt] are non-zero for LRP, [~,~,lams] is zero
fprintf('LRP:: Nl: %1.2f; Iters: %i; Data ID: %s; lambda_X: %2.2e; lambda_T: %2.2e; lambda_s: %2.2e\n', N.l, N.outer,S.d, lam(1), lam(2), lam(3));

[Xlrp, Tlrp, Clrp] = CLRoptm(d,D,E,X,T,lam,Xp,Tp,N,f,S);

%% k-t PSF reconstruction
%Uses same priors as LRP recon

lam = [0 inf 0];      % [lamx, lamt, ~] are non-zero for LRP, [~,~,lams] is zero
f.psf     = isinf(lam(2));  % set a psf flag

fprintf('k-t PSF:: Nl: %1.2f; Iters: %i; Data ID: %s; lambda_X: %2.2e; lambda_T: %2.2e; lambda_s: %2.2e\n', N.l, N.outer,S.d, lam(1), lam(2), lam(3));

[Xpsf, Tpsf, Cpsf] = CLRoptm(d,D,E,X,T,lam,Xp,Tp,N,f,S);

%% Truth generation
[Xtru,Ttru] = load_prior_n_data(S.d,'t');

%% Save outputs
save(S.xt,'Xktf','Tktf','Xtik','Ttik','Xsmo','Tsmo','Xlrp','Tlrp','Xpsf','Tpsf');

%% CCS Comparison

if ~f.k %Doing a CCS comparison with prospective data isn't as useful for interpretability since the 'ground truth' has varying temporal resolution (we relied more on z-stat analysis, which would require the user to download fsl)
[~,Xccs(1)] = cssa(Xtru,Xktf);
[~,Xccs(2)] = cssa(Xtru,Xtik);
[~,Xccs(3)] = cssa(Xtru,Xsmo);
[~,Xccs(4)] = cssa(Xtru,Xlrp);
[~,Xccs(5)] = cssa(Xtru,Xpsf);

Xccs = Xccs./N.r;

[~,Tccs(1)] = cssa(Ttru,Tktf);
[~,Tccs(2)] = cssa(Ttru,Ttik);
[~,Tccs(3)] = cssa(Ttru,Tsmo);
[~,Tccs(4)] = cssa(Ttru,Tlrp);
[~,Tccs(5)] = cssa(Ttru,Tpsf);

Tccs = Tccs./N.r;
end

%% Plotting

if f.k, Nrows = 2;
else;   Nrows = 3;
end

% Brain plots
figure
subplot(Nrows,6,1:2)
D=Xktf*Tktf';
show(rot90(reshape(abs(D(:,1)),sqrt(N.x),sqrt(N.x)),3),[])
title('k-t FASTER')

subplot(Nrows,6,3:4)
D=Xtik*Ttik';
show(rot90(reshape(abs(D(:,1)),sqrt(N.x),sqrt(N.x)),3),[])
title('Tikhonov')

subplot(Nrows,6,5:6)
D=Xsmo*Tsmo';
show(rot90(reshape(abs(D(:,1)),sqrt(N.x),sqrt(N.x)),3),[])
title('Smoothness')

subplot(Nrows,6,7:8)
D=Xlrp*Tlrp';
show(rot90(reshape(abs(D(:,1)),sqrt(N.x),sqrt(N.x)),3),[])
title('LRP')

subplot(Nrows,6,9:10)
D=Xpsf*Tpsf';
show(rot90(reshape(abs(D(:,1)),sqrt(N.x),sqrt(N.x)),3),[])
title('k-t PSF')

subplot(Nrows,6,11:12)
D=Xtru*Ttru';
show(rot90(reshape(abs(D(:,1)),sqrt(N.x),sqrt(N.x)),3),[])
title('Truth')


%CCS plots

if ~f.k % No CCS comparison with prospective data

subplot(Nrows,6,13:15)
bar(Xccs','k');
axis([0.5 5.5 0 1])
set(gca,'XTickLabel',{'k-t FASTER','Tikhonov','Smoothness','LRP','k-t PSF'},'XTickLabelRotation',45)
title('X CCS')

subplot(Nrows,6,16:18)
bar(Tccs','k');
axis([0.5 5.5 0 1])
set(gca,'XTickLabel',{'k-t FASTER','Tikhonov','Smoothness','LRP','k-t PSF'},'XTickLabelRotation',45)
title('T CCS')
set(gca,'YAxisLocation','right')
end





