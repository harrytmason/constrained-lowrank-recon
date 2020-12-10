function [U,V,C] = init_UV(init,Nr,d,E,Em,Up,Vp,S)
%% init_UV.m
%
% An updated initialisation function for U and V
% Rarely meant to be accessed on it's own.
%
% Inputs:   init    - the string for initialisation
%           N.r     - the number of components
%           d       - the k-space data
%           E       - the nufft_xfm structure
%           Em      - the nufft_xfm structure of the temporal mean
%           Up      - The prior estimate of the spatial components (N.x * N.r)
%           Vp      - The prior estimate of the temporal components (N.t * N.r)
%           S       - a structure of strings
%
% Outputs: U (the spatial component) - N.x * N.r (N.x is number of pixels)
%          V (the temporal component) - N.t * N.r
%          C (the costs - usually empty, unless loading a previously run reconstruction)

C=[];
Nx = E.msize(1);
Nt = E.msize(2);

switch init
    
    case{'0','zeros'}
        U = zeros(Nx,Nr)+eps; V = zeros(Nt,Nr)+eps;
        
    case{'00'}
        U = zeros(Nx,Nr);     V = zeros(Nt,Nr);
        
    case{'1','ones'}
        U =  ones(Nx,Nr);     V = ones(Nt,Nr);% Never have both recons as one - always make sure one is starting below norm of dataset
        
    case{'01'}
        U =  zeros(Nx,Nr);     V = ones(Nt,Nr);% Never have both recons as one - always make sure one is starting below norm of dataset
        
    case{'1e6','mill'}
        U =  ones(Nx,Nr)*1e6; V = ones(Nt,Nr)*1e6;% Never have both recons as one - always make sure one is starting below norm of dataset
        
    case{'r','randn'}
        U = randn(Nx,Nr);     V = randn(Nt,Nr);
        
    case{'leg'} %initialise v with the legendre polynomial
        for i1 = 1:Nr
            p(i1,:) = legendreP(i1-1,linspace(-1,1,Nt));
        end
        V=orth(p');
        U=zeros(Nx,Nr);
        
    case{'or','orthrandn'}
        U = orth(randn(Nx,Nr));     V = orth(randn(Nt,Nr));    
         
    case{'truth'} %load truth (for debugging)
        [U,      V]      = load_prior_n_data(S.d,'l',Nr);
        U=U*0;
        
    case{'p','prior','uvpsf'}
        U = Up(:,1:Nr);       V = Vp(:,1:Nr); % [U, V] = load_prior_n_data(pID,'p');
        
    case{'m','mean'}
        U = Up(:,1:Nr);       V = Vp(:,1:Nr);
        U = mean(U(:)) + U*0; V = mean(V(:)) + V*0;%U = U/norm(U); V = V/norm(V);
        
    case{'t','temporal'}
        [U,V] = init_UV('0',Nr,d,E);
        [U(:,1), V(:,1)] = lsvdn(repmat(Em'*reshape(d,Em.dsize),1,Nt),1);
        
    case{'tnew','tlsqr'}
        [U,V] = init_UV('0',Nr,d,E);
        try
            [U(:,1), V(:,1)] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,Nt),1);
        catch
            [U,V,C] = init_UV('t',Nr,d,E,Em);
            fprintf('Iterative temporal mean failed (one rank), producing cheap gridded temporal mean instead\n')
        end
        U = U/(Nt^0.25); V = V/(Nt^0.25); %Normalising factor
        
    case{'told','tfull'} %old temporal mean
        [U,      V]      = lsvdn(repmat(Em'*reshape(d,Em.dsize),1,Nt),Nr);
        
    case{'tnewfull','tlsqrfull'}
        try
            [U, V] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,Nt),Nr);
            %U(:,2:16) = U(:,2:16) * 10^mean(-6-log10(sum(U(:,2:16,1))));
            %V(:,2:16) = V(:,2:16) * 10^mean(-6-log10(sum(V(:,2:16,1))));
        catch
            [U,V,C] = init_UV('told',Nr,d,E,Em);
            fprintf('Iterative temporal mean failed (full rank), producing cheap gridded temporal mean instead\n')
        end
       % U = U/(Nt^0.25); V = V/(Nt^0.25);
        
    case{'tr','trlsqr'}
        try
            [U, V] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,Nt),Nr);
            U(:,2:16) = 0;
            V(:,2:16) = orth(randn(Nt,Nr-1));
        catch
            [U,V,C] = init_UV('tnewfull',Nr,d,E,Em);
            fprintf('It temporal mean failed w/ randn V failed, producing Iterative temporal mean (full rank) instead\n')
        end
        U = U/(Nt^0.25); V = V/(Nt^0.25);
        
    case{'trpr','trprlsqr'}
        try
            [U, ~] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,Nt),Nr);
            U(:,2:16) = 0;
            V = Vp;
        catch
            [U,V,C] = init_UV('tnewfull',Nr,d,E,Em);
            fprintf('It temporal mean w/ Vprior failed, producing Iterative temporal mean (full rank) instead\n')
        end
        %U = U/(Nt^0.25); V = V/(Nt^0.25);
        
    case{'trr'}
        try
            [U, V] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,Nt),Nr);
            U(:,2:16) = orth(randn(Nx,Nr-1));
            V(:,2:16) = orth(randn(Nt,Nr-1));
        catch
            [U,V,C] = init_UV('tnewfull',Nr,d,E,Em);
            fprintf('It temporal mean failed w/ randn V failed, producing Iterative temporal mean (full rank) instead\n')
        end
        U = U/(Nt^0.25); V = V/(Nt^0.25);
    
    case{'trrwin','trrw'}
        try
            [U, V] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,Nt),Nr);
            U(:,2:16) = orth(randn(Nx,Nr-1));
            V(:,2:16) = orth(randn(Nt,Nr-1));
        catch
            [U,V,C] = init_UV('tnewfull',Nr,d,E,Em);
            fprintf('It temporal mean failed w/ randn V failed, producing Iterative temporal mean (full rank) instead\n')
        end
        U = U/(Nt^0.25); V = V/(Nt^0.25);
        U = Window2Ddata(U);    
        
    case{'trwin','trw'}
        try
            [U, V] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,Nt),Nr);
            U(:,2:16) = 0;
            V(:,2:16) = orth(randn(Nt,Nr-1));
        catch
            [U,V,C] = init_UV('tnewfull',Nr,d,E,Em);
            fprintf('It temporal mean failed w/ randn V failed, producing Iterative temporal mean (full rank) instead\n')
        end
        U = U/(Nt^0.25); V = V/(Nt^0.25);
        U = Window2Ddata(U);
        
    case{'trpre','preloaded'}
        try
            load('VNt300Nr0216.mat','Vor');
            [U, V] = lsvdn(repmat(Em.iter(reshape(d,Em.dsize), @pcg, 1e-6,1e2),1,Nt),Nr);
            U(:,2:16) = 0;
            V(1:Nt,2:Nr) = Vor(1:Nt,1:Nr-1);
        catch
            [U,V,C] = init_UV('tr',Nr,d,E,Em);
            fprintf('It temporal mean failed w/ loaded randn V failed, producing without pre-loaded V\n')
        end
        U = U/(Nt^0.25); V = V/(Nt^0.25);
        
    case{'lsqr','l'}
        [U, V] = lsvdn(E.iter(d, @pcg, 1e-6,1e2),Nr);
        
    case{'ltm','lsqrtemporal'}
        [U, V] = init_UV('lsqr',Nr,d,E);
        [U(:,1),V(:,1)] = init_UV('tnewfull',1,d,E,Em);
        
    case{'svd'}
        [U,      V]      = lsvdn(reshape(E'*d,Nx,Nt), Nr);
        
    case{'psf','vpsf'}
        %U = zeros(Nx,Nr)+eps;
        V = Vp;  %U starts as temporal mean, V doesn't change.
        [U,~] = init_UV('tr',Nr,d,E,Em);
        
    case{'upsf'}
        %U = zeros(Nx,Nr)+eps;
        U = Up;  %U starts as temporal mean, V doesn't change.
        [~,V] = init_UV('tr',Nr,d,E,Em);
        
    case{'prev','-1','-p','-psf'}
        if exist(strcat(S.uv(1:end-4),'Conv.mat'),'file')
                        fprintf('Converged version of file exists, stopping recon now\n')
            this_will_stop_the_program_running
        else
              try    load(S.uv,'U','V','C')
                            fprintf('Starting from previous version of file\n')
              catch; [U,V] = init_UV('p',Nr,d,E,Em,Up,Vp); %If file doesn't exist, start from priors
                     disp('No previous save file found')
              end
              if sum(V(:))==0 %don't start with identically zero prior in time.
                  [U,V] = init_UV('tr',Nr,d,E,Em);
              end
        end
        
   case{'new','n'} %Only want to run if no previous version exists
        if exist(S.uv,'file')
            this_will_prevent_the_program_running
        else
            [U,V] = init_UV('trpre',Nr,d,E,Em);
            disp('No previous save file found')
        end
        
   case{'newp','np'} %Only want to run if no previous version exists (initiate with priors
        if exist(S.uv,'file')
            this_will_prevent_the_program_running
        else
            [U,V] = init_UV('p',Nr,d,E,Em,Up,Vp);
            disp('No previous save file found')
        end
        
    case{'prevTemp','-2','-t'}
        try load(S.temp);
        catch; [U,V] = init_UV('tnewfull',Nr,d,E,Em); disp('No prior version found, starting from scratch');
        end  %Continue from partway through a simulation. If it had already completed, return. If no optimisation is saved, start from scratch.
        
    otherwise  
        if exist(init,'file')
            load(init,'U','V','Up','Vp');
            if exist('U','var') %default is to accept U and V values
                disp('Starting from prexisting reconstruction')
            elseif exist('Up','var')
                U = Up; V=Vp;
                disp('Starting from prexisting priors')
            else
                disp('No valid option selected, stopping recon')
                this_line_will_stop_the_program
            end
        else
        fprintf('No legitimate initial setting has been chosen. Check init_UV.m for valid options.\nOption selected was %s.\n',init)
        end
end