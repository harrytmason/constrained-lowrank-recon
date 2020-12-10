function [U, V, C] = PERRIoptm(d, D, E, U,V, lam,Up,Vp, N, f,S, Win)

%% The optimisation loop for k-t PERRI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function carries out alternating optimisation of the U (Spatial) and
% V (Temporal) subproblems, with various styles of alternation.
%
% Low rank conditions are enforced, as is regularisation with respect to a
% prior.
%
% Used in prior generation and the in the final reconstruction. Capable of
% k-t FASTER, k-t PSF and k-t PERRI reconstructions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:  D:   nufft_adj(k-space data)                                  (Nk , Nt)
%          E:   nufft structure + more                                   (xfm_nufft structure)
%          U:   current guess of u, the spatial subspace                 (Nx , Nr)
%          V:   current guess of v, the temporal subspace                (Nt , Nr)
%          lam: the weighting on the priors                              [lambda_u, lambda_v]
%          Up:  the prior of u                                           (Nx , Nr)
%          Vp:  the prior of v                                           (Nt , Nr)
%          N.tol: tolerance of problem                                   scalar
%          N.r:   number of ranks in u and v                             scalar
%          N.in:  number of iterations of subproblem                     scalar
%          N.out: number of switches between u and v per rank            scalar
%          T:     underlying Truth (not necessary)                       (Nx , Nt)
%          f.rmet: optimisation method (full, irpf, rbyr)                string
%          f.upsf: whether we're doing a u psf recon (so not updating u) boolean
%          f.vpsf: whether we're doing a v psf recon (so not updating v) boolean
%          f.print: whether to print out information as you go           boolean
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs: U:   final guess of u, the spatial subspace                  (Nx , Nr)
%          V:   final guess of v, the temporal subspace                 (Nt , Nr)
%          C:   Seven columns of optimisation information               (Nr * Nout, 7)
%            C(n,1): The nth residual of the U subproblem
%            C(n,2): The number of times the U subproblem iterated in the nth pass
%            C(n,3): The nth residual of the V subproblem
%            C(n,4): The number of times the V subproblem iterated in the nth pass
%            C(n,-19): The combined u v cost function after n iterations
%            C(n,6): The error with respect to a ground truth after n iterations
%            C(n,7): The time taken after n iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Required Functions: PERRIcost.m, the underlying cost function
%                     PERRIusub.m, the expression of the u (spatial) subproblem
%                     PERRIvsub.m, the expression of the v (temporal) subproblem
%                     PERRItrueErr.m, the true state of the solution
%                     Dr. Mark Chiew's xfm_nufft structures (represented by E in the code)
%                     The irt toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by Harry Mason, M.Eng, University of Oxford
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALISATION

if nargin < 6;      lamu   = 0;
                    lamv   = 0;
                    lams   = 0;
else;               lamu   = lam(1);
                    if numel(lam)==1
                        lamv = lam(end);
                        lams = 0;
                    elseif numel(lam)>1
                        lamv = lam(2);
                        if numel(lam)>2
                            lams = lam(3);
                        else
                            lams = 0;
                        end
                    end
end
if nargin < 7;      Up     = zeros(size(U));    end
if nargin < 8;      Vp     = zeros(size(V));    end
if nargin < 9;      utol   = 1e-15;
                    vtol   = 1e-15;
                    Nr     = size(U,2);
                    Nin    = 1;
                    Nout   = 1;
                    Nx     = size(U,1);
                    Nt     = size(V,1);
                    Nl     = 1;
else;               utol   = N.tol;     % Accessing from a structure takes more time, so we extract all the values here.
                    vtol   = N.tol;
                    Nr     = N.r;
                    Nin    = N.inner;
                    Nout   = N.outer;
                    Nx     = N.x;
                    Nt     = N.t;
                    Nl     = N.l;
end
if nargin < 10;     frmet  = 'rbyr';
                    fupsf  = false;
                    fvpsf  = false;
                    fprint = true;
                    fk     = false;
else;               frmet  = f.rmet;
                    fupsf  = f.upsf;
                    fvpsf  = f.vpsf;
                    fprint = f.print;
                    fk     = f.k;

end
if nargin < 11;     Suv = 'UVNo0000.mat'; fsave = 0; %%Debug option to save output for intermediate iterations
                    Sd  = '';
else; fsave = 1;       
    if isstruct(S); Suv = S.uv;
                    Sd  = S.d;
    else
                    Suv = S;
                    Sd = '';
    end  
end
if nargin < 12;     Win = ones(Nx,Nr); end

svInd = [1 2 5 10:10:90 100 Nout]; %Which iterations to save on

% Ensure priors have correct weighting
if strcmp(S.p(1),'g')
try                                 %   This bit is also in PERRIcore.m, but some redundancy for safety isn't always bad
    [Up,Vp] = lsvdn(Up*Vp',Nr,'S'); % makes Up = Up*sqrt(Sp)
    if sum(isnan(Up(:)))>0||sum(isnan(Vp(:)))>0; Up(:) = 0;Vp(:) = 0; end
catch
    Up(:) = 0;Vp(:) = 0;            % just make them zero if there's an error (likely due to them being zeros initially anyway)
end
end

mask = load_prior_n_data(Sd,'m'); % mask = repmat(mask,1,Nr); %Uncomment for backwards compatability


%% PRINT HEADER

if fprint; fprintf('%3s\t%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\n','Rank','Cost Func U','Cost Func V','Cost Func','CF Grad','ItU','ItV','Time','It Loop'); end
tic
ftol3 = 0;
ftol4 = 0;
ftol55 = 0;

%% FLAG CHECKS
if fupsf && fvpsf %% it's going to be the same each time, so I don't care about doing it Nr*Nout times...
    C(:,5) = PERRIcost(d,E,U,V,lamu,lamv,Up,Vp);
    C(:,6) = 0;
    C(:,7) = toc;
    if fsave %For convenience, more than anything
            for ii = svInd
                svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                save(svStr,'U','V','C')
            end
    end
    if fprint;     fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(end,1),C(end,2),C(end,5),C(end,6),C(end,3),C(end,4),round(C(end,7)),1); end
else
    
    frewe = (sum(lamu+lamv)==0) && ~fupsf && ~ fvpsf; %If we have any non-zero lambda, we don't need to reweight (psf flags as security for any infinities)
    fdebug = 0; %While debugging, there's a few metrics I want to check. 

%% RECONSTRUCTION
switch frmet

    case{'full'}
        %% Full rank recon
        C = zeros(Nout,7);
         
      %% PSF loop (no subproblem switching needed)
      if      fvpsf || fupsf 
            
            % Reconstruction
            if fupsf;   [V, C(:,2), C(:,4)] = PERRIvsub(D,E,U,V,lamv,Vp,vtol,Nin,Nt,Nr);   end %lamu = inf, so just do the v subproblem
            if fvpsf;   [U, C(:,1), C(:,3)] = PERRIusub(D,E,U,V,lamu,Up,utol,Nin*10,Nx,Nr); end %lamv = inf, so just do the u subproblem
            
            % Calculate cost function
                        C(:,5) = PERRIcost(d,E,U,V,lamu,lamv,Up,Vp);
                        C(:,6) = 0;
                        C(:,7) = toc;
            
            %Save output
            if fsave %For convenience, more than anything
                    for ii = svInd
                        svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                        save(svStr,'U','V','C')
                        save(Suv,'U','V','C') %Also save the most up to date version
                        svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
                        save(svStr,'U','V','C')
                    end
            end
   
            ii = Nout;
      else % The non-PSF loop (requires iterations
        

            
        % Load the truth for canonocal correlation checks
        if fdebug; [Ut,Vt]  = load_prior_n_data(strcat(Sd,'w'),'l',16);
                    Ut      = Ut.*mask;
                    if fk; Ut = U; Vt = V;
                    end
                    c       = zeros(Nout,4);
        end

        U = U.*mask;
        frewe=0;
        nUV=zeros(1000,2);
        
        %% Iterations
        for ii = 1:Nout
             
            % Canoncorr debug
            if fdebug;      [uccfull ucc(ii)] = cssa(U,Ut,1);   
                            [vccfull vcc(ii)] = cssa(V,Vt);     
            end
             
            
            % Spatial subproblem
                            [U, C(ii,1), C(ii,3)] = PERRIusub(D,E,U,V,lamu,Up,utol,Nin*10,Nx,Nr);   %put back in *10  
            if frewe ;      [U,V] = lsvdn(U*V',Nr,"V");      end %Reweight if kt FASTER
            
            
            % Temporal subproblem
                            [V, C(ii,2), C(ii,4)] = PERRIvsub(D,E,U,V,lamv,Vp,vtol,50,   Nt,Nr); 
            if frewe ;      [U,V] = lsvdn(U*V',Nr,"U");                                              end
                        
            % Calculate cost function (unbelievably no flags here)
                            C(ii,[5 7]) = [PERRIcost(d,E,U,V,lamu,lamv,Up,Vp)     toc];
            if ii>1;        C(ii,6) = (C(ii,5) - C(ii-1,5))/C(ii,5);                                 end
            
            % Constituent parts of cost function for debugging
            if fdebug;      c(ii,1) = norm(vec(d-E*(U*V')), "fro" )^2;                              
                            c(ii,2) = lamu^2*norm(U - Up, "fro" )^2;                                
                            c(ii,3) = lamv^2*norm(V-Vp,"fro")^2;                                    
                            c(ii,4) = c(ii,1)+c(ii,2)+c(ii,3);   
            if ii>1;        c(ii,5) = (C(ii,5) - C(ii-1,5))/norm(d(:));                              end
            end

            if ii==10;
                1
            end

            % Print out current status
            if mod(ii,1)==0; if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end; end
            
            nUV(ii,:) = [norm(U(:)) norm(V(:))];
            
            %[~,em] = componentThreshold(V,Sd,Nl,0.4,0)
            
            save(Suv,'U','V','C'); %Always save current version
            % Save intermediate iterations
            if sum(ii == svInd) && fsave
                svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                save(svStr,'U','V','C')
            end
            

            if abs(C(ii,6))<1e-5 && ii>=10 %At least five iterations, and a relative change in gradient of less than 1e-5
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                for jj = (1+ii):1:Nout
                    if sum(jj == svInd) && fsave
                        svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(jj))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                        save(svStr,'U','V','C')
                    end
                end
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
                U = Window2Ddata(U);
                save(svStr,'U','V','C')
                fprintf('Reconstruction reached a tolerance with a gradient of %1.3e at iteration %d\n',C(ii,6),ii)
                return
            elseif abs(C(ii,6))<5e-5 && ~ftol55 && ii>1 %At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.5-e5loT','once')); %Include a "Tol5e-5" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol55 = 1;
            elseif abs(C(ii,6))<1e-4 && ~ftol4 && ii>1%At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.4-e1loT','once')); %Include a "Tol1e-4" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol4 = 1;
            elseif abs(C(ii,6))<1e-3 && ~ftol3 && ii>1%At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.3-e1loT','once')); %Include a "Tol1e-3" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol3 = 1;
            end
                
        end
        
      end
        % if psf, pour everything into Nin.
        if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end

   
    case{'nful'}
        %% Full rank recon with null priors
        C = zeros(Nout,7);
        [good_vec] = componentThreshold(Vp,Sd,Nl);
        Up(:,~good_vec)=0;  Vp(:,~good_vec)=0;
        Up = null(Up');     Vp = null(Vp');
        Upn = (Up*Up');      Vpn = (Vp*Vp');
        
        	U = U.*mask;
        %% PSF loop (no subproblem switching needed)
        if      fvpsf || fupsf 
            
            % Reconstruction
            if fupsf;   [V, C(:,2), C(:,4)] = PERRIvsubnull(D,E,U,V,lamv,Vpn,vtol,Nin,Nt,Nr);   end %lamu = inf, so just do the v subproblem
            if fvpsf;   [U, C(:,1), C(:,3)] = PERRIusubnull(D,E,U,V,lamu,Upn,utol,Nin*10,Nx,Nr); end %lamv = inf, so just do the u subproblem
            
            % Calculate cost function
                        C(:,5) = PERRIcostnull(d,E,U,V,lamu,lamv,Up,Vp);
                        C(:,6) = 0;
                        C(:,7) = toc;
            
            %Save output
            if fsave %For convenience, more than anything
                    for ii = svInd
                        svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                        save(svStr,'U','V','C')
                        save(Suv,'U','V','C') %Also save the most up to date version
                    end
            end
   
            ii = Nout;
        else % The non-PSF loop (requires iterations
            
        % Load the truth for canonocal correlation checks
        if fdebug; [Ut,Vt]  = load_prior_n_data(strcat(Sd,'w'),'l',16); 
                    c       = zeros(Nout,4);
        end
        
        %% Iterations
        for ii = 1:Nout
            
            % Canoncorr debug
            if fdebug;      [~,ucc(ii)] = cssa(U,Ut,1);           
                            [~,vcc(ii)] = cssa(V,Vt);                             
            end
            
            % Spatial subproblem
                            [U, C(ii,1), C(ii,3)] = PERRIusubnull(D,E,U,V,lamu,Upn,utol,Nin*10,Nx,Nr);     
            if frewe ;      [U,V] = lsvdn(U*V',Nr,'V');                                              end %Reweight if kt FASTER
                     
            % Temporal subproblem
                            [V, C(ii,2), C(ii,4)] = PERRIvsubnull(D,E,U,V,lamv,Vpn,vtol,Nin,   Nt,Nr); 
            if frewe ;      [U,V] = lsvdn(U*V',Nr,'U');                                              end
                        
            % Calculate cost function (unbelievably no flags here)
                            C(ii,[5 7]) = [PERRIcostnull(d,E,U,V,lamu,lamv,Up,Vp)     toc];
            if ii>1;        C(ii,6) = (C(ii,5) - C(ii-1,5))/C(ii,5);                                 end
            
            % Constituent parts of cost function for debugging
            if fdebug;      c(ii,1) = norm(vec(d-E*(U*V')), 'fro' )^2;                              
                            c(ii,2) = lamu^2*norm(Up'*U, 'fro' )^2;                                
                            c(ii,3) = lamv^2*norm(Vp'*V,'fro')^2;                                    
                            c(ii,4) = c(ii,1)+c(ii,2)+c(ii,3);   
            if ii>1;        c(ii,5) = (C(ii,5) - C(ii-1,5))/norm(d(:));                              end
            end

            % Print out current status
            if mod(ii,1)==0; if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end; end
        
            save(Suv,'U','V','C'); %Always save current version
            % Save intermediate iterations
            if sum(ii == svInd) && fsave
                svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                save(svStr,'U','V','C')
            end
            

            if abs(C(ii,6))<1e-5 && ii>=10 %At least two iterations, and a relative change in gradient of less than 1e-5
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                for jj = (1+ii):1:Nout
                    if sum(jj == svInd) && fsave
                        svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(jj))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                        save(svStr,'U','V','C')
                    end
                end
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
                U = Window2Ddata(U);
                save(svStr,'U','V','C')
                fprintf('Reconstruction reached a tolerance with a gradient of %1.3e at iteration %d\n',C(ii,6),ii)
                return
            elseif abs(C(ii,6))<5e-5 && ~ftol55 && ii>1 %At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.5-e5loT','once')); %Include a "Converged" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol55 = 1;
            elseif abs(C(ii,6))<1e-4 && ~ftol4 && ii>1 %At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.4-e1loT','once')); %Include a "Converged" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol4 = 1;
            elseif abs(C(ii,6))<1e-3 && ~ftol3 && ii>1 %At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.3-e1loT','once')); %Include a "Converged" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol3 = 1;
            end
                
        end
        
        end
        % if psf, pour everything into Nin.
        if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end
    case{'nall'}
        %% Full rank recon with null priors
        C = zeros(Nout,7);
        %[good_vec] = componentThreshold(Vp,Sd,Nl);
        %Up(:,~good_vec)=0;  Vp(:,~good_vec)=0;
        Up = null(Up');     Vp = null(Vp');
        Upn = (Up*Up');      Vpn = (Vp*Vp');
        
        	U = U.*mask;
        %% PSF loop (no subproblem switching needed)
        if      fvpsf || fupsf 
            
            % Reconstruction
            if fupsf;   [V, C(:,2), C(:,4)] = PERRIvsubnull(D,E,U,V,lamv,Vpn,vtol,Nin,Nt,Nr);   end %lamu = inf, so just do the v subproblem
            if fvpsf;   [U, C(:,1), C(:,3)] = PERRIusubnull(D,E,U,V,lamu,Upn,utol,Nin*10,Nx,Nr); end %lamv = inf, so just do the u subproblem
            
            % Calculate cost function
                        C(:,5) = PERRIcostnull(d,E,U,V,lamu,lamv,Up,Vp);
                        C(:,6) = 0;
                        C(:,7) = toc;
            
            %Save output
            if fsave %For convenience, more than anything
                    for ii = svInd
                        svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                        save(svStr,'U','V','C')
                        save(Suv,'U','V','C') %Also save the most up to date version
                    end
            end
   
            ii = Nout;
        else % The non-PSF loop (requires iterations
            
        % Load the truth for canonocal correlation checks
        if fdebug; [Ut,Vt]  = load_prior_n_data(strcat(Sd,'w'),'l',16); 
                    c       = zeros(Nout,4);
        end
        
        %% Iterations
        for ii = 1:Nout
            
            % Canoncorr debug
            if fdebug;      [~,ucc(ii)] = cssa(U,Ut,1);           
                            [~,vcc(ii)] = cssa(V,Vt);                             
            end
            
            % Spatial subproblem
                            [U, C(ii,1), C(ii,3)] = PERRIusubnull(D,E,U,V,lamu,Upn,utol,Nin*10,Nx,Nr);     
            if frewe ;      [U,V] = lsvdn(U*V',Nr,'V');                                              end %Reweight if kt FASTER
                     
            % Temporal subproblem
                            [V, C(ii,2), C(ii,4)] = PERRIvsubnull(D,E,U,V,lamv,Vpn,vtol,Nin,   Nt,Nr); 
            if frewe ;      [U,V] = lsvdn(U*V',Nr,'U');                                              end
                        
            % Calculate cost function (unbelievably no flags here)
                            C(ii,[5 7]) = [PERRIcostnull(d,E,U,V,lamu,lamv,Up,Vp)     toc];
            if ii>1;        C(ii,6) = (C(ii,5) - C(ii-1,5))/C(ii,5);                                 end
            
            % Constituent parts of cost function for debugging
            if fdebug;      c(ii,1) = norm(vec(d-E*(U*V')), 'fro' )^2;                              
                            c(ii,2) = lamu^2*norm(Up'*U, 'fro' )^2;                                
                            c(ii,3) = lamv^2*norm(Vp'*V,'fro')^2;                                    
                            c(ii,4) = c(ii,1)+c(ii,2)+c(ii,3);   
            if ii>1;        c(ii,5) = (C(ii,5) - C(ii-1,5))/norm(d(:));                              end
            end

            % Print out current status
            if mod(ii,1)==0; if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end; end
        
            save(Suv,'U','V','C'); %Always save current version
            % Save intermediate iterations
            if sum(ii == svInd) && fsave
                svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                save(svStr,'U','V','C')
            end
            

            if abs(C(ii,6))<1e-5 && ii>=10 %At least two iterations, and a relative change in gradient of less than 1e-5
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                for jj = (1+ii):1:Nout
                    if sum(jj == svInd) && fsave
                        svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(jj))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                        save(svStr,'U','V','C')
                    end
                end
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
                U = Window2Ddata(U);
                save(svStr,'U','V','C')
                fprintf('Reconstruction reached a tolerance with a gradient of %1.3e at iteration %d\n',C(ii,6),ii)
                return
            elseif abs(C(ii,6))<5e-5 && ~ftol55 && ii>1 %At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.5-e5loT','once')); %Include a "Converged" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol55 = 1;
            elseif abs(C(ii,6))<1e-4 && ~ftol4 && ii>1 %At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.4-e1loT','once')); %Include a "Converged" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol4 = 1;
            elseif abs(C(ii,6))<1e-3 && ~ftol3 && ii>1 %At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.3-e1loT','once')); %Include a "Converged" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol3 = 1;
            end
                
        end
        
        end
        % if psf, pour everything into Nin.
        if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end

   
    case{'irpf'}
        %% Incremented Rank PowerFactorization recon
        
        C = zeros(Nout*Nr,7);
        U = U.*mask;
         %% PSF loop (no subproblem switching needed)
      if      fvpsf || fupsf 
            
          for r=1:Nr
            hh = Nout*(r-1)+1;
            ii = Nout*r;
            % Reconstruction
            if fupsf;   [V(:,1:r), C(:,2), C(:,4)] = PERRIvsub(D,E,U(:,1:r),V(:,1:r),lamv,Vp(:,1:r),vtol,Nin,Nt,r);   end %lamu = inf, so just do the v subproblem
            if fvpsf;   [U(:,1:r), C(:,1), C(:,3)] = PERRIusub(D,E,U(:,1:r),V(:,1:r),lamu,Up(:,1:r),utol,Nin*10,Nx,r); end %lamv = inf, so just do the u subproblem
            
            % Calculate cost function
                        C(hh:ii,5) = PERRIcost(d,E,U(:,1:r),V(:,1:r),lamu,lamv,Up(:,1:r),Vp(:,1:r));
                        C(hh:ii,6) = 0;
                        C(hh:ii,7) = toc;
            
            %Save output
             if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',r,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end
                        
                 svStr = fliplr(regexprep(fliplr(Suv),'\d{3}rN',fliplr(sprintf('Nr%03i',round(r))),'once')); %Flipping is to ensure only last instance of NrXX is matched in the string
             save(svStr,'U','V','C')            
             save(Suv,'U','V','C') %Also save the most up to date version
               
          end
            
            svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
            save(svStr,'U','V','C')
                   
      else % The non-PSF loop (requires iterations
        
        
          for r = 1:Nr
              hh = Nout*(r-1)+1;
              
              for ii = hh:r*Nout
                  
                  [U(:,1:r), C(ii,1), C(ii,3)] = PERRIusub(D,E,U(:,1:r),V(:,1:r),lamu,Up(:,1:r),utol,Nin*10,Nx,r);
                  % if frewe ;      [U(:,1:r),V(:,1:r)] = lsvdn(U(:,1:r)*V(:,1:r)',r,"V");      end %Reweight if kt FASTER
                  
                  [V(:,1:r), C(ii,2), C(ii,4)] = PERRIvsub(D,E,U(:,1:r),V(:,1:r),lamv,Vp(:,1:r),vtol,Nin,Nt,r);
                  % if frewe ;      [U(:,1:r),V(:,1:r)] = lsvdn(U(:,1:r)*V(:,1:r)',r,"V");      end %Reweight if kt FASTER
                  
                  C(ii,[5 7]) = [PERRIcost(d,E,U(:,1:r),V(:,1:r),lamu,lamv,Up(:,1:r),Vp(:,1:r))     toc];
                  if ii>hh
                      C(ii,6) = (C(ii,5) - C(ii-1,5))/C(ii,5);
                  end
                  
                  svStr = fliplr(regexprep(fliplr(Suv),'\d{2}rN',fliplr(sprintf('Nr%02',round(r))),'once')); %Flipping is to ensure only last instance of NrXX is matched in the string
                    save(svStr,'U','V','C')
                  save(Suv,'U','V','C') %Also save the most up to date version
              if fprint;     fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',r,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii-hh+1); end
              
              
              if abs(C(ii,6))<1e-5 && ii>=(10+hh) %At least five iterations, and a relative change in gradient of less than 1e-5
                  C(ii+1:end,:) = repmat(C(ii,:),Nout*Nr-ii,1);
                  
                  
                  svStr = fliplr(regexprep(fliplr(svStr),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
                  U = Window2Ddata(U);
                  save(svStr,'U','V','C')
                  
                  fprintf('Reconstruction reached a tolerance with a gradient of %1.3e at iteration for rank %i %d\n',C(ii,6),ii-hh,r)
                  break
              elseif abs(C(ii,6))<5e-5 && ~ftol55 && ii>(1+hh) %At least ten iterations, and a relative change in gradient of less than 1e-5
                  svStr = fliplr(regexprep(fliplr(svStr),'tam.','tam.5-e5loT','once')); %Include a "Tol5e-5" save at the end
                  C(ii+1:end,:) = repmat(C(ii,:),(Nout*Nr)-ii,1);
                  save(svStr,'U','V','C')
                  ftol55 = 1;
              elseif abs(C(ii,6))<1e-4 && ~ftol4 && ii>(1+hh)%At least ten iterations, and a relative change in gradient of less than 1e-5
                  svStr = fliplr(regexprep(fliplr(svStr),'tam.','tam.4-e1loT','once')); %Include a "Tol1e-4" save at the end
                  C(ii+1:end,:) = repmat(C(ii,:),(Nout*Nr)-ii,1);
                  save(svStr,'U','V','C')
                  ftol4 = 1;
              elseif abs(C(ii,6))<1e-3 && ~ftol3 && ii>(1+hh)%At least ten iterations, and a relative change in gradient of less than 1e-5
                  svStr = fliplr(regexprep(fliplr(svStr),'tam.','tam.3-e1loT','once')); %Include a "Tol1e-3" save at the end
                  C(ii+1:end,:) = repmat(C(ii,:),(Nout*Nr)-ii,1);
                  save(svStr,'U','V','C')
                  ftol3 = 1;
              end
              %early break conditions
              end
                 svStr = fliplr(regexprep(fliplr(Suv),'\d{3}rN',fliplr(sprintf('Nr%03i',round(r))),'once')); %Flipping is to ensure only last instance of NrXX is matched in the string
              save(svStr,'U','V','C')
              save(Suv,'U','V','C') %Also save the most up to date version
              
              
          end

        
        
      end
    case{'rbyr'}
        %% Rank-by-rank recon
        C = zeros(Nout*Nr,7);
        tmp= 0;
        U = U.*mask;
         if      fvpsf || fupsf 
            
          for r=1:Nr
            hh = Nout*(r-1)+1;
            ii = Nout*r;
            % Reconstruction
            if fupsf;   [V(:,r), C(:,2), C(:,4)] = PERRIvsub(D-tmp,E,U(:,r),V(:,r),lamv,Vp(:,r),vtol,Nin,Nt,r);   end %lamu = inf, so just do the v subproblem
            if fvpsf;   [U(:,r), C(:,1), C(:,3)] = PERRIusub(D-tmp,E,U(:,r),V(:,r),lamu,Up(:,r),utol,Nin*10,Nx,r); end %lamv = inf, so just do the u subproblem
            
            
            tmp= E.mtimes2(U(:,1:r)*V(:,1:r)');
            
            % Calculate cost function
                        C(hh:ii,5) = PERRIcost(d,E,U(:,1:r),V(:,1:r),lamu,lamv,Up(:,1:r),Vp(:,1:r));
                        C(hh:ii,6) = 0;
                        C(hh:ii,7) = toc;
            
            %Save output
             if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',r,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end
                        
                 svStr = fliplr(regexprep(fliplr(Suv),'\d{3}rN',fliplr(sprintf('Nr%03i',round(r))),'once')); %Flipping is to ensure only last instance of NrXX is matched in the string
             save(svStr,'U','V','C')            
             save(Suv,'U','V','C') %Also save the most up to date version
               
          end
            
            svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
            save(svStr,'U','V','C')
                   
         else % The non-PSF loop (requires iterations
             
             for r = 1:Nr
                 hh = (r-1)*Nout+1;
                 for ii = hh:r*Nout
                     
                     
                     [U(:,r), C(ii,1), C(ii,3)] = PERRIusub(D-tmp,E,U(:,r),V(:,r),lamu,Up(:,r),utol,Nin,Nx,1); 
    %if frewe ;      [U(:,r),V(:,r)] = lsvdn(U(:,r)*V(:,r)',1,"V");      end %Reweight if kt FASTER
                     
                     [V(:,r), C(ii,2), C(ii,4)] = PERRIvsub(D-tmp,E,U(:,r),V(:,r),lamv,Vp(:,r),vtol,Nin,Nt,1); 
    %if frewe ;      [U(:,r),V(:,r)] = lsvdn(U(:,r)*V(:,r)',1,"U");      end %Reweight if kt FASTER
                     
                     
                     C(ii,[5 7]) = [PERRIcost(d,E,U(:,1:r),V(:,1:r),lamu,lamv,Up(:,1:r),Vp(:,1:r))      toc];
                     
                     if ii>(hh)
                         C(ii,6) = (C(ii,5) - C(ii-1,5))/C(ii,5);
                     end
                     
                 
                 svStr = fliplr(regexprep(fliplr(Suv),'\d{3}rN',fliplr(sprintf('Nr%03i',round(r))),'once')); %Flipping is to ensure only last instance of NrXX is matched in the string
              save(svStr,'U','V','C')
              save(Suv,'U','V','C') %Also save the most up to date version
                      if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',r,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii-hh+1); end
             
              
              
              if abs(C(ii,6))<1e-5 && ii>=(10+hh) %At least five iterations, and a relative change in gradient of less than 1e-5
                  C(ii+1:end,:) = repmat(C(ii,:),Nout*Nr-ii,1);
                  
                  
                  svStr = fliplr(regexprep(fliplr(svStr),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
                  U = Window2Ddata(U);
                  save(svStr,'U','V','C')
                  
                  fprintf('Reconstruction reached a tolerance with a gradient of %1.3e at iteration for rank %i %d\n',C(ii,6),ii-hh,r)
                  break
              elseif abs(C(ii,6))<5e-5 && ~ftol55 && ii>(1+hh) %At least ten iterations, and a relative change in gradient of less than 1e-5
                  svStr = fliplr(regexprep(fliplr(svStr),'tam.','tam.5-e5loT','once')); %Include a "Tol5e-5" save at the end
                  C(ii+1:end,:) = repmat(C(ii,:),(Nout*Nr)-ii,1);
                  save(svStr,'U','V','C')
                  ftol55 = 1;
              elseif abs(C(ii,6))<1e-4 && ~ftol4 && ii>(1+hh)%At least ten iterations, and a relative change in gradient of less than 1e-5
                  svStr = fliplr(regexprep(fliplr(svStr),'tam.','tam.4-e1loT','once')); %Include a "Tol1e-4" save at the end
                  C(ii+1:end,:) = repmat(C(ii,:),(Nout*Nr)-ii,1);
                  save(svStr,'U','V','C')
                  ftol4 = 1;
              elseif abs(C(ii,6))<1e-3 && ~ftol3 && ii>(1+hh)%At least ten iterations, and a relative change in gradient of less than 1e-5
                  svStr = fliplr(regexprep(fliplr(svStr),'tam.','tam.3-e1loT','once')); %Include a "Tol1e-3" save at the end
                  C(ii+1:end,:) = repmat(C(ii,:),(Nout*Nr)-ii,1);
                  save(svStr,'U','V','C')
                  ftol3 = 1;
              end
              %early break conditions
                    end
                 
                 tmp= E.mtimes2(U(:,1:r)*V(:,1:r)');
                 
                 
                 svStr = fliplr(regexprep(fliplr(Suv),'\d{3}rN',fliplr(sprintf('Nr%03i',round(r))),'once')); %Flipping is to ensure only last instance of NrXX is matched in the string
              save(svStr,'U','V','C')
              save(Suv,'U','V','C') %Also save the most up to date version
              
                 
                 %[U,V] = lsvdn(U*V',Nr,'U');
             end
         end
         
       
    
         
         
    case{'ftsm'} %full method with temporal smoothing
        %% Full rank recon
        C = zeros(Nout,7);
         
      %% PSF loop (no subproblem switching needed)
      if      fvpsf || fupsf 
            
            % Reconstruction
            if fupsf;   [V, C(:,2), C(:,4)] = PERRIvsubsmooth(D,E,U,V,lamv,Vp,lams,vtol,Nin,Nt,Nr);   end %lamu = inf, so just do the v subproblem
            if fvpsf;   [U, C(:,1), C(:,3)] = PERRIusub(      D,E,U,V,lamu,Up,     utol,Nin*10,Nx,Nr); end %lamv = inf, so just do the u subproblem
            
            % Calculate cost function
                        C(:,5) = PERRIcostsmooth(d,E,U,V,lamu,lamv,Up,Vp,lams);
                        C(:,6) = 0;
                        C(:,7) = toc;
            
            %Save output
            if fsave %For convenience, more than anything
                    for ii = svInd
                        svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                        save(svStr,'U','V','C')
                        save(Suv,'U','V','C') %Also save the most up to date version
                        svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
                        save(svStr,'U','V','C')
                    end
            end
   
            ii = Nout;
      else % The non-PSF loop (requires iterations
        

            
        % Load the truth for canonocal correlation checks
        if fdebug; [Ut,Vt]  = load_prior_n_data(strcat(Sd,'w'),'l',16);
                    Ut      = Ut.*mask;
                    if fk; Ut = U; Vt = V;
                    end
                    c       = zeros(Nout,4);
        end

        U = U.*mask;
        
        %% Iterations
        for ii = 1:Nout
             
            % Canoncorr debug
            if fdebug;      [uccfull ucc(ii)] = cssa(U,Ut,1);   
                            [vccfull vcc(ii)] = cssa(V,Vt);     
            end
             
            
            % Spatial subproblem
                            [U, C(ii,1), C(ii,3)] = PERRIusub(      D,E,U,V,lamu,Up,    utol,Nin*10,Nx,Nr);   %put back in *10  
            if frewe ;      [U,V] = lsvdn(U*V',Nr,"V");      end %Reweight if kt FASTER
            
            
            % Temporal subproblem
                            [V, C(ii,2), C(ii,4)] = PERRIvsubsmooth(D,E,U,V,lamv,Vp,lams,vtol,Nin,  Nt,Nr); 
            if frewe ;      [U,V] = lsvdn(U*V',Nr,"U");                                              end
                        
            % Calculate cost function (unbelievably no flags here)
                            C(ii,[5 7]) = [PERRIcostsmooth(d,E,U,V,lamu,lamv,Up,Vp,lams)     toc];
            if ii>1;        C(ii,6) = (C(ii,5) - C(ii-1,5))/C(ii,5);                                 end
            
            % Constituent parts of cost function for debugging
            if fdebug;      c(ii,1) = norm(vec(d-E*(U*V')), "fro" )^2;                              
                            c(ii,2) = lamu^2*norm(U - Up, "fro" )^2;                                
                            c(ii,3) = lamv^2*norm(V-Vp,"fro")^2;                                    
                            c(ii,4) = c(ii,1)+c(ii,2)+c(ii,3);   
            if ii>1;        c(ii,5) = (C(ii,5) - C(ii-1,5))/norm(d(:));                              end
            end

            % Print out current status
            if mod(ii,1)==0; if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end; end
        
            %[~,em] = componentThreshold(V,Sd,Nl,0.4,0)
            
            save(Suv,'U','V','C'); %Always save current version
            % Save intermediate iterations
            if sum(ii == svInd) && fsave
                svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(ii))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                save(svStr,'U','V','C')
            end
            

            if abs(C(ii,6))<1e-5 && ii>=10 %At least five iterations, and a relative change in gradient of less than 1e-5
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                for jj = (1+ii):1:Nout
                    if sum(jj == svInd) && fsave
                        svStr = fliplr(regexprep(fliplr(Suv),'\d{4}oN',fliplr(sprintf('No%04d',round(jj))),'once')); %Flipping is to ensure only last instance of NoXXXX is matched
                        save(svStr,'U','V','C')
                    end
                end
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.vnoC','once')); %Include a "Converged" save at the end
                U = Window2Ddata(U);
                save(svStr,'U','V','C')
                fprintf('Reconstruction reached a tolerance with a gradient of %1.3e at iteration %d\n',C(ii,6),ii)
                return
            elseif abs(C(ii,6))<5e-5 && ~ftol55 && ii>1 %At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.5-e5loT','once')); %Include a "Tol5e-5" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol55 = 1;
            elseif abs(C(ii,6))<1e-4 && ~ftol4 && ii>1%At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.4-e1loT','once')); %Include a "Tol1e-4" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol4 = 1;
            elseif abs(C(ii,6))<1e-3 && ~ftol3 && ii>1%At least ten iterations, and a relative change in gradient of less than 1e-5
                svStr = fliplr(regexprep(fliplr(Suv),'tam.','tam.3-e1loT','once')); %Include a "Tol1e-3" save at the end
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                save(svStr,'U','V','C')
                ftol3 = 1;
            end
                
        end
        
      end
        % if psf, pour everything into Nin.
        if fprint;    fprintf('%-3i\t%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t%-5i\t%-5i\t%-3i\n',Nr,C(ii,1),C(ii,2),C(ii,5),C(ii,6),C(ii,3),C(ii,4),round(C(ii,7)),ii); end

end
end
end
