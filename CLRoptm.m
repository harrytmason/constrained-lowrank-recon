function [X, T, C] = CLRoptm(d, D, E, X,T, lam,Xp,Tp, N, f,S)

%% The optimisation loop for Constrained Low-Rank Optimisation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function carries out alternating optimisation of the X (Spatial) and
% T (Temporal) subproblems, with various styles of alternation.
%
% Low rank conditions are enforced, as is regularisation with respect to a
% prior.
%
% Used in prior generation and also in final reconstructions. Capable of
% k-t FASTER, k-t PSF and Constrained Low-Rank reconstructions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:  D:   nufft_adj(k-space data)                                  (Nk , Nt)
%          E:   nufft structure + more                                   (xfm_nufft structure)
%          X:   current guess of X, the spatial subspace                 (Nx , Nr)
%          T:   current guess of T, the temporal subspace                (Nt , Nr)
%          lam: the weighting on the priors                              [lamx, lamt, lams]
%          Xp:  the prior of X                                           (Nx , Nr)
%          Tp:  the prior of T                                           (Nt , Nr)
%          N.stol: tolerance of subproblems                               scalar
%          N.tol:  convergence tolerance                                  scalar
%          N.r:    the rank of matrices X and T                           scalar
%          N.in:   number of iterations of subproblem                     scalar
%          N.out:  number of switches between X and T subproblems         scalar
%          f.psf: whether we're doing a psf recon (so not updating T)     boolean
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs: X:   final guess of X, the spatial subspace                  (Nx , Nr)
%          T:   final guess of T, the temporal subspace                 (Nt , Nr)
%          C:   Seven columns of optimisation information               (Nr * Nout, 7)
%            C(n,1): The nth residual of the X subproblem
%            C(n,2): The number of times the X subproblem iterated in the nth pass
%            C(n,3): The nth residual of the T subproblem
%            C(n,4): The number of times the T subproblem iterated in the nth pass
%            C(n,5): The total cost function after n iterations
%            C(n,6): The error with respect to a ground truth after n iterations
%            C(n,7): The time taken after n iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Required Functions: CLRcost.m, the underlying cost function
%                     CLRxsub.m, the expression of the x (spatial) subproblem
%                     CLRtsub.m, the expression of the t (temporal) subproblem
%                     Dr. Mark Chiew's xfm_nufft structures (represented by E in the code)
%                     The irt toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by Harry Mason, M.Eng, University of Oxford
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALISATION

if nargin < 6;      lamx   = 0;
                    lamt   = 0;
                    lams   = 0;
else;               lamx   = lam(1);
                    if numel(lam)==1
                            lamt = lam(end);
                            lams = 0;
                    elseif numel(lam)>1
                            lamt = lam(2);
                        if numel(lam)>2
                            lams = lam(3);
                        else
                            lams = 0;
                        end
                    end
end
if nargin < 7;      Xp     = zeros(size(X));    end
if nargin < 8;      Tp     = zeros(size(T));    end
if nargin < 9;      xtol   = 1e-15;   % sub problem tolerence
                    ttol   = 1e-15;
                    Ntol   = 1e-5;
                    Nr     = size(X,2);
                    Nin    = 10;
                    Nout   = 10;
                    Nx     = size(X,1);
                    Nt     = size(T,1);
else;               xtol   = N.stol;     % Accessing from a structure takes more time, so we extract all the values here.
                    ttol   = N.stol;
                    Ntol   = N.tol;
                    Nr     = N.r;
                    Nin    = N.inner;
                    Nout   = N.outer;
                    Nx     = N.x;
                    Nt     = N.t;
end
if nargin < 10;     fpsf  = false;
else;               fpsf  = f.psf;
end
if nargin < 11;     Sxt = 'XTNo0000.mat'; 
                    Sd  = '';
else; fsave = 1;       
    if isstruct(S); Sxt = S.xt;
                    Sd  = S.d;
    else
                    Sxt = S;
                    Sd = '';
    end  
end

svInd = [1:5 10:10:100 Nout]; %Which intermediate iterations to print out on

mask = load_prior_n_data(Sd,'m'); 

%% PRINT HEADER

fprintf('%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-5s\n','Cost Func X','Cost Func T','Cost Func','CF Grad','Time (s)','It Loop'); 
tic

%% FLAG CHECKS

    frewe = (sum(lamx+lamt)==0) && ~fpsf; %If we have any non-zero lambda, we don't need to reweight (psf flags as security for any infinities)
    C = zeros(Nout,7);
    X = X.*mask;
    
%% RECONSTRUCTION   
         
      %% PSF loop (no need to alternate subproblems)
      if fpsf
            
            % Reconstruction
            [X, C(:,1), C(:,3)] = CLRxsub(D,E,X,T,lamx,Xp,     xtol,Nin*10,Nx,Nr);
            
            % Calculate cost function
                        C(:,5) = CLRcost(d,E,X,T,lamx,lamt,lams,Xp,Tp);
                        C(:,6) = 0;
                        C(:,7) = toc;
                        
            fprintf('%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t\t%-3i\n',C(1,1),C(1,2),C(1,5),C(1,6),round(C(1,7)),1);

            %Save output
            % save(Sxt,'X','T','C') %If you want to run this function outside of the test script, uncomment this line to save the output 
   
      else % The non-PSF loop (requires iterations)
       
        %% Iterations
        for ii = 1:Nout

            % Spatial subproblem
                            [X, C(ii,1), C(ii,3)] = CLRxsub(D,E,X,T,lamx,Xp,    xtol,Nin*10,Nx,Nr);   %put back in *10  
            if frewe ;      [X,T] = lsvdn(X*T',Nr,"T");      end %Reweight if kt FASTER
            
            % Temporal subproblem
                            [T, C(ii,2), C(ii,4)] = CLRtsub(D,E,X,T,lamt,Tp,lams,ttol,Nin,  Nt,Nr); 
            if frewe ;      [X,T] = lsvdn(X*T',Nr,"X");                                              end
                        
            % Calculate cost function 
                            C(ii,[5 7]) = [CLRcost(d,E,X,T,lamx,lamt,lams,Xp,Tp)     toc];
            if ii>1;        C(ii,6) = (C(ii,5) - C(ii-1,5))/C(ii,5);  end % The convergence test (normalised cost function gradient)                                end            

            % Print out current status
            if sum(ismember(svInd,ii)); fprintf('%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t\t%-3i\n',C(ii,1),C(ii,2),C(ii,5),C(ii,6),round(C(ii,7)),ii);  end
                   
            %Save intermediate iterations
            if sum(ii == svInd) && fsave
                % save(Sxt,'X','T','C') %If you want to run this outside of the test script, uncomment this line to save the output
            end
            
            % Convergence catch
            if abs(C(ii,6))<Ntol && ii>=2 %At least two iterations, and a relative change in gradient of less than Ntol
                C(ii+1:end,:) = repmat(C(ii,:),Nout-ii,1);
                X = Window2Ddata(X);
                % save(Sxt,'X','T','C') %If you want to run this outside of the test script, uncomment this line to save the output
            if ~sum(ismember(svInd,ii));    fprintf('%-10.6e\t%-10.6e\t%-10.6e\t%-10.6e\t%-5i\t\t%-3i\n',C(ii,1),C(ii,2),C(ii,5),C(ii,6),round(C(ii,7)),ii); end 
                fprintf('Reconstruction reached a tolerance of %1.3e with a gradient of %1.3e at iteration %d\n\n',Ntol, C(ii,6),ii)
                return
            end
                
        end
        
      end
       
fprintf('\n');
end


