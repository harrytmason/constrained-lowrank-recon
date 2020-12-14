function [X,T] = lsvdn(D,Nr,ID);
%% Produce version of X and T where X and T are normalised by equal distribution of singular values
%
%  Inputs: D: Dataset to perform lsvd on
%          Nr: Number of components to reduce D to
%
%  Outputs: X, T: The decomposition of D. X*T' = D
%
%  Requires: Mark Chiew's lsvd.m
%
%  Created by Harry Mason, University of Oxford, 24/05/2018
%
%% Initialisation

if nargin<3; ID = 'X'; end
if nargin<2; Nr = 16 ; end

%% Decomposition

[U,S,V] = lsvd(D,Nr);

switch ID
    case{'X','XS'}
        X = U*S;
        T = V;
        
    case{'S','XST'}
        X = U*sqrt(S);
        T = V*sqrt(S);
        
    case{'T','TS'}
        X = U;
        T = V*S;
        
    otherwise;  fprintf('No legitimate initial setting has been chosen. Check lsvdn.n for valid options.\n')
end

end