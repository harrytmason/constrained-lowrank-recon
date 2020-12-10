function lambda = genLam(Urange,Vrange,spacing,zero_edge, inf_edge, Srange, UeqV)
%% genLam.m
%
% Purpose: A simple function to create a 2 column matrix where each lambda
% pair is a unique combination within the xrange and yrange, and are
% ordered sensibly.
%
% Inputs:   Urange: a 2-element array with the 1e(xrange(1)) > 1e(xrange(2))
%           Vrange: as above, but the range of lambda_V rather than lambda_U
%           spacing: the spacing to use between log powers of U_range and V_range. Can be a scalar, or a vector to refer to [Uspace Vspace Sspace]
%           zero edge: boolean of whether to put zeros around the lowest edge of U and V (if 2D, the first element is for u, the second is for v)
%           inf edge: boolean of whether to put infs around the highest edge of U and V (if 2D, the first element is for u, the second is for v)
%           Srange: Like Urange and V range, but for smoothing
%           UeqV: A flag to restrict the same values of U and V in output lambda (i.e. for Tikhonov, Vrange input will be ignored) 
%
% Outputs:  lambda. Col 1 = lambda_U values, Col 2 = lambda_V values, Col 3 = lambda_S values
%
%%% Function by Harry Mason, 2018, University of Oxford

%% Initialisation
%
if nargin<7; UeqV      = 0; end
if nargin<6; Srange    = [];end
if nargin<5; inf_edge  = 0; end
if nargin<4; zero_edge = 0; end
if nargin<3; spacing   = 1; end

%% Creating lamU and lamV and lamS
%
if numel(Urange)>0
    lamU = 10.^(Urange(1):spacing(1):Urange(end));
else
    lamU = [];
end
%
if numel(Vrange)>0
    
    if numel(spacing  )==1; spacing   = [spacing   spacing  ]; end %Little mod to allow 1D inputs
    if numel(inf_edge )==1; inf_edge  = [inf_edge  inf_edge ]; end
    if numel(zero_edge)==1; zero_edge = [zero_edge zero_edge]; end
    
    lamV = 10.^(Vrange(1):spacing(2):Vrange(end));
    
else
    lamV = []; 
    if numel(spacing  )==1; spacing   = [spacing   0 spacing  ]; end %Little mod to allow 1D inputs
    if numel(inf_edge )==1; inf_edge  = [inf_edge  0 inf_edge ]; end
    if numel(zero_edge)==1; zero_edge = [zero_edge 0 zero_edge]; end
end
%
if numel(Srange)>0
lamS = 10.^(Srange(1):spacing(end):Srange(end));
else
    lamS = [];
end


if isnan(sum(lamU)); lamU = []; end
if isnan(sum(lamV)); lamV = []; end

%% Adding leading zeros?
%
if zero_edge(1);    lamU = [0    lamU];      end
if zero_edge(2);    lamV = [0    lamV];      end
if zero_edge(end);  lamS = [0    lamS];      end

%% Add trailing zeros?
%
if inf_edge(1);     lamU = [     lamU inf];  end
if inf_edge(2);     lamV = [     lamV inf];  end
if inf_edge(end);   lamS = [     lamS inf];  end

%% Final dimension check (easy than adding if conditions above, go ahead and sue me)
%
if numel(Urange)==0;  lamU=0;  end
if numel(Vrange)==0;  lamV=0;  end
if numel(Srange)==0;  lamS=[];  end

%% Final number of lambda_u, lambda_v, and lambda_s values
%
Ulen = length(lamU);
Vlen = length(lamV);
Slen = length(lamS);

%% Create n*2 array, where each row is a unique combination of lambda_u and lambda_v
%
if UeqV
lambda=  [lamU' lamU']; Vlen=1;
else
    
lambda = [repmat(lamU,1,Vlen); reshape(repmat(lamV,Ulen,1),[],1)']';
end

%% if numel(lamS)>0, then add a third column for smoothing

if Slen>0
lambda = [repmat(lambda,Slen,1) reshape(repmat(lamS,Ulen*Vlen,1),[],1)];
end
%% End! Simple right?