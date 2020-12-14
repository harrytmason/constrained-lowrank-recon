function [Out1, Out2] = load_prior_n_data(ID, out_type, Nr)
%% A function for loading different prior sets and data sets
%
% A single function to handle all loading of information
%
% Inputs:
%       ID - A unique identifier of the dataset. Current IDs listed below.
%            Could be string to directory or numeric id
%       out_type - Choose the form of the output ('d'ata, 't'ruth, 'ds'ata shell, 't'ruth, 'dID'ata ID, etc.)
%       Nr - The number of components to include (end rank) in X and T.
%
% Outputs:
%       Out1, Out 2. The form is dependent on the input
%
% Created by Harry Mason, FMRIB, University of Oxford, 2020
%
%%%%%%%%

%% Initialisation
if nargin<1; ID = '1';         end
if nargin<2; out_type = 'd'; end
if nargin<3; Nr = 16;    end

%% Generate flags/Parse ID
ID_orig = ID;    %As we adapt the ID, keep an original copy.
winFlag = false; %do we want to window the corners of k-space?

if ischar(ID)
if ~isfile(ID) %if it's a valid file string, no parsing
    if strcmpi(ID(end),'w'); ID = ID(1:end-1); winFlag = true; end %e.g. ID = 2w would load a windowed 2nd dataset
    if ~isnan(str2double(ID)); ID = floor(str2double(ID)); end %convert to number
end
elseif isnumeric(ID)
    ID = floor(ID);
end

%% Loading
switch out_type
    case{'data','d'}
        %% Just load data
        % 1st output = x * y * t 3D data
        % 2nd output = (x*y) * t 2D data
        
        [ID,fk] = load_prior_n_data(ID,'dID'); %Replace ID with filename
        Ld = tryload(ID); %load data
        
        % Sanitise loaded variables
        if  isfield(Ld,'X') %if loading from data which is stored as X and T
            if isfield(Ld,'S')
                X = Ld.X * Ld.S;
                T = Ld.T;
            else
                X = Ld.X;
                T = Ld.T;
            end
            numRows = factor(size(X,1));
            data = reshape(X*T',prod(numRows(1:2:end)),prod(numRows(2:2:end)),[]);
        else
            data = Ld.data;
        end
              
        if winFlag %for retrospective - prospective already windowed
            if ~fk
                [~,mask] = load_prior_n_data(ID_orig,'m');
                data = Window2Ddata(data).*mask; %Applies window in k-space
            end
        end
        
        Nd = size(data);
        
        Out1 = data;
        Out2 = reshape(data,[],Nd(end));
        
        
    case{'truth','t','lsvd','l'}
        %% Load data in form of X_truth (Out1) and T_truth (Out2)
        [dID,proFlag] = load_prior_n_data(ID,'dID'); %Replace ID with filename
        if ~proFlag
        [~, data] = load_prior_n_data(ID_orig,'d');
        
        [Out1, Out2] = lsvdn(data,Nr,'S');
        else
            try
                [~, data] = load_prior_n_data(sprintf('true%ddata.mat',ID),'d');
                [Out1, Out2] = lsvdn(data,Nr,'S');
            catch
            Out1 = 0;
            Out2 = 0;
            disp('No underlying truth has been generated')
            end
        end
        
    case{'data_ID','dID','dataStr'}
        %% Output string to access data
        % First output = data string
        % Second ouput = pro flag (is data prospective). Only dataset 3 triggers so far
        if strcmpi(ID,'phantom');  ID = 0; end % allow dataset descriptions
        if strcmpi(ID,'FMRIB');    ID = 1; end % as inputs
        if ~ischar(ID)
            if     ID==0 || strcmpi(ID,'phantom'); ID = 'data0.mat';%ID = 'phantom';
            
            elseif ID==1; ID = 'Datasets/reducedRetroData.mat';              %Retrospective Dataset A from the paper, but interpolated to be small
            elseif ID==2; ID = 'Datasets/retroDataA.mat';                             %Retrospective Dataset A from the paper
            elseif ID==3; ID = 'Datasets/prospectiveDataSlice1.mat';         %Slice 1 of the prospective dataset from the paper
            end
        end

        proFlag = contains(ID,'prospective');
        Out1 = ID;
        Out2 = proFlag;
        
    case{'k', 'kspace'}
        %% Load k-space indices
        % Out1 is the indicies in k-space, Out2 is a flag for if those indicies exist
        try if strcmpi(ID(end),'r'); ID = str2double(ID(1:end-1));end; catch; end %e.g. ID = 2r would load a reduced 2nd dataset
        
        [ID,~] = load_prior_n_data(ID,'dID'); %Replace ID with filename
        Ld   = tryload(ID);
        if isfield(Ld,'k')
            Out1 = Ld.k;
            Out2 = true; %k-space indecies are included
        else
            disp('No precalculated k indicies available, calculating by hand');
            Out1 = 0;
            Out2 = false;
        end
        
    case{'s', 'sens' , 'sensitivity'}
        %% Load sensitivity matrix
        % Out1 is the sensitivity matrix, 
        % Out2 is a flag for if the sensitivity matrix exists
        
        [ID,fsID] = load_prior_n_data(ID,'sID'); %Replace ID with filename
        if ~fsID
%            disp('No sensitivity variable available');
            Out1 = [];
            Out2 = false;
        else
            Ld = tryload(ID);
            if isfield(Ld,'sens')
                Out1 = Ld.sens;
                Out2 = true; %sFlag! is there a sensitivity matrix
            else
                disp('No sensitivity variable available');
                Out1 = [];
                Out2 = false;
            end
        end
            
        
    case{'sens_ID','sID','sensStr'}
        %% Output string to access prior
        % First output = prior string
        % Second ouput = does it exist in current workspace?
        
        if ~isnan(str2double(ID)) ;ID = str2double(ID); end
        if ~ischar(ID)
            % ID selection
            if     ID==0 || strcmpi(ID,'phantom'); ID = 'data0.mat'; %ID = 'phantom'
            elseif ID==1; ID = 'Datasets/reducedRetroData.mat'; %reduced size coils
            elseif ID==2; ID = 'Datasets/retroDataASensitivities.mat'; %Coils work for retrospective as well as invivo dataset
            elseif ID==3; ID = 'Datasets/prospectiveDataSlice1.mat';         %Slice 1 of the prospective dataset from the paper
            end
                       
        else
            if isfile(fullfile(ID))
            else; disp('Requested file containing Sensitivity Matrix doesn''t exist');
            end
           
        end
        
        Out1 = ID;
        Out2 = isfile(fullfile(ID));
        
        
    case{'m', 'mask' }
        %% Load sensitivity mask (which elements are non-zero in the sensitivity matrix)
        % Out1 is the sensitivity matrix, 
        % Out2 is a flag for if the sensitivity matrix exists
        
        [sens,fs] = load_prior_n_data(ID,'s'); %Replace ID with filename
        [Nx]   = load_prior_n_data(ID,'sz');
        
        if ~fs; Out1 = vec(ones(Nx(1:2))); %if there's no coil sensitivity info to mask
        else;   Out1 = vec(sum(abs(sens),4)~=0);
        end
        
        [~,fk] = load_prior_n_data(ID,'dID');

        if fk;  Out2 = reshape(Out1,Nx(1)/2,Nx(1)/2); 
        else;   Out2 = reshape(Out1,Nx(1),  Nx(2)  );
        end
        
    case{'size','Nd','sz'}
        %% Size of data
        % 1st output = [sqrt(Nx) sqrt(Nx) Nt] (if retrospective); [Nk,Nl,Nc] (if prospective)
        % 2nd output = [Nx Nt] (if retrospective); [Nk*Nl,Nc] (if prospective)
        
        [d,d2] = load_prior_n_data(ID,'d'); %Replace ID with filename
        Out1 = size(d);
        Out2 = size(d2);
        
    case{'TR'}
        %% TR of dataset
        % 1st output = TR
        % 2nd Output = max freq
        
        [ID,fk] = load_prior_n_data(ID,'dID');
        if fk
            Out1 = 0.05; %TR per scan line - each line took roughly 0.05 seconds
        else
            Out1 = 1;
            if contains(ID,'reduced')
                Out1 = 6;
            end
        end
        Out2 = 1/Out1; % Maximum frequency
        
    case{'angles','a'}
        %% Angles of k-space lines
        % 1st output = List of angles
        % 2nd output = flag of if list of angles exist
        
        [ID,~] = load_prior_n_data(ID,'dID');
        Ld = tryload(ID);
        if isfield(Ld,'angles')
            Out1 = Ld.angles;
            Out2 = 1;
        else
            Out1 = 0;
            Out2 = 0;
        end
        
    case{'Sd','ID'}
        %% Parse the ID
        % 1st output = Clean ID
        % 2nd output = 0 (not used)
        
        if isnumeric(ID)
            Out1 = num2str(ID);
        else
            Out1 = ID;
        end

            Out2 = 0;        
        
    otherwise
        disp('No correct option selected as second input')
end


function a=tryload(varargin)
% If a lot of reconstructions were running at once, this gives a two chances to
% load a file - one after a slight randomised loading delay
pause(0.5)
if exist(varargin{1},'file') %only try loading again if the file actually exists though

loadStr = strcat("load('",varargin{1});
if nargin>1
    for i1 = 2:nargin
        loadStr = strcat(loadStr,''',''',varargin{i1});
    end
end
loadStr = strcat(loadStr,"')");

try a=eval(loadStr);
catch
    pause(randi([5 50]));
    a=eval(loadStr);
end
end
end

end
        
