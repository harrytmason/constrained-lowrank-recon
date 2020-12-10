function [Out1, Out2] = load_prior_n_data(ID, out_type, Nr)
%% A function for loading different prior sets and data sets
%
% A single function to handle all loading of information
%
% Inputs:
%       ID - A unique identifier of the dataset. Current IDs listed below.
%            Could be string to directory or numeric id
%       out_type - Choose the form of the output ('d'ata, 'p'rior, 't'ruth, 'ds'ata shell, 'ps'rior shell, 'l'svd decomposition, 'dID'ata ID, 'pID'rior ID  etc.)
%       Nr - The number of components to include (end rank) in U and V.
%
% Outputs:
%       Out1: If data, this will be the fully shaped output (3D)
%             If prior, this will be U_prior
%       Out2: If data, this will be the space by time output (2D)
%             If prior, this will be V_prior
%         If shell types, output zeros in shape or specified above
%         If strings types, the first output will be the string, the second output will be a flag
%
% Created by Harry Mason, FMRIB, University of Oxford, 2017
%%%%%%%%

%% Initialisation
if nargin<1; ID = 4;         end
if nargin<2; out_type = 'uv'; end
if nargin<3; Nr = 16;    end

%% Generate flags/Parse ID
ID_orig = ID;    %As we adapt the ID, keep an original copy.
redFlag = false; %Are we working with a reduced dataset?
winFlag = false; %do we want to window the corners of k-space?

if ischar(ID)
if ~isfile(ID) %if it's a valid file string, no parsing
    if strcmpi(ID(end),'r'); ID = ID(1:end-1); redFlag = true; end %e.g. ID = 2r would load a reduced 2nd dataset
    if strcmpi(ID(end),'w'); ID = ID(1:end-1); winFlag = true; end %e.g. ID = 2w would load a windowed 2nd dataset
    if ~isnan(str2double(ID)); ID = floor(str2double(ID)); end %convert to number
end
elseif isnumeric(ID)
    ID = floor(ID);
end

%% Loading
switch out_type
    case{'prior','p','uv','UV'}
        %% Just load a previous recon
        % Out1 = Spatial Information
        % Out2 = Temporal Information
        
        %% Flags for space prior and chronos prior
          if strcmpi(ID(1),'s') 
              fs = 1; fc = 0; 
              ID = regexprep(ID,'^s(pace)?(Priors)?','generatedPriors');
          elseif strcmpi(ID(1),'c')
              fs = 0; fc = 1;
              ID = regexprep(ID,'^c(hronos)?(Priors)?','generatedPriors');
          else
              fc = 1; fs = 1;
          end
        
        if ~isfile(ID)
            [ID] = load_prior_n_data(ID,'pID'); % Let's get the actual filename we'll be working with.
        end
        
        
        % Check ID is valid
        try Ld = tryload(ID);
        catch
            try Ld = tryload(strcat('DataHoldingPen/',ID));
                disp('loading from DHP');
            catch
                display( strcat('ID not found. ID given was', 32, ID_orig, ' and current directory is ', 32, pwd, '. Please double check input and location.') ); return
            end
        end
        % Sanitise Outputs
        if     isfield(Ld,'u_prior'); Out1 = Ld.u_prior; 
                                      Out2 = Ld.v_prior;
        elseif isfield(Ld,'U_prior'); Out1 = Ld.U_prior; 
                                      Out2 = Ld.V_prior;
        elseif isfield(Ld,'u');       Out1 = Ld.u;       
                                      Out2 = Ld.v;
        elseif isfield(Ld,'U');       Out1 = Ld.U;       
                                      Out2 = Ld.V;
        elseif isfield(Ld,'Up');      Out1 = Ld.Up;      
                                      Out2 = Ld.Vp;
        else
            display(strcat('No valid recon found in specified file',32,ID_orig))
            Out1 = 0; Out2 = 0;
        end

        Out1 = Out1*fs; 
        Out2 = Out2*fc;
        
    case{'prior_ID','pID','priorStr'}
        %% Output string to access prior
        % First output = prior string
        % Second ouput = does it exist in current workspace?
        
        if ~isfile(ID) %if it's a valid file string, no parsing
        if ~ischar(ID)
            % No return
            if      ID == -1 || nargin == 0; Out1 = 0; Out2 = 0; disp('No input selected'); return; end
            % ID selection
            if      ID == 0;   ID = 'Data_for_Harry'; %d050 also loaded
            elseif  ID == 1;   ID = 'NewPrior_1.mat';
            elseif  ID == 2;   ID = 'NewPrior_2.mat';
            elseif  ID == 3;   ID = 'NewPrior_3.mat';
            elseif  ID == 20;  ID = 'Hann_UVPrior20.mat'; %20 spoke prior for dataset 3 (task FMRI)
            elseif  ID == 25;  ID = 'priors_25pct.mat';
            elseif  ID == 40;  ID = 'Hann_UVPrior.mat'; %40 spoke prior for dataset 3 (task FMRI)
            elseif  ID == 50;  ID = 'priors_50pct.mat';
            elseif  ID == 420; ID = 'ktFASTER_Prior_20.mat';
            elseif  ID == 425; ID = 'ktFASTER_Prior_25.mat';
            elseif  ID == 450; ID = 'ktFASTER_Prior_50.mat';
            end
            
        elseif strcmpi(ID(1),'g')||strcmpi(ID,'generated') %Priors generated through k-t FASTER
            ID = regexprep(ID,'^g(enerated)?(Priors)?', 'generatedPriors');
            
        elseif strcmpi(ID(1),'l')||strcmpi(ID,'generated') %Priors generated through optimisation with no decomposition
            ID = regexprep(ID,'^l(sqr)?(Priors)?',      'lsqrPriors');
            
        elseif strcmpi(ID(1),'t')||strcmpi(ID,'truth') %Priors generated by blurring the ground truth
            ID = regexprep(ID,'^t(ruth)?(Priors)?',     'truthPriors');
            
        elseif strcmpi(ID(1),'r')||strcmpi(ID,'randomised') %Priors generated randomly
            ID = regexprep(ID,'^r(andomised)?(Priors)?','randomisedPriors');
            
        elseif strcmp(ID(1),'d')||strcmp(ID,'density') %Priors generated with a density-weighted recon. strcmp rather than strcmpi in case of Data_for_Harry.mat entry (although it should catch the full filename first)
            ID = regexprep(ID,'^d(ensity)?(Priors)?',   'densityPriors');
            
        elseif strcmp(ID(1),'0')||strcmp(ID(1),'z')||strcmp(ID,'zero') %Priors generated with a density-weighted recon. strcmp rather than strcmpi in case of Data_for_Harry.mat entry (although it should catch the full filename first)
            ID = regexprep(ID,'^[z0](ero)?(Priors)?',   'zeroPriors');
            
        elseif strcmp(ID(1),'s') %Spatial Priors generated with k-t FASTER, Temporal is zeros
            ID = regexprep(ID,'^[s](pace)?(Priors)?',   'spacePriors');
            
        elseif strcmp(ID(1),'c') %Temporal Priors generated with k-t FASTER, Spatial is zeros
            ID = regexprep(ID,'^[c](hronos)?(Priors)?', 'chronosPriors');
            
        end
        end
        
        ID = regexprep(ID,'.*(?<!.mat)$',strcat(ID,'.mat')); %allow for simpler inputs that don't have .mat endings.
        if isfile(fullfile(ID))
        else; disp('Priors not yet generated');
        end
        
        Out1 = ID;
        Out2 = isfile(fullfile(ID));
        
    case{'prior_shell','ps'}
        %% Load priors, create zeros in shape of priors
        [Out1, Out2] = load_prior_n_data(ID, 'p');
        Out1 = Out1*0;
        Out2 = Out2*0;
        
    case{'C','c'}
        [ID,~] = load_prior_n_data(ID,'dID'); %Replace ID with filename
        Ld     = tryload(ID,'C');
        Out1   = Ld.C;
        Out2   = size(Out1);
        
    case{'data','d'}
        %% Just load data
        % 1st output = x * y * t 3D data
        % 2nd output = (x*y) * t 2D data
        
        [ID,fk] = load_prior_n_data(ID,'dID'); %Replace ID with filename
        if strcmpi(ID,'phantom')
            %  data = vec(phantom(100));
            data = repmat(phantom,1,1); %phantom. Condensing this part since I rarely use it
            Out1 = data;
            Out2 = vec(data);
            return
            %for d1 = 1:12; data(:,d1) = data(:,d1)+sin(pi*d1*0.1)+cos(pi*d1*0.2);      end
        else; Ld = tryload(ID);
        end %load data
        
        % Sanitise loaded variables
        if isfield(Ld,'d') %Dynamic Sheep Logan Phantom (ID 5)
            data = Ld.d;
        elseif isfield(Ld,'data') %Dynamic Sheep Logan Phantom (ID 5)
            data = Ld.data;
            
        elseif isfield(Ld,'U_truth')
            if isfield(Ld,'S')
                Ut = Ld.U_truth * Ld.S;
                Vt = Ld.V_truth;
            else
                Ut = Ld.U_truth;
                Vt = Ld.V_truth;
            end
            numRows = factor(size(Ut,1));
            data = reshape(Ut*Vt',prod(numRows(1:2:end)),prod(numRows(2:2:end)),[]);
        elseif isfield(Ld,'U')
            if isfield(Ld,'S')
                U = Ld.U * Ld.S;
                V = Ld.V;
            end
            numRows = factor(size(U,1));
            data = reshape(U*V',prod(numRows(1:2:end)),prod(numRows(2:2:end)),[]);
        end
        
        % Reduce dataset size
        if redFlag
            data = data(floor(end/4)+1:ceil(end*0.75),floor(end/4)+1:ceil(end*0.75),1:ceil(end/8));
            % data = data(floor(end/32)+1:ceil(end*3/32),floor(end/32)+1:ceil(end*3/32),1:ceil(end/8));
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
        
        
    case{'data_shell','ds'}
        %% Load data, create zeros in shape of data
        [Out1, Out2] = load_prior_n_data(ID, 'd');
        Out1 = Out1*0;
        Out2 = Out2*0;
        
    case{'ones','1','d1','psf'}
        %% Load data, create ones in shape of data
        [Out1, Out2] = load_prior_n_data(ID, 'd');
        Out1 = Out1*0+1;
        Out2 = Out2*0+1;
        
    case{'truth','t','lsvd','l'}
        %% Load data in form of U_truth (Out1) and V_truth (Out2)
        [dID,proFlag] = load_prior_n_data(ID,'dID'); %Replace ID with filename
        if ~proFlag
        [~, data] = load_prior_n_data(ID_orig,'d');
        try addpath('MC_utils/misc'); catch; end %try is so that function can work with runtime, which doesn't like addpath
        
        [Out1, Out2] = lsvdn(data,Nr,'S');
        else
            try
                [Out1, Out2] = load_prior_n_data(sprintf('true%ddata.mat',ID),'uv');
                if Nr< size(Out1,2)
                    Out1 = Out1(:,1:Nr);
                    Out2 = Out2(:,1:Nr);
                end
            catch
            Out1 = 0;
            Out2 = 0;
            disp('Prospective, no underlying truth has been generated')
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
            elseif ID==1; ID = 'Truth_for_Harry.mat';
            elseif ID==2; ID = 'orig01_z40_t512.mat';
            elseif ID==3; ID = 'task_data_for_harry.mat';       %In k-space
            elseif ID==4; ID = 'im_4_Harry_abs.mat';            %images_For_Harry.mat contains complex and absolute versions. This is absolute only
            elseif ID==5; ID = 'DynamicSheppLoganBlurred.mat';
            elseif ID==6; ID = 'invivo_sens_harry.mat';         %Prospective k-space data with coils
            elseif ID==7; ID = 'im_4_Harry_abs.mat';            %images_For_Harry.mat contains complex and absolute versions. This is absolute only
            elseif ID==8; ID = 'reduced4data.mat';              %images_For_Harry.mat, but interpolated to small
            elseif ID==9; ID = 'reduced4dataequalS.mat';        %reduced4data.mat, but with the S values equalised
            elseif ID==10; ID = 'reduced4datalinearS.mat';      %reduced4data.mat, but with the S values in linear scale
            elseif ID==11; ID = 'reduced4datalogS.mat';         %reduced4data.mat, but with the S values in log scale
            elseif ID==12; ID = 'im_4_Harry_com.mat';           %Complex version of dataset 4
            elseif ID==13; ID = 'pro_data_30052019.mat';        %Task dataset acquired 30-May-2019. Mark's version with original variable names in ..._original.mat, slice 25 of overall dataset
            elseif ID==14; ID = 'long_pro_data_30052019.mat';    %The longer version of the same task, used as truth for task dataset pro_data_30052019. Mark's version in long_for_harry.mat
            elseif ID==15; ID = 'im_4_Harry_full.mat';            %A version of 4 which contains a fully sampled Cartesian k-space
            elseif ID==16; ID = 'im_4_Harry_cart32.mat';          %A version of 4 which contains a 32 line Cartesian k-space
            elseif ID==17; ID = 'im_4_Harry_cart20.mat';          %A version of 4 which contains a 20 line Cartesian k-space
            elseif ID>=18 && ID<=24; ID = sprintf('pro_data_30052019_slice_%02i.mat',ID); %pro_data_30052019_all_slices for full dataset
            elseif ID==26; ID = sprintf('pro_data_30052019_slices_18-25.mat'); %Concatenated slices 18-25.
            elseif ID==27; ID = 'reduced4data.mat';              %no coil reduced data
            elseif ID==28; ID = 'im_4_Harry_abs.mat';            %no coil retro data images_For_Harry.mat contains complex and absolute versions. This is absolute only
            elseif ID==29; ID = 'FMRIBletters.mat';            %no coil retro data images_For_Harry.mat contains complex and absolute versions. This is absolute only
            elseif ID==30; ID = 'toyFMRIBsmall.mat';            %Toy FMRIB dataset 50x50x300
            elseif ID==31; ID = 'im_4_Harry_full_2.mat';            %A version of 4 which contains a fully sampled Cartesian k-space
            elseif ID==32; ID = 'toyFMRIBsmall4.mat';            %Toy FMRIB dataset 64x64x300
            elseif ID==33; ID = 'im_4_Harry_cart50.mat'; 
            elseif ID==34; ID = 'im_4_Harry_cart40.mat';
            elseif ID==35; ID = 'im_4_Harry_cart38.mat'; 
            elseif ID==36; ID = 'im_4_Harry_cart36.mat';  
            elseif ID==37; ID = 'im_4_Harry_cart34.mat'; 
            elseif ID==38; ID = 'im_4_Harry_cart32b.mat';  %same as ID 16, just a different pseudorandom pattern outside the central sampling lines.
            elseif ID==39; ID = 'im_4_Harry_cart38b.mat';  %same as ID 35, but with 14 central sampling lines rather than 10

            end
        end
        if ~strcmpi(ID,'phantom')
            if ~isfile(ID); display( strcat('Data ID not found. ID given was', 32, ID, ' and current directory is ', 32, pwd, '. Please double check input and location.') ); end
        end
        proFlag = strcmpi(ID,'task_data_for_harry.mat')||strcmpi(ID,'invivo_sens_harry.mat')||contains(ID,'pro_data');
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
            elseif ID==1; ID = 'Truth_for_Harry.mat';
            elseif ID==2; ID = 'orig01_z40_t512.mat';
            elseif ID==3; ID = 'invivo_sens_harry.mat'; %Just for coils, we can use this dataset
            elseif ID==4; ID = 'invivo_sens_harry.mat'; %Coils work for retrospective as well as invivo dataset
            elseif ID==5; ID = 'sens_128.mat';
            elseif ID==6; ID = 'invivo_sens_harry.mat'; %Prospective k-space data with coils
            elseif ID==7; ID = 'invivo_sens_harry_3.mat'; %Prospective k-space data with coils
            elseif ID==8; ID = 'reduced4data.mat'; %reduced size coils
            elseif ID==9; ID = 'reduced4dataequalS.mat'; %reduced4data.mat, but with the S values equalised (so sensitivities are equal
            elseif ID==10; ID = 'reduced4datalinearS.mat'; %reduced4data.mat, but with the S values linearly ramped
            elseif ID==11; ID = 'reduced4datalogS.mat'; %reduced4data.mat, but with the S values logarithmically ramped
            elseif ID==12; ID = 'invivo_sens_harry.mat'; %Same as 4, but data it's paired with is complex
            elseif ID==13; ID = 'pro_data_30052019.mat';        %Task dataset acquired 30-May-2019. Mark's version with original variable names in ..._original.mat
                %'sens_for_harry_32'; is the full 32 coils for this dataset
            elseif ID==14; ID = 'long_pro_data_30052019.mat';   %The longer version of the same task, used as truth for task dataset pro_data_30052019. Mark's version in long_for_harry.mat
            elseif ID==15; ID = 'NoCoils';% 'im_4_Harry_full.mat';%'NoCoils';%            %A version of 4 which contains a fully sampled Cartesian k-space
            elseif ID==16; ID = 'NoCoils';%'im_4_Harry_cart32.mat';%'NoCoils';%           %A version of 4 which contains a 32 line Cartesian k-space
            elseif ID==17; ID = 'NoCoils';%'im_4_Harry_cart20.mat';%'NoCoils';%'NoCoils';%           %A version of 4 which contains a 20 line Cartesian k-space

            elseif ID>=18 && ID<=24; ID = sprintf('pro_data_30052019_slice_%02i.mat',ID); %pro_data_30052019_all_slices for full dataset
            elseif ID==27; ID = 'NoCoils';
            elseif ID==28; ID = 'NoCoils';  
            elseif ID==29; ID = 'NoCoils';
            elseif ID==30; ID = 'NoCoils';
            elseif ID==31; ID = 'im_4_Harry_full_2.mat';            %A version of 4 which contains a fully sampled Cartesian k-space
            elseif ID==32; ID = 'NoCoils';
            elseif ID==33; ID = 'NoCoils';
            elseif ID==34; ID = 'NoCoils';
            elseif ID==35; ID = 'NoCoils';
            elseif ID==36; ID = 'NoCoils';
            elseif ID==37; ID = 'NoCoils';
            elseif ID==38; ID = 'NoCoils';
            elseif ID==39; ID = 'NoCoils';
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
            if contains(ID,'reduced4data')
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
        % 2nd output = The representation of the decimal * 1000, if it exists (start frame in a reduced DOF recon)
        
        if isnumeric(ID)
            Out1 = num2str(ID);
        else
            Out1 = ID;
        end
        
        if ischar(ID_orig)
            if ~isfile(ID_orig) %if it's a valid file string, no parsing
                if strcmpi(ID_orig(end),'r'); ID_orig = ID_orig(1:end-1); end %e.g. ID = 2r would load a reduced 2nd dataset
                if strcmpi(ID_orig(end),'w'); ID_orig = ID_orig(1:end-1); end %e.g. ID = 2w would load a windowed 2nd dataset
                if ~isnan(str2double(ID_orig)); Out2 = str2double(ID_orig)-floor(str2double(ID_orig)); end %convert to number
            end
        elseif isnumeric(ID_orig)
            Out2 = ID-floor(ID_orig);
        else
            Out2 = 0;
        end
            Out2 = round(Out2*1000); 
        
        
    otherwise
        disp('No correct option selected as second input')
end


function a=tryload(varargin)
% Because a lot of reconstructions run at once, this gives a two chances to
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
        
