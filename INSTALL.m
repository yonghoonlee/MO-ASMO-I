%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Automatically adds project files to your MATLAB path, downloads the required
% MATLAB Central submissions, checks available solvers, and opens an example.
%===============================================================================
function INSTALL
% Installation code adopted by INSTALL_DT_QP_Project function of DT_QP_Project
% by Daniel R. Herber from https://github.com/danielrherber/dt-qp-project
% (with modification by Yong Hoon Lee)
% See the license/INSTALL.License for detail.
    % AddSubmissionContents(mfilename) % add contents to path
    RequiredWebFiles % download required web files
    RequiredWebZips % download required web zips
    AddSubmissionContents(mfilename('fullpath')) % add contents to path
    %CheckSolvers % check available solvers
    %OpenThisFile('BrysonDenham') % open an example
    CloseThisFile(mfilename) % close this file
end
%===============================================================================
function RequiredWebFiles
    disp('Obtaining required web files')
    ind = 0; % initialize index
    files = struct('url',[],'folder',[]); % initialize structure
    %---------------------------------------------------------------------------
    % File 1
    ind = ind + 1; % increment
    files(ind).url = 'http://dmpeli.math.mcmaster.ca/Matlab/Math4Q3/Lecture2-1/LagrangeInter.m';
    files(ind).folder = 'LagrangeInter';
    %---------------------------------------------------------------------------
    % File 2
    ind = ind + 1; % increment
    files(ind).url = 'http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/lepoly.m';
    files(ind).folder = 'lepoly';
    %---------------------------------------------------------------------------
    full_fun_path = which(mfilename('fullpath')); % obtain full function path
    outputdir = fullfile(fileparts(full_fun_path),'external');
    DownloadWebFiles(files,outputdir); % download
end
%===============================================================================
function RequiredWebZips
    disp('Obtaining required web zips')
    ind = 0; % initialize index
    zips = struct('url',[],'folder',[],'test',[]); % initialize structure
    %---------------------------------------------------------------------------
    % Zip 1
	ind = ind + 1; % increment
	zips(ind).url = 'https://github.com/altmany/export_fig/archive/master.zip';
	zips(ind).folder = 'MFX23629_export_fig';
	zips(ind).test = 'export_fig';
    %---------------------------------------------------------------------------
    % Zip 2
	ind = ind + 1;
    zips(ind).url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/61768/versions/2/download/zip/hcparula.zip';
	zips(ind).folder = 'MFX61768_hcparula';
	zips(ind).test = 'hcparula';
    %---------------------------------------------------------------------------
    % Zip 3
	ind = ind + 1;
	zips(ind).url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/51986/versions/9/download/zip/Colormaps.zip';
	zips(ind).folder = 'MFX51986_colormaps';
	zips(ind).test = 'plasma';
    %---------------------------------------------------------------------------
    % Zip 4
    ind = ind + 1;
    zips(ind).url = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/40397/versions/7/download/zip/mfoldername_v2.zip';
    zips(ind).folder = 'MFX40397_mfoldername';
    zips(ind).test = 'mfoldername';
    %---------------------------------------------------------------------------
    % Zip 5
	ind = ind + 1;
	zips(ind).url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/60678/versions/2/download/zip/ypea126-nsga-iii.zip';
	zips(ind).folder = 'MFX60678_nsga3';
	zips(ind).test = 'nsga3';
    %---------------------------------------------------------------------------
    % zip 6
	ind = ind + 1;
	zips(ind).url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/10429/versions/9/download/zip/NSGA-II.zip';
	zips(ind).folder = 'MFX10429_nsga2';
	zips(ind).test = 'nsga2';
    %---------------------------------------------------------------------------
    % zip 7
	ind = ind + 1;
	zips(ind).url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/8773/versions/38/download/zip/Multiprod_2009.zip';
	zips(ind).folder = 'MFX8773_multiprod';
	zips(ind).test = 'multiprod';
    %---------------------------------------------------------------------------
    full_fun_path = which(mfilename('fullpath')); % obtain full function path
    outputdir = fullfile(fileparts(full_fun_path),'external');
    DownloadWebZips(zips,outputdir); % download and unzip
end
%===============================================================================
function AddSubmissionContents(name)
    disp('Adding submission contents to path');
    fullfuncdir = which(name); % current file
    submissiondir = fullfile(fileparts(fullfuncdir)); % current folder
    addpath(genpath(submissiondir)); % add folders and subfolders to path
end
%===============================================================================
function CloseThisFile(name)
    disp(['Closing ', name])
    h = matlab.desktop.editor.getAll; % get editor information
    for k = 1:numel(h) % go through all open files in the editor
        if ~isempty(strfind(h(k).Filename,name)) % if this is the file
            h(k).close % close this file
        end
    end
end
%===============================================================================
function OpenThisFile(name)
    disp(['--- Opening ', name])

    try
        % open the file
        open(name);
    catch % error
        disp(['Could not open ', name])
    end

    disp(' ')
end
%===============================================================================
function DownloadWebFiles(files,outputdir)

    % store the current directory
    olddir = pwd;
    
    % create a folder for outputdir
    if ~exist(outputdir, 'dir')
        mkdir(outputdir); % create the folder
    else
        addpath(genpath(outputdir)); % add folders and subfolders to path
    end
    
    % change to the output directory
    cd(outputdir)
    
    % go through each file
    for k = 1:length(files)
        
        % get data
        url = files(k).url;
        folder = files(k).folder;
        [~,nameurl,exturl] = fileparts(url);
        name = [nameurl,exturl];
        
        % first check if the test file is in the path
        if exist(name,'file') == 0
            
            try
                % download file
                outfilename = websave(name,url);
            
                % create a folder utilizing name as the foldername name
                if ~exist(fullfile(outputdir,folder), 'dir')
                    mkdir(fullfile(outputdir,folder));
                end

                % move the file
                movefile(outfilename,fullfile(outputdir,folder))

                % output to the command window
                disp(['Downloaded ',folder,'/',name])

            catch % failed to download
                % output to the command window
                disp(['Failed to download ',folder,'/',name])
                
                % remove the html file
                delete([name,'.html'])
            end
            
        else
            % output to the command window
            disp(['Already available ',name])
        end
    end
    
    % change back to the original directory
    cd(olddir)
end
%===============================================================================
function DownloadWebZips(zips,outputdir)

    % store the current directory
    olddir = pwd;
    
    % create a folder for outputdir
    if ~exist(outputdir, 'dir')
        mkdir(outputdir); % create the folder
    else
        addpath(genpath(outputdir)); % add folders and subfolders to path
    end
    
    % change to the output directory
    cd(outputdir)

    % go through each zip
    for k = 1:length(zips)

        % get data
        url = zips(k).url;
        folder = zips(k).folder;
        test = zips(k).test;

        % first check if the test file is in the path
        if exist(test,'file') == 0

            try
                % download zip file
                zipname = websave(folder,url);

                % save location
                outputdirname = fullfile(outputdir,folder);

                % create a folder utilizing name as the foldername name
                if ~exist(outputdirname, 'dir')
                    mkdir(outputdirname);
                end

                % unzip the zip
                unzip(zipname,outputdirname);

                % delete the zip file
                delete([folder,'.zip'])

                % output to the command window
                disp(['Downloaded and unzipped ',folder])
            
            catch % failed to download
                % output to the command window
                disp(['Failed to download ',folder])
                
                % remove the html file
                delete([folder,'.html'])
                
            end
            
        else
            % output to the command window
            disp(['Already available ',folder])
        end
    end
    
    % change back to the original directory
    cd(olddir)
end
%===============================================================================
function CheckSolvers
	disp('--- Checking available solvers')
	disp(' ')

    % to be completed 
end
%===============================================================================