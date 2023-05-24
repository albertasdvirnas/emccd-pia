function [] = emccd_pia_gui()


    
    %% Add source files to path and initialize settings
    mFilePath = mfilename('fullpath');
    mfolders = split(mFilePath, {'\', '/'});
    [fd, fe] = fileparts(mFilePath);
    % read settings txt
    setsTable  = readtable(fullfile(fd,'emccd_settings.txt'),'Format','%s%s%s','ReadVariableNames', false);

    warning(''); % empty warning
    addpath(genpath(fullfile(mfolders{1:end - 1})));
    
    [~, lwid] = lastwarn;
    
    if strcmp(lwid, 'MATLAB:mpath:nameNonexistentOrNotADirectory')
        error('Unexpected error when asserting source folder path.')
    end

    %% Generate UI
    mFilePath = mfilename('fullpath');
    mfolders = split(mFilePath, {'\', '/'});
    versionEMCCD = importdata(fullfile(mfolders{1:end-1},'VERSION'));
        
    % create tabbed figure
    hFig = figure('Name', ['EMCCD-PIA GUI v',num2str(versionEMCCD)], ...
        'Units', 'normalized', ...
        'OuterPosition', [0 0 1 1], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none', ...
        'ToolBar', 'figure' ...
    );

    hPanel = uipanel('Parent', hFig);
    h = uitabgroup('Parent',hPanel);
    t1 = uitab(h, 'title', 'EMCCD-PIA');
    tsHCC = uitabgroup('Parent',t1);
    hPanelImport = uitab(tsHCC, 'title', 'Main tab');

    dotImport = uicontrol('Parent', hPanelImport, 'Style', 'edit','String',fullfile(pwd,'testfolder'),'Units', 'normal', 'Position', [0 0.9 0.5 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    set(dotImport, 'Min', 0, 'Max', 25)% limit to 10 files via gui;

    dotButton = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'Browse folder'},'Callback',@selection,'Units', 'normal', 'Position', [0.6 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    dotButtonFile = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'Browse file'},'Callback',@selection2,'Units', 'normal', 'Position', [0.7 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    dotButtonUI = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'uigetfiles'},'Callback',@selection3,'Units', 'normal', 'Position', [0.8 0.9 0.1 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]
    runButton = uicontrol('Parent', hPanelImport, 'Style', 'pushbutton','String',{'Run'},'Callback',@run,'Units', 'normal', 'Position', [0.7 0.2 0.2 0.05]);%, 'Max', Inf, 'Min', 0);  [left bottom width height]

    % Calculate moments

    % 
end

