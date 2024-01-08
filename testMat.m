clc
clear variables
close all
initCobraToolbox(false)
addpath(fullfile('.','MATLAB'))

% Define arena
arena = defineArena(20, 20, 'seed', 10);

% Add organism
load(fullfile('.','test','models','Alistipes_putredinis_DSM_17216.mat'))
name = 'Alistipes_putrdinis_DSM_17216';
arena = addOrganism(arena, model, name, 3, 'biomassT0', 0.5, 'limitGrowth', true);

load(fullfile('.','test','models','Bacteroides_uniformis_ATCC_8492.mat'))
name = 'Bacteroides_uniformis_ATCC_8492';
arena = addOrganism(arena, model, name, 10, 'biomassT0', 0.5, 'limitGrowth', true);

load(fullfile('.','test','models','Barnesiella_intestinihominis_YIT_11860.mat'))
name = 'Barnesiella_intestinihominis_YIT_11860';
arena = addOrganism(arena, model, name, 10, 'biomassT0', 0.5, 'limitGrowth', true);

% Add media
media = readtable(fullfile('.','test','media','M1.xlsx'));
mets = table2cell(media(:,'REACTION'));
cant = cell2mat(table2cell(media(:,'mM')));
arena = addSubs(arena, mets, cant, 'unit', 'mM');

% Load coefficient matrix
bacMat = readtable(fullfile('.','test','coeffMat','cellCoeff.xlsx'), 'ReadVariableNames', false);
bacMat = table2array(bacMat);

nutMat = readtable(fullfile('.','test','coeffMat','nutCoeff.xlsx'), 'ReadVariableNames', false);
nutMat = table2array(nutMat);

% Simulate
data = simEnv(arena, 3, 'bacCoeff', bacMat, 'nutCoeff', nutMat, 'secObj', 'none');

% Plot microbe abundance
plotCellAbundance(data, 'rel', true)

% Find cross-feeding interactions
cross = findCrossFeeding(data, 3);