function plotCellAbundance(data, varargin)

% plotCellAbundance.m
%
% Script to plot cell type evolution through time.
%
% INPUTS:
%
%   - data:  Simulated environment structure.
%
% OPTIONAL INPUTS:
%
%   - org:  Cell array of organism names to be considered.
%   - rel:  Logical value to plot the counts as relative or absolute.
%           Default: false.
%
% EXAMPLE:
%
%   plotCellAbundance(data, {'Abiotrophia_defectiva_ATCC_49176'}, true)
%
% .. Authors: 
%       - Telmo Blasco, 08/11/2023, University of Navarra, TECNUN School of Engineering.

% Manage arguments
parser = inputParser;
addRequired(parser, 'data', @(x) isstruct(x))
addParameter(parser, 'org', [], @(x) isempty(x) || iscell(x))
addParameter(parser, 'rel', false, @(x) islogical(x))

% Extract argument values
parse(parser, data, varargin{:});
data = parser.Results.data;
org = parser.Results.org;
rel = parser.Results.rel;

% Number of iterations
numIt = length(data.simList);

% Number of organisms
n = length(data.orgID);

% Extract organism names
orgName = cellfun(@strrep, data.orgName, repmat({'_'},n,1), repmat({' '},n,1), 'Un', 0);

% Extract cell abundances across time
cellAbundance = cell(n,1);
for i = 1:n
    cellAbundance{i,1} = cellfun(@(x) sum(ismember(x(:,1),data.orgID(i))), data.simList);
end

if rel
    % Extract total number of cells across time
    totalAbundance = cellfun(@(x) size(x,1), data.simList);
    
    % Transform abundances into relative counts
    cellAbundance = cellfun(@(x) x./totalAbundance, cellAbundance, 'Un', 0);
end

% Filter by organism name
if ~isempty(org)
    idx = ismember(data.orgName,org);
    if any(idx)
        cellAbundance = cellAbundance(idx);
        orgName = orgName(idx);
        n = length(orgName);
    else
        error('Input organism names not found in the arena')
    end
end

% Plot cell abundance across time
figure()
for i = 1:n
    plot(0:(numIt-1),cellAbundance{i,1})
    hold on
end

% Plot settings
xlabel('Time [h]')
title('Organism Abundance')
set(gca,'XTick',(0:1:(numIt-1)))
if rel
    ylim([0,1])
    ylabel('Relative counts')
else
    ylabel('Absolute counts')
end
legend(orgName, 'Location', 'eastoutside')
grid on

end