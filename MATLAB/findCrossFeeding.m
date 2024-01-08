function res = findCrossFeeding(data, time, varargin)

% findCrossFeeding.m
%
% Script to extract cross-feeding interactions in a given time point.
%
% INPUTS:
%
%   - data:  Simulation data object.
%   - time:  Time point.
%
% OPTIONAL INPUTS:
%
%   - mets:     Set of metabolites to consider.
%   - cutOff:   Threshold for numeric accuracy.
%
% OUTPUTS:
%
%   - res: Different cross-feeding interactions.
%
% EXAMPLE:
%
%   res = findCrossFeeding(data, 3, 'mets', {'EX_ac(e)'; 'EX_for(e)'}, ...
%                          'cutOff', 1e-06)
%
% .. Authors: 
%       - Telmo Blasco, 08/06/2023, University of Navarra, TECNUN School of Engineering.

% Manage arguments
parser = inputParser;
addRequired(parser, 'data', @(x) isstruct(x))
addRequired(parser, 'time', @(x) isnumeric(x))
addParameter(parser, 'mets', data.exRxns, @(x) iscell(x))
addParameter(parser, 'cutOff', 1e-06, @(x) isnumeric(x) && x>0)

% Extract argument values
parse(parser, data, time, varargin{:});
data = parser.Results.data;
time = parser.Results.time;
mets = parser.Results.mets;
cutOff = parser.Results.cutOff;

% Check time step
time = time + 1;
if length(data.mfluxList) < time
    error('Simulation did not arrive at the desired time step')
end

% Filter mets
mets = mets(ismember(mets, data.exRxns));
nM = length(mets);

% Generate flux matrix
mflux = data.mfluxList{time};
mfluxMat = zeros(nM, length(data.orgID));
for i = 1:length(data.orgID)
    metsOrg = mets(ismember(mets, data.orgExch{i}));
    idxOrg = cellfun(@(x) find(ismember(data.orgExch{i}, x)), metsOrg);
    idxFlux = cellfun(@(x) find(ismember(mets, x)), metsOrg);
    mfluxMat(idxFlux, i) = mflux{i}(idxOrg);
end

% Generate interaction matrix
res = [];
for i = 1:nM
    if any(mfluxMat(i,:)>cutOff) && any(mfluxMat(i,:)<-cutOff)
        prods = find(mfluxMat(i,:)>cutOff);
        cons = find(mfluxMat(i,:)<-cutOff);
        [pMesh, cMesh] = meshgrid(prods, cons);
        
        newIt = [data.orgName(pMesh(:)), data.orgName(cMesh(:)), repmat(mets(i), length(prods)*length(cons), 1), ...
            num2cell(mfluxMat(i,pMesh(:))'), num2cell(mfluxMat(i,cMesh(:))')];
        res = [res; newIt];
    end
end

end