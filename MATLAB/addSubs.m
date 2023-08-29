function arena = addSubs(arena, rxns, cant, varargin)

% addSubs.m
%
% Script to add specific substances to the environment.
%
% INPUTS:
%
%   - arena:        Environment structure.
%   - rxns:         Name of the exchange reactions of the culture media.
%   - cant:         Numeric vector indicating the substances concentration.
%
% OPTIONAL INPUTS:
%
%   - unit:         Chemical unit of the amount of substances to be added.
%                   Options: mM, mmol/cm2, mmol/arena, mmol/cell,
%                   fmol/cell. Default: mM.
%   - diffSpeed:    Number indicating the diffusion rate of the substance.
%                   Default (Glucose): 0.02412 cm^2/h.
%
% OUTPUTS:
%
%   - arena: Environment structure with the information of the culture media.
%
% EXAMPLE:
%
%   - arena = addSubs(arena, {'EX_h2o(e)'; 'Ex_glc_D(e)'}, [0.5; 0.01], ...
%             'unit', 'mM')
%
% .. Authors: 
%       - Telmo Blasco, 16/12/2022, University of Navarra, TECNUN School of Engineering.

% Manage arguments
parser = inputParser;
addRequired(parser, 'arena', @(x) isstruct(x))
addRequired(parser, 'rxns', @(x) iscell(x))
addRequired(parser, 'cant', @(x) isnumeric(x) && all(x)>=0)
addParameter(parser, 'unit', 'mM', @(x) ischar(x))
addParameter(parser, 'diffSpeed', 0.02412, @(x) isnumeric(x) && all(x)>0)

% Extract argument values
parse(parser, arena, rxns, cant, varargin{:});
arena = parser.Results.arena;
rxns = parser.Results.rxns;
cant = parser.Results.cant;
unit = parser.Results.unit;
diffSpeed = parser.Results.diffSpeed;

% Check if organisms are already in the arena
if ~isfield(arena,'models')
    error('Organisms need to be defined first to determine what substances can be exchanged')
end

% Check substrate concentrations
if (length(cant) ~= length(rxns)) && (length(cant)~=1)
    error('The parameter cant should be of the same size of rxns or equal to 1')
end
if length(cant) == 1
    cant = repmat(cant, length(rxns), 1);
end

% Check unit values
if ~any(strcmp(unit,{'mM','mmol/cm2','mmol/arena','mmol/cell','fmol/cell'}))
    error('Wrong unit parameter. Options are: {mM, mmol/cm2, mmol/arena, mmol/cell, fmol/cell}')
end

% Transform units to fmol/cell
switch unit
    case 'mM'
        conv = (10^12) * 0.01 * arena.scale;
    case 'mmol/cm2'
        conv = (10^12) * arena.scale;
    case 'mmol/arena'
        conv = (10^12) / (arena.x * arena.y);
    case 'mmol/cell'
        conv = 10^12;
    case 'fmol/cell'
        conv = 1;
end
cant = cant * conv;
cantMM = cant / ((10^12) * 0.01 * arena.scale);

% Introduce culture media information in the arena
if isfield(arena,'mediaMetID')
    idx = ismember(rxns, arena.mediaMets);
    if any(idx==1)
        index = cellfun(@(x) find(ismember(arena.mediaMets,rxns(x))), num2cell(find(idx)));
        arena.mediaMM(index) = arena.mediaMM(index) + cantMM(index);
    end
    if any(idx==0)
        arena.mediaMetID = [arena.mediaMetID; (arena.mediaMetID(end) + 1 :arena.mediaMetID(end) + sum(~idx))'];
        arena.mediaMets = [arena.mediaMets; rxns(~idx)];
        arena.mediaMM = [arena.mediaMM; cantMM(~idx)];
    end
else
    arena.mediaMetID = (1:length(rxns))';
    arena.mediaMets = rxns;
    arena.mediaMM = cantMM;
end

% Generate substrate ID
exRxns = unique(cat(1,arena.orgExch{:}),'stable');
n = length(exRxns);
exRxnID = (1:n)';

% Generate substrate metabolite name
exRxnsMet = unique(cat(1,arena.orgExchMets{:}),'stable');

% Check available substrates in the arena
idx = ismember(rxns, exRxns);
if sum(idx) == 0
    error('None of the input substrates are present in the arena')
end
if sum(idx) ~= length(rxns)
    fprintf('Some of the input substrates are not present in the arena\n')
end
if length(diffSpeed) == length(rxns)
    diffSpeed = diffSpeed(idx);
end
rxns = rxns(idx);
cant = cant(idx);

% Check diffusion speed
if (length(diffSpeed) ~= length(rxns)) && (length(diffSpeed)~=1)
    error('The parameter diffSpeed should be of the same size of rxns or equal to 1')
end
if length(diffSpeed) == 1
    diffSpeed = repmat(diffSpeed, n, 1);
end

% Generate diffusion matrix
diffMat = cell(n,1);
for i = 1:n
    if ismember(exRxns{i},rxns)
        diffMat{i} = repmat(cant(ismember(rxns,exRxns{i})), arena.x, arena.y);
    else
        diffMat{i} = zeros(arena.x, arena.y);
    end
end

% Generate matrix of diffusion speed
Dgrid = cell(n,1);
for i = 1:n
    Dgrid{i,1} = struct();
    Dgrid{i,1}.xMid = repmat(diffSpeed(i), arena.y, arena.x);
    Dgrid{i,1}.yMid = repmat(diffSpeed(i), arena.y, arena.x);
    Dgrid{i,1}.xInt = repmat(diffSpeed(i), arena.y + 1, arena.x);
    Dgrid{i,1}.yInt = repmat(diffSpeed(i), arena.y, arena.x + 1);
end

% Add information to arena
if isfield(arena,'diffMat')
    arena.diffMat = cellfun(@(x,y) x + y, arena.diffMat, diffMat, 'UniformOutput', false);
else
    arena.exRxnID = exRxnID;
    arena.exRxns = exRxns;
    arena.exRxnsMet = exRxnsMet;
    arena.diffMat = diffMat;
    arena.diffSpeed = diffSpeed;
    arena.Dgrid = Dgrid;
end

end