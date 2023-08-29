function arena = addOrganism(arena, model, name, amount, varargin)

% addOrganism.m
%
% Script to add cell organisms to the environment.
%
% INPUTS:
%
%   - arena:    Environment structure.
%   - model:    Genome scale metabolic model of the organism
%               (COBRA format).
%   - name:     Name of the organism.
%   - amount:   Number of cells of the organism.
%
% OPTIONAL INPUTS:
%
%   - exPattern:        Identifier of the exchange reactions. 
%                       Default: 'EX_'.
%   - setAllExInf:      Boolean variable indicating if lower bounds of
%                       exchange reactions should be set to "-Inf".
%                       Default: true.
%   - posX:             X coordinates of the cells to be allocated.
%   - posY:             Y coordinates of the cells to be allocated.
%   - minWeight:        Growth limit at which organism dies. 
%                       Default (E.coli): 0.083 pg.
%   - maxWeight:        Maximal dry weight of a single organism. 
%                       Default (E.coli): 1.172 pg. 
%   - cellWeightMean:   Mean of starting biomass. 
%                       Default (E.coli): 0.489 pg.
%   - cellWeightDev:    Standard deviation of starting biomass. 
%                       Default (E.coli): 0.132 pg.
%   - cellArea:         Surface that the organism occupies. 
%                       Default (E.coli): 4.42 µm^2. 
%   - speed:            Integer value indicating the radius in which the
%                       cell can move. Default: 2.
%   - biomassT0:        Starting biomass of the organism.
%   - deathRate:        Factor by which biomass is reduced if no growth is
%                       possible. Default (E.coli): 0.21 pg.
%   - growType:         Functional type for growth. Default: 'exponential'.
%                       Options: {'linear','exponential'}.
%   - limitGrowth:      Boolean variable indicating if growth should be 
%                       limited by the maximum weight of the organism.
%                       Default: true.
%   - predator:         Cell array with the names of the organisms which
%                       can kill this one. 
%   - lyse:             Boolean variable indicating in the organism
%                       should lyse after death. Default: false.
%   - chem:             Cell array with the list of exchange reactions
%                       which are chemotaxis attractant for the organism.
%
% OUTPUTS:
%
%   - arena: Environment structure with the information of the organism.
%
% EXAMPLE:
%
%   arenaOrg = addOrganism(arena, model_Ecoli, 'Escherichia_coli', 2, ...
%              'posX', [2;5], 'posY', [1;3], 'speed', 3, 'limitGrowth', ...
%              false, 'growType', 'linear', 'lyse', true)
%
% .. Authors: 
%       - Telmo Blasco, 16/12/2022, University of Navarra, TECNUN School of Engineering.

% Manage arguments
parser = inputParser;
addRequired(parser, 'arena', @(x) isstruct(x))
addRequired(parser, 'model', @(x) isstruct(x))
addRequired(parser, 'name', @(x) ischar(x))
addRequired(parser, 'amount', @(x) isnumeric(x) && x>0 && x == round(x))
addParameter(parser, 'exPattern', 'EX_', @(x) ischar(x) || isempty(x))
addParameter(parser, 'setAllExInf', true, @(x) islogical(x))
addParameter(parser, 'posX', [], @(x) isnumeric(x))
addParameter(parser, 'posY', [], @(x) isnumeric(x))
addParameter(parser, 'minWeight', 0.083, @(x) isnumeric(x) && x>0)
addParameter(parser, 'maxWeight', 1.172, @(x) isnumeric(x) && x>0)
addParameter(parser, 'cellWeightMean', 0.489, @(x) isnumeric(x) && x>0)
addParameter(parser, 'cellWeightDev', 0.132, @(x) isnumeric(x) && x>0)
addParameter(parser, 'cellArea', 4.42, @(x) isnumeric(x) && x>0)
addParameter(parser, 'speed', 2, @(x) isnumeric(x) && x>0 && x == round(x))
addParameter(parser, 'biomassT0', [], @(x) isnumeric(x) && x>0)
addParameter(parser, 'deathRate', 0.21, @(x) isnumeric(x) && x>0)
addParameter(parser, 'growType', 'exponential', @(x) ischar(x))
addParameter(parser, 'limitGrowth', true, @(x) islogical(x))
addParameter(parser, 'predator', [], @(x) iscell(x) || isempty(x))
addParameter(parser, 'lyse', false, @(x) islogical(x))
addParameter(parser, 'chem', [], @(x) iscell(x) || isempty(x))

% Extract argument values
parse(parser, arena, model, name, amount, varargin{:});
arena = parser.Results.arena;
model = parser.Results.model;
name = parser.Results.name;
amount = parser.Results.amount;
exPattern = parser.Results.exPattern;
setAllExInf = parser.Results.setAllExInf;
posX = parser.Results.posX;
posY = parser.Results.posY;
minWeight = parser.Results.minWeight;
maxWeight = parser.Results.maxWeight;
cellWeightMean = parser.Results.cellWeightMean;
cellWeightDev = parser.Results.cellWeightDev;
cellArea = parser.Results.cellArea;
speed = parser.Results.speed;
biomassT0 = parser.Results.biomassT0;
deathRate = parser.Results.deathRate;
growType = parser.Results.growType;
limitGrowth = parser.Results.limitGrowth;
predator = parser.Results.predator;
lyse = parser.Results.lyse;
chem = parser.Results.chem;

% Initialize organism data
orgData = zeros(amount, 4);

% Check organism ID
if isfield(arena,'models')
    flag = true;
    if any(ismember(name,arena.orgName))
        orgID = arena.orgID(ismember(arena.orgName,name));
    else
        orgID = length(arena.models) + 1;
    end
else
    flag = false;
    orgID = 1;
end
orgData(:,1) = orgID;

% Check free grid cells
if amount > sum(sum(arena.occupyMat==0))
    error('More individuals than space on the grid')
end

% Check number of cells per grid cell
if round(arena.scale/(cellArea*10^(-8))) < 1
    error('Physical arena size (Lx, Ly) too small. Maximal amount of cells in one grid would be zero')
end 

% Check cell weights
if maxWeight <= minWeight*2
    error('Maximal weight needs to be bigger than two times minimal weight')
end

% Check growth type
if ~any(strcmp(growType,{'linear','exponential'}))
    error('Wrong growType parameter. Options are: {linear, exponential}')
end

% Check positions
if (~isempty(posX) && ~isempty(posY)) && (length(posX) ~= length(posY))
    error('Introduced X and Y coordinates have different length')
elseif (~isempty(posX) && ~isempty(posY)) && (length(posX) ~= amount)
    error('Introduced X or Y coordinates differ in the amount of cells')
end

% Define cell positions
if isempty(posX.*posY)
    idx = randsample(find(arena.occupyMat==0),amount);
    arena.occupyMat(idx) = orgID;
    [x,y] = ind2sub(size(arena.occupyMat),idx);
    orgData(:,3) = x;
    orgData(:,4) = y;
else
    idx = sub2ind(size(arena.occupyMat),posX,posY);
    if any(arena.occupyMat(idx)~=0)
        error('Introduced X and Y coordinates are already occupied by other cells')
    end
    arena.occupyMat(idx) = orgID;
    orgData(:,3) = posX;
    orgData(:,4) = posY;
end

% Define initial biomass
if isempty(biomassT0)
    orgData(:,2) = normrnd(cellWeightMean, cellWeightDev, amount, 1);
else
    orgData(:,2) = repmat(biomassT0, amount, 1);
end

% Define exchange reactions
exRxns = model.rxns(findExcRxns(model));
exRxns = exRxns(~ismember(exRxns,model.rxns(logical(model.c))));
if ~isempty(exPattern)
    exRxns = exRxns(~cellfun(@isempty, regexp(exRxns,exPattern)));
end

% Define exchange reactions metabolites
exMets = cellfun(@(x) model.mets(logical(model.S(:,ismember(model.rxns,x)))), exRxns);

% Set lower bounds of exchange reactions to minimum value
if setAllExInf
    model.lb(ismember(model.rxns,exRxns)) = -1000;
end

% Calculate biomass reactions substrates
bioSubs = [];
if lyse || ~isempty(predator)
    bioS = full(model.S(any(model.S(:, logical(model.c)) < 0, 2),logical(model.c)));
    bioMets = model.mets(any(model.S(:, logical(model.c)) < 0, 2));
    bioMets = cellfun(@strrep, bioMets, repmat({'[c]'},length(bioMets),1), repmat({'[e]'},length(bioMets),1), 'Un', 0);
    
    % Filter those precursos that are not present in the environment
    idx = ismember(bioMets,exMets);
    bioSubs = [bioMets(idx), num2cell(bioS(idx))];
end

% Introduce organism specifications in the model
model.minWeight = minWeight;
model.maxWeight = maxWeight;
model.cellWeightMean = cellWeightMean;
model.cellWeightDev = cellWeightDev;
model.cellArea = cellArea;
model.speed = speed;
model.deathRate = deathRate;
model.growType = growType;
model.limitGrowth = limitGrowth;
model.predator = predator;
model.lyse = lyse;
model.bioSubs = bioSubs;
model.chem = chem;

% Introduce organism information in the arena
if flag
    if ~any(ismember(name,arena.orgName))
        arena.models = [arena.models; {model}];
        arena.orgID = [arena.orgID; orgID];
        arena.orgName = [arena.orgName; {name}];
        arena.orgExch = [arena.orgExch; {exRxns}];
        arena.orgExchMets = [arena.orgExchMets; {exMets}];
    end
    arena.orgData = [arena.orgData; orgData];
else
    arena.models = {model};
    arena.orgID = orgID;
    arena.orgName = {name};
    arena.orgExch = {exRxns};
    arena.orgExchMets = {exMets};
    arena.orgData = orgData;  
end
    
end