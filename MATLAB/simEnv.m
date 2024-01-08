function arena = simEnv(arena, time, varargin)

% simEnv.m
%
% Script to simulate the environment evolution.
%
% INPUTS:
%
%   - arena:  Environment structure.
%   - time:   Number of time steps.
%
% OPTIONAL INPUTS:
%
%   - cutOff:     Threshold for numeric accuracy.
%   - secObj:     Secondary objective for bi-level LP. 
%                 Options: {'none','pFBA'}.
%   - bacCoeff:   Cell type X Nutrient coefficient matrix.
%   - nutCoeff:   Cell type X Nutrient coefficient matrix.
%   - verbose:    Logical variable to print status messages.
%
% OUTPUTS:
%
%   - arena: Updated environment structure.
%
% EXAMPLE:
%
%   arena = simEnv(arena, 6, 'secObj', 'none', 'bacCoeff', bacMat, ...
%           'nutCoeff', nutMat, 'verbose', false)
%
% .. Authors: 
%       - Telmo Blasco, 14/11/2023, University of Navarra, TECNUN School of Engineering.

% Manage arguments
parser = inputParser;
addRequired(parser, 'arena', @(x) isstruct(x))
addRequired(parser, 'time', @(x) isnumeric(x) && x>0 && x == round(x))
addParameter(parser, 'cutOff', 1e-06, @(x) isnumeric(x) && x>0)
addParameter(parser, 'secObj', 'none', @(x) ischar(x))
addParameter(parser, 'bacCoeff', [], @(x) isnumeric(x))
addParameter(parser, 'nutCoeff', [], @(x) isnumeric(x))
addParameter(parser, 'verbose', true, @(x) islogical(x))

% Extract argument values
parse(parser, arena, time, varargin{:});
arena = parser.Results.arena;
time = parser.Results.time;
cutOff = parser.Results.cutOff;
bacCoeff = parser.Results.bacCoeff;
nutCoeff = parser.Results.nutCoeff;
secObj = parser.Results.secObj;
verbose = parser.Results.verbose;

% Check culture media
if ~isfield(arena,'diffMat')
    error('No media present in the arena')
end

% Check secondary objective value
if ~any(strcmp(secObj,{'none','pFBA'}))
    error('Wrong secObj parameter. Options are: {none, pFBA}')
end

% Secondary objective recommendation
if strcmp(secObj,'none') && any(cellfun(@(x) x.limitGrowth, arena.models))
    warning('If growth is limited by maximum weight (limitGrowth = true), parsimonius FBA is recommended (secObj = pBA)')
end

% Check coefficient matrix
if ~isempty(bacCoeff) && (size(bacCoeff,1) ~= length(arena.models) || size(bacCoeff,2) ~= length(arena.models))
    error('Cell coefficient matrix should have the same number of rows and columns as the number of cell types')
else
    arena.bacCoeff = bacCoeff;
end
if ~isempty(nutCoeff) && (size(nutCoeff,1) ~= length(arena.models) || size(nutCoeff,2) ~= length(arena.mediaMets))
    error('Nutrient coefficient matrix should have the same number of rows as the number of cell types and columns as the number of metabolites in the media')
else
    arena.nutCoeff = nutCoeff;
end

% Model parameters checking within arena
if any(~cellfun(@(x) isempty(x.chem), arena.models))
    idx = find(~cellfun(@(x) isempty(x.chem), arena.models));
    for i = 1:length(idx)
        arena.models{idx(i)}.chem = arena.models{idx(i)}.chem(ismember(arena.models{idx(i)}.chem, arena.exRxns));
    end
end

% Build output information
arena.relAbundance = cell(time+1,1);
arena.nutAbundance = cell(time+1,1);
arena.exAbundance = cell(time+1,1);
arena.simList = cell(time+1,1);
arena.mfluxList = cell(time+1,1);

% Initital iteration cell position information
arena.simList{1,1} = arena.orgData;

% Initital iteration cell information
arena.relAbundance{1,1} = [arena.orgID, zeros(length(arena.orgID),1)];
ab = tabulate(arena.orgData(:,1));
arena.relAbundance{1,1}(cellfun(@(x) find(ismember(arena.orgID,x)), num2cell(ab(:,1))),2) = ab(:,3);

% Initital iteration nutrient information
arena.nutAbundance{1,1} = [arena.mediaMetID, arena.mediaMM];
arena.exAbundance{1,1} = [arena.exRxnID, cellfun(@(x) sum(sum(x))/((10^12) * 0.01 * arena.scale * arena.x * arena.y), arena.diffMat)];

% Initial iteration exchanges information
arena.mfluxList{1,1} = cellfun(@(x) zeros(length(x),1), arena.orgExch, 'Un', 0);

% Simulate through time
for i = 1:time
    
    % Number of organisms
    n = size(arena.orgData,1);
    
    % Shuffle organism information
    arena.orgData = arena.orgData(randsample(1:n, n),:);
    
    % Define exchanges information
    arena.mflux = cellfun(@(x) zeros(length(x),1), arena.orgExch, 'Un', 0);

    % Simulate through organism
    if verbose
        fprintf('##########################################################\n')
        sentence = fprintf('Organisms 0 / %d', n);
    end
    tic
    for j = 1:n
        if verbose
            fprintf(repmat('\b',1,sentence))
            sentence = fprintf('Organisms %d / %d\n', j, n);
        end
        arena = simBac(arena, j, i, 'cutOff', cutOff, 'secObj', secObj);
    end
    simTime = toc;
    
    % Remove dead cells
    idx = arena.orgData(:,2) == -Inf;
    if any(idx)
        arena.occupyMat(sub2ind(size(arena.occupyMat),arena.orgData(idx,3),arena.orgData(idx,4))) = 0;
        arena.orgData(idx,:) = [];
    end
    
    % Extract unvariant exchange metabolites
    idx = (cellfun(@(x) length(unique(x)), arena.diffMat) ~= 1);
    nM = length(idx);
    
    % Perform diffusion
    if verbose
        sentence = fprintf('Substances 0 / %d', length(idx));
    end
    tic
    for j = 1:nM
        if verbose
            fprintf(repmat('\b',1,sentence))
            sentence = fprintf('Substances %d / %d\n', j, length(idx));
        end
        
        if idx(i)
            arena.diffMat{j,1} = diffuseR(arena.diffMat{j,1});
        end
    end
    diffTime = toc;
    
    % Extract nutrient and exchange abundances
    arena.exAbundance{i+1,1} = [arena.exRxnID, cellfun(@(x) sum(sum(x))/((10^12) * 0.01 * arena.scale * arena.x * arena.y), arena.diffMat)];
    for j = 1:arena.mediaMetID(end)
        idx = ismember(arena.exRxns, arena.mediaMets(j));
        if any(idx)
            arena.mediaMM(j) = arena.exAbundance{i+1,1}(idx,2);
        end
    end
    arena.nutAbundance{i+1,1} = [arena.mediaMetID, arena.mediaMM];
    
    % Update abundance variables
    arena.relAbundance{i+1,1} = [arena.orgID, zeros(length(arena.orgID),1)];
    ab = tabulate(arena.orgData(:,1));
    arena.relAbundance{i+1,1}(cellfun(@(x) find(ismember(arena.orgID,x)), num2cell(ab(:,1))),2) = ab(:,3);
    
    % Update cell position information
    arena.simList{i+1,1} = arena.orgData;
    
    % Update exchanges information
    arena.mfluxList{i+1,1} = arena.mflux;
    
    % Show progress
    if verbose
        fprintf('Iteration: %d\t Organisms: %d\t Biomass: %.3f pg\n',i,size(arena.orgData,1),round(sum(arena.orgData(:,2)),3))
        fprintf('----------------------------------------------------------\n')
        if size(arena.orgData,1) > 0
            tmp = tabulate(arena.orgData(:,1));
            bio = cellfun(@(x) sum(arena.orgData(ismember(arena.orgData(:,1),x),2)), num2cell(tmp(:,1)));
            for k = 1:size(tmp,1)
                fprintf('Organism: %s\t Count: %d\t Biomass: %.3f\n',arena.orgName{k,1},tmp(k,2),bio(k))
            end
        end
        fprintf('----------------------------------------------------------\n')
        fprintf('Time total: %.3f\t Diffusion: %.3f\t (%.3f%%)\n',simTime+diffTime, diffTime, (diffTime/(diffTime+simTime))*100)
        fprintf('----------------------------------------------------------\n')
    end
end

end

function arena = simBac(arena, j, time, varargin)

% simBac.m
%
% Script to simulate the interactions of cells at each time step.
%
% INPUTS:
%
%   - arena:  Environment structure.
%   - j:      Cell number.
%   - time:   Time step.
%
% OPTIONAL INPUTS:
%
%   - cutOff:   Threshold for numeric accuracy.
%   - secObj:   Secondary objective for bi-level LP. 
%               Options: {'none','pFBA'}.
%
% OUTPUTS:
%
%   - arena: Updated environment structure.
%
% EXAMPLE:
%
%   arena = simBac(arena, 1, 10, 'cutOff', 1e-08, 'secObj', 'pFBA')
%
% .. Authors: 
%       - Telmo Blasco, 14/11/2023, University of Navarra, TECNUN School of Engineering.

% Manage arguments
parser = inputParser;
addRequired(parser, 'arena', @(x) isstruct(x))
addRequired(parser, 'j', @(x) isnumeric(x) && x>0 && x == round(x))
addRequired(parser, 'time', @(x) isnumeric(x) && x>0 && x == round(x))
addParameter(parser, 'cutOff', 1e-06, @(x) isnumeric(x) && x>0)
addParameter(parser, 'secObj', 'none', @(x) ischar(x))

% Extract argument values
parse(parser, arena, j, time, varargin{:});
arena = parser.Results.arena;
j = parser.Results.j;
time = parser.Results.time;
cutOff = parser.Results.cutOff;
secObj = parser.Results.secObj;

% Define organism ID, model and exchange reactions
orgID = arena.orgData(j,1);
model = arena.models{orgID};
exRxns = arena.orgExch{orgID};
n = length(exRxns);

% Extract organism location
posX = arena.orgData(j,3);
posY = arena.orgData(j,4);

% Check predators in the Moore neighborhoud
if ~isempty(model.predator)
    nHoods = checkHood(arena, posX, posY, arena.occupyMat, 'step', 1, 'presence', true);
    species = arena.orgName(arena.occupyMat(sub2ind(size(arena.occupyMat),nHoods(:,1), nHoods(:,2))));
    
    if any(ismember(species,model.predator))
        % Cell lysis
        bioSubs = model.bioSubs(ismember(model.bioSubs(:,1),arena.exRxnsMet),:);
        idx = cellfun(@(x) find(ismember(arena.exRxnsMet,x)), bioSubs(:,1));
        for i = 1:length(idx)
            arena.diffMat{idx(i),1}(posX,posY) = arena.diffMat{idx(i),1}(posX,posY) + round(abs(bioSubs{i,2})*model.minWeight,6);
        end
        
        % Cell death
        arena.orgData(j,2) = -Inf;
        return
    end
end

% Calculate growth factor
factor = 1;
if ~isempty(arena.bacCoeff) && ~isempty(arena.nutCoeff)
    factor = (sum(arena.relAbundance{time,1}(:,2).*full(arena.bacCoeff(orgID,:))') + sum(arena.nutAbundance{time,1}(:,2).*full(arena.nutCoeff(orgID,:))') + arena.relAbundance{time,1}(orgID,2)) / arena.relAbundance{time,1}(orgID,2);
elseif ~isempty(arena.bacCoeff) && isempty(arena.nutCoeff)
    factor = (sum(arena.relAbundance{time,1}(:,2).*full(arena.bacCoeff(orgID,:))') + arena.relAbundance{time,1}(orgID,2)) / arena.relAbundance{time,1}(orgID,2);
elseif isempty(arena.bacCoeff) && ~isempty(arena.nutCoeff)
    factor = (sum(arena.nutAbundance{time,1}(:,2).*full(arena.nutCoeff(orgID,:))') + arena.relAbundance{time,1}(orgID,2)) / arena.relAbundance{time,1}(orgID,2);
end
if factor < 0
    factor = 0;
end

% Set reaction indexes
idx = cellfun(@(x) find(ismember(arena.exRxns,x)), exRxns);
exInd = cellfun(@(x) find(ismember(model.rxns,x)), exRxns);

% Save original bounds
lobnd = model.lb;

% Extract culture media
lb = -cellfun(@(x) arena.diffMat{x,1}(posX,posY), num2cell(idx));
lb(lb>0) = 0;

% Costrain according to flux definition: mmol/(gDW*hr)
lobnd(exInd) = model.lb(exInd) * (arena.orgData(j,2) / model.cellWeightMean) * time;
if any(lb < lobnd(exInd))
    ind = find(lb < lobnd(exInd));
    lb(ind) = lobnd(exInd(ind));
end

% Set culture media
model = changeRxnBounds(model, exRxns, lb, 'l');

% Limit growth
if model.limitGrowth
    growthLimit = (model.maxWeight * 1.5) - arena.orgData(j,2);
    if growthLimit > 0
        model.ub(logical(model.c)) = growthLimit * factor;
    else
        model.ub(logical(model.c)) = cutOff;
    end
end

% Perform FBA
fbasol = optimizeCbModel(model, 'max');
fluxes = fbasol.v(exInd);

% Secondary objective
if fbasol.f > cutOff
    switch secObj

        % Parsimonius FBA
        case 'pFBA'
            [~, ~, tmpM, tmpSol] = pFBA(model, 'geneoption', 0, 'skipclass', -1);
            tmpSol.v(~cellfun(@isempty,regexp(tmpM.rxns, '_b$'))) = - tmpSol.v(~cellfun(@isempty,regexp(tmpM.rxns, '_b$')));
            tmpM.rxns = regexprep(regexprep(tmpM.rxns, '_b$', ''), '_f$', '');
            tmpM.rxns(abs(tmpSol.v) < cutOff) = [];
            tmpSol.v(abs(tmpSol.v) < cutOff) = [];
            fluxes = cellfun(@(x) tmpSol.v(ismember(tmpM.rxns,x)), exRxns, 'Un', 0);
            fluxes(cellfun(@isempty,fluxes)) = {0};
            fluxes = cell2mat(fluxes);
    end
end

% Update cell grid substrates concentration
if fbasol.f > cutOff
    for i = 1:n
        arena.diffMat{idx(i),1}(posX,posY) = round(arena.diffMat{idx(i),1}(posX,posY) + fluxes(i), round(-log10(cutOff)));
    end
end

% Update exchanges fluxes information
arena.mflux{orgID} = arena.mflux{orgID} + fluxes; 

% Check growth
switch model.growType
    case 'linear'
        if fbasol.f > 0
            arena.orgData(j,2) = arena.orgData(j,2) + fbasol.f;
        else
            arena.orgData(j,2) = arena.orgData(j,2) - (model.deathRate * time);
        end
    case 'exponential'
        if fbasol.f > 0
            arena.orgData(j,2) = (arena.orgData(j,2) * fbasol.f) + arena.orgData(j,2);
        else
            arena.orgData(j,2) = arena.orgData(j,2) - (arena.orgData(j,2) * model.deathRate * time);
        end
end

% Cell division
while arena.orgData(j,2) > model.maxWeight
    % Check free grid cells
    nHoods = checkHood(arena, posX, posY, arena.occupyMat, 'step', 1, 'presence', false);
    if isempty(nHoods)
        break
    else
        % Select random position
        pos = nHoods(randsample(1:size(nHoods,1),1),:);
        
        % Divide original biomass
        arena.orgData(j,2) = arena.orgData(j,2) / 2;
        
        % Generate new daughter cell
        arena.orgData = [arena.orgData; [arena.orgData(j,1), arena.orgData(j,2), pos(1), pos(2)]];
        arena.occupyMat(pos(1), pos(2)) = arena.orgData(j,1);
    end
end

% Cell death
if arena.orgData(j,2) < model.minWeight
    arena.orgData(j,2) = -Inf;
end

% Cell lysis
if (arena.orgData(j,2) == -Inf) && model.lyse
    bioSubs = model.bioSubs(ismember(model.bioSubs(:,1),arena.exRxnsMet),:);
    idx = cellfun(@(x) find(ismember(arena.exRxnsMet,x)), bioSubs(:,1));
    for i = 1:length(idx)
        arena.diffMat{idx(i),1}(posX,posY) = arena.diffMat{idx(i),1}(posX,posY) + round(abs(bioSubs{i,2})*model.minWeight,6);
    end
end

% Cell movement
if (arena.orgData(j,2) ~= -Inf) && (model.speed ~= 0)
    nHoods = checkHood(arena, posX, posY, arena.occupyMat, 'step', model.speed, 'presence', false);
    if ~isempty(nHoods) && isempty(model.chem)
        % Select random position
        pos = nHoods(randsample(1:size(nHoods,1),1),:);

        % Move cell
        arena.occupyMat(arena.orgData(j,3), arena.orgData(j,4)) = 0;
        arena.occupyMat(pos(1), pos(2)) = orgID;
        arena.orgData(j,[3,4]) = pos;
    elseif ~isempty(nHoods) && ~isempty(model.chem) 
        
        % Get metabolite attractant values
        chemIdx = cellfun(@(x) find(ismember(arena.exRxns,x)), model.chem)';
        val = cellfun(@(x) arena.diffMat{x}(sub2ind([arena.x, arena.y], nHoods(:,1), nHoods(:,2))), num2cell(chemIdx), 'Un', 0);
        val = sum(reshape(cat(1,val{:}), size(nHoods,1), length(model.chem)),2);
        
        % Select random position
        nHoods = nHoods(val == max(val),:);
        pos = nHoods(randsample(1:size(nHoods,1),1),:);
        
        % Move cell
        arena.occupyMat(arena.orgData(j,3), arena.orgData(j,4)) = 0;
        arena.occupyMat(pos(1), pos(2)) = orgID;
        arena.orgData(j,[3,4]) = pos;
    end
end

end

function nHoods = checkHood(arena, posX, posY, mat, varargin)

% checkHood.m
%
% Script to explore the Moore neighbourhood of a given position.
%
% INPUTS:
%
%   - arena:  Environment structure.
%   - posX:   X coordinate of the position.
%   - posY:   Y coordinate of the position.
%   - mat:    Matrix with the information of the grid cells.
%
% OPTIONAL INPUTS:
%
%   - step:        Radious of the Moore neighbourhood to explore.
%                  Default = 1.
%   - presence:    If true the function returns positions in the arena
%                  occupied by cells, otherwise, returns free positions.
%                  Default: true.
%
% OUTPUTS:
%
%   - nHoods: Positions of interest in the Moore neighbourhood.
%
% EXAMPLE:
%
%   nHoods = checkHood(arena, 10, 21, arena.occupyMat, 'step', 2, ...
%                      'presence', false)
%
% .. Authors: 
%       - Telmo Blasco, 08/06/2023, University of Navarra, TECNUN School of Engineering.

% Manage arguments
parser = inputParser;
addRequired(parser, 'arena', @(x) isstruct(x))
addRequired(parser, 'posX', @(x) isnumeric(x) && x>0 && x == round(x))
addRequired(parser, 'posY', @(x) isnumeric(x) && x>0 && x == round(x))
addRequired(parser, 'mat', @(x) ismatrix(x))
addParameter(parser, 'step', 1, @(x) isnumeric(x) && x>0 && x == round(x))
addParameter(parser, 'presence', true, @(x) islogical(x))

% Extract argument values
parse(parser, arena, posX, posY, mat, varargin{:});
arena = parser.Results.arena;
posX = parser.Results.posX;
posY = parser.Results.posY;
mat = parser.Results.mat;
step = parser.Results.step;
presence = parser.Results.presence;

% Extract neighborhood index
xp = (posX-step):(posX+step);
yp = (posY-step):(posY+step);

% Remove index outside the arena
xp = xp(logical((xp>0) .* (xp<=arena.x)));
yp = yp(logical((yp>0) .* (yp<=arena.y)));

% Extract all combinations
[M,N] = ndgrid(xp,yp);
nHoods = [M(:),N(:)];

% Remove original positions
nHoods(logical((nHoods(:,1)==posX) .* (nHoods(:,2)==posY)),:) = [];

% Extract occupied or free positions
if presence
    ind = full(mat(sub2ind(size(mat),nHoods(:,1),nHoods(:,2))) ~= 0);
else
    ind = full(mat(sub2ind(size(mat),nHoods(:,1),nHoods(:,2))) == 0);
end
nHoods = nHoods(ind,:);

end

function mat = diffuseR(mat)

% diffuseR.m
%
% Script to perform the diffusion of a compound in the Moore neighbourhood.
%
% INPUTS:
%
%   - mat:  Diffusion matrix of a compound.
%
% OUTPUTS:
%
%   - mat:  Updated diffusion matrix of a compound.
%
% EXAMPLE:
%
%   diffMatPost = diffuseR(diffMatPrev)
%
% .. Authors: 
%       - Telmo Blasco, 08/06/2023, University of Navarra, TECNUN School of Engineering.

% Size of diffusion matrix
[nr, nc] = size(mat);
smat = mat;

% Check by Moore neighbourhood
for nX = 1:nr
    for nY = 1:nc
        
        % Define neighbours
        nHoods = [[nX,nY];[nX,nY-1];[nX,nY+1];[nX-1,nY];[nX-1,nY-1];[nX-1,nY+1];[nX+1,nY];[nX+1,nY-1];[nX+1,nY+1]];
        nHoods = nHoods((nHoods(:,1) > 0 & nHoods(:,1) <= nr) & (nHoods(:,2) > 0 & nHoods(:,2) <= nc),:);
        
        % Define minimum value
        neighbours = smat(sub2ind([nr,nc], nHoods(:,1), nHoods(:,2)));
        minC = min(neighbours);
        
        % Update values
        if smat(nX,nY) > minC
            nMin = find(neighbours == minC);
            if length(nMin) ~= 1
                nMin = nMin(randsample(length(nMin), 1));
                nHoods = nHoods(nMin,:);
            end
            meanVal = mean([mat(nX,nY), mat(sub2ind([nr,nc], nHoods(1,1), nHoods(1,2)))]);
            mat(nX,nY) = meanVal;
            mat(sub2ind([nr,nc],nHoods(1,1),nHoods(1,2))) = meanVal;
        end
    end
end

end