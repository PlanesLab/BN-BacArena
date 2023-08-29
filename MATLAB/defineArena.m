function arena = defineArena(x,y,varargin)

% defineArena.m
%
% Script to define an environment.
%
% INPUTS:
%
%   - x:  Horizontal size of the environment.
%   - y:  Vertical size of the environment.
%
% OPTIONAL INPUTS:
%
%   - Lx:     Horizontal grid size in cm.
%   - Ly:     Vertical grid size in cm.
%   - seed:   Random number seed to be reproducible.
%
% OUTPUTS:
%
%   - arena: Environment structure.
%
% EXAMPLE:
%
%   arena = defineArena(20, 20, 'Lx', 0.01, 'Ly', 0.01, 'seed', 1200)
%
% .. Authors: 
%       - Telmo Blasco, 08/06/2023, University of Navarra, TECNUN School of Engineering.

% Manage arguments
parser = inputParser;
addRequired(parser, 'x', @(x) isnumeric(x) && x>0 && x == round(x))
addRequired(parser, 'y', @(x) isnumeric(x) && x>0 && x == round(x))
addParameter(parser, 'Lx', [], @(x) isempty(x) || (isnumeric(x) && x>0))
addParameter(parser, 'Ly', [], @(x) isempty(x) || (isnumeric(x) && x>0))
addParameter(parser, 'seed', randsample(1:10000,1), @(x) isnumeric(x) && x>0 && x == round(x))

% Extract argument values
parse(parser, x, y, varargin{:});
x = parser.Results.x;
y = parser.Results.y;
Lx = parser.Results.Lx;
Ly = parser.Results.Ly;
seed = parser.Results.seed;

% Set seed
rng(seed)

% Define cell grid size
if isempty(Lx)
    Lx = (0.025/100) * x;
end
if isempty(Ly)
    Ly = (0.025/100) * y;
end

% Define grid geometry
gridGeometry = struct();
gridGeometry.xUp = 0;
gridGeometry.xDown = Ly;
gridGeometry.xMid = (Ly/y/2):(Ly/y):(Ly-(Ly/y/2));
gridGeometry.xInt = (0:(Ly/y):Ly)';
gridGeometry.dx = repmat(Ly/y, x, y);
gridGeometry.dxAux = [repmat(Ly/y/2, 1, x); repmat(Ly/y, y-1, x); repmat(Ly/y/2, 1, x)];
gridGeometry.xN = y;
gridGeometry.yUp = 0;
gridGeometry.yDown = Lx;
gridGeometry.xMid = (Lx/x/2):(Lx/x):(Lx-(Lx/x/2));
gridGeometry.yInt = (0:(Lx/x):Lx)';
gridGeometry.dy = repmat(Lx/x, y, x);
gridGeometry.dyAux = [repmat(Lx/x/2, y, 1), repmat(Lx/x, y, x-1), repmat(Lx/x/2, y, 1)];
gridGeometry.yN = x;

% Introduce information in the arena
arena = struct();
arena.x = x;
arena.y = y;
arena.Lx = Lx;
arena.Ly = Ly;
arena.scale = (Lx*Ly)/(x*y);
arena.occupyMat = sparse(x,y);
arena.gridGeometry = gridGeometry;

end