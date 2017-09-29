function lattice = fi_prettyLattice( sortinginds, FM, showFMvals, showStats, showPaths )
% Copyright 2017 Derek T. Anderson and Muhammad Islam
% Mississippi State University
% 7/5/2017
%	We would like to acknowledge others who have contributed in one form or another to this library (code, theory, etc.)
%	Timothy C. Havens
%	Anthony Pinar
%	Christian Wagner
%	James M. Keller
%	Melissa Anderson
%	Daniel Wescott
%	Chee Seng Chan
%	Lequn Hu
%	Ryan Smith
%	Charlie Veal
%	Alina Zare
%	Fred Petry
%	Paul Elmore
%
% This file is part of FuzzyIntegralComputerVisionLibrary.
%     FuzzyIntegralComputerVisionLibrary is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     FuzzyIntegralComputerVisionLibrary is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with FuzzyIntegralComputerVisionLibrary.  If not, see <http://www.gnu.org/licenses/>.
%
% If you use this library, please reference it as such
%
%    @misc{ FuzzyIntegralComputerVisionLibrary,
%		author = "Derek Anderson and James Keller and Chee Seng Chan",
%		title = "Fuzzy Integral and Computer Vision Library",
%		url = {http://derektanderson.com/FuzzyLibrary/},
%		urldate = {2017-07-03}
%		year = "2017" };
%
% Program use:
%   sortinginds - this variable is an output of the QPmatrices function
%   FM          - the binary fuzzy measure found via
%               quadprog in the DeFIMKL algorithm
%   showFMvlas  - a boolean that, if true, shows the FM values on the
%               visualization (default = 0)
%   showStats   - a boolean that, if true, shows how many nodes are used on
%               each layer of the visualization (default = 0)

% these will adjust the node sizes
maxFMarea = 1000; % 200 is default, 2000 will make them large enough to put labels inside
minFMarea = 0.1; %ish

if nargin < 5
    showPaths = 1;
end

if nargin < 4
    showStats = 0;
end

if nargin < 3
    showFMvals = 0;
end

si = sortinginds;
% si( 10:end, : ) = [];

lattice = getpaths( si, FM );
%%

lattice = makeLattice(lattice, 0.3);
x = lattice.newx;
y = lattice.yy;
% plot paths
if showPaths
    for aa = 1:size(si,1)
        p1 = plot(x,y(aa,:),'w','LineWidth',1.5);
        hold on
        p1.Color(4) = 0.25;
    %     p1.Color(4) = 1;
    end
end

% plot nodes
for ii = 1:length(lattice.yloc)
%     sprintf('Layer %d of %d',ii,length(lattice.yloc))
    ya = lattice.yloc{ii};
    xa = lattice.x(ii).*ones( size( ya ) );
    point_area = lattice.FM{ ii }.^2;
    scaled_areas = maxFMarea.*(point_area+minFMarea);
    scatter( xa, ya, scaled_areas, 'black','fill'), hold on
    scatter( xa, ya, scaled_areas, 'white'), hold on
    if showStats
        text(xa( 1 ), -1.2, lattice.stats{ ii }, 'Color', 'white')
    end
    if showFMvals
        %
        for aa = 1:length(ya)
            if mod(aa,2)==0
                lvl = 0.01;
            else
                lvl = -0.01;
            end
            text(xa(1), ya(aa)-sqrt(scaled_areas(aa)/pi)/100', num2str(lattice.FM{ii}(aa)',2), 'Color', 'white')
        end
    end
end
set(gca,'Color',[0 0 0]);
xlim([-1.2 1.2])
ylim([-1.55 1.55])
camroll(90)
set(gca,'xtick',[],'ytick',[])
end

function [lattice] = makeLattice(lattice, scaleFactor)

ndensities = lattice.ndensities;
levelSize = lattice.levelSize;
yinds = lattice.yinds;
x = lattice.x;
FM = lattice.FM;

yloc{ 1 } = 0;
yloc{ ndensities + 1 } = 0;
stats{ 1 } = '';
stats{ ndensities + 1 } = '';

for ii = 2:ndensities
    thislevelsize = levelSize( ii - 1 );
    scaling = ( levelSize( ii - 1 ) / max( levelSize ) ) ^ scaleFactor;
    yloc{ ii } = -linspace( -scaling, scaling, levelSize( ii - 1 ) );
    
    nnodes_used( ii ) = length( unique( lattice.yinds( :, ii ) ) );
    node_inds_unused{ ii } = ~ismember(1:thislevelsize, unique(lattice.yinds(:,ii)));
    stats{ ii } = [num2str( nnodes_used( ii ) ) '/' num2str( levelSize( ii - 1 ) )];
end
stats{ 1 } = [num2str( sum( nnodes_used ) + 2 ) '/' num2str( sum( lattice.levelSize ) + 2 )];
stats{ ndensities + 1 } = [num2str( 100 * (sum( nnodes_used ) + 2) / (sum( lattice.levelSize ) + 2) ) '%'];


ys = zeros( size(yinds) );
for ii = 1:ndensities + 1
    ys( :, ii ) = yloc{ ii }( yinds( :, ii ) );
end

lattice.newx = linspace(-1,1,1000);
lattice.yy = spline(x,ys,lattice.newx);
lattice.ys = ys;
lattice.yloc = yloc;
lattice.stats = stats;
lattice.unusednodes = node_inds_unused;
end

function [lattice] = getpaths( sortingInds, FM )
% sortingInds is a matrix where each row represents an instance, which is
% the singleton ordering
si = sortingInds;

[nsamples, ndensities] = size( si );

for ii = 1:ndensities-1
    levelSize( ii ) = choose( ndensities, ii );
end

y = linspace( -1.5, 1.5, 2*max( levelSize ) );
x = linspace( -1, 1, ndensities + 1 );

for ii = 1:nsamples
%     sprintf('Paths are %0.2f percent complete',100*ii/nsamples)
    ordering = si( ii, : );
    yinds( :, ii ) = getys( ordering, ndensities );
end
yinds = [ones( 1, nsamples ); yinds; ones( 1, nsamples )];
yinds = fliplr( yinds )';

lattice.ndensities = ndensities;
lattice.levelSize = levelSize;
lattice.yinds = yinds;
lattice.x = x;
lattice.FM = FM_lattice( ndensities, FM );
lattice.FM_old = FM;
lattice.si = si;

end

function ys = getys( o, ndensities )
ys = zeros( ndensities - 1, 1 );
ys( 1 ) = o( 1 );
for ii = 2:ndensities-1
    ys( ii ) = getYind( ndensities, ii, o( 1:ii ) );
end
end

function yind = getYind( ndensities, layerNumber, targetArray )
layer = flipud( combnk( 1:ndensities, layerNumber ) );
if layer(1) ~= 1
    layer = flipud( layer );
end
   yind = find( ismember(layer, perms(targetArray), 'rows') == 1 );
end

function newFM = FM_lattice( ndensities, FM )

newFM{ 1 } = 0;
newFM{ ndensities + 1 } = FM( end );

for layer = 1:ndensities - 1
%     sprintf('FMs are %0.2f percent complete',100*layer/(ndensities-1))
    combos = nchoosek( 1:ndensities, layer );
    inds = cumsum(2.^( combos-1 ), 2); % had a + 1 here
    newFM{ layer + 1 } = FM( inds( :, end ) );
end
end

function y=choose(n,k)
% Computes the binomial coefficient "n choose k"
% y=factorial(n)/(factorial(k)*factorial(n-k))
% y = n!/(k!-(n-k)!)

y=factorial(n)/(factorial(k)*factorial(n-k));
end