function [nodes, edge_source, edge_target, labelsOut] = FM_to_Lattice(nVar, input)
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

	nodes = num2str(1:nVar);
	nInput = log2(nVar+1);

	C = monotonicity_constr_lattice(nInput);

	edge_source = num2str(cell2mat(arrayfun (@(x) find(C(x,:)==1), 1:size(C,1),'Un',0)));
	edge_target = num2str(cell2mat(arrayfun (@(x) find(C(x,:)==-1), 1:size(C,1),'Un',0)));
	%[~, edge_source] =ind2sub(size(C), C==1);

	%[~, edge_target] = find(C == -1);

	labels = binaryPos2LexPos(1:nVar, nInput);

	labelsOut = '';
	if (nargin == 2)
	  for i=1:length(labels)
		temp = ['g' num2str(labels(i)) '=' num2str(input(i))];
		labelsOut = [labelsOut ' ' temp];
	  end
	  else
		for i=1:length(labels)
		temp = ['x' num2str(labels(i))];
		labelsOut = [labelsOut ' ' temp];
	  end
	end
	
end