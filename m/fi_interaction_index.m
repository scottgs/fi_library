function [ InteractionIndex ] = fi_interaction_index( FMbinary )
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
%	Inputs:
%		FMbinary = fuzzy measure
%	Outputs:
%		InteractionIndex = interaction index terms
%
% Program description:
%   compute the interaction index

	% number of FM variables
	nVars = length(FMbinary);

	% append FM with  value for empty set 
	FMbinary = [0 ; FMbinary(:)];
	% number of inputs
	N = ceil(log2(nVars));

	InteractionIndex = nan(N);

	% binary representation of variable Index
	vars = de2bi(1:nVars, N);

	% layer index of variables
	varsLayers = sum(vars, 2);
	for i = 1:N
	  for j=i+1:N
		% set of variables that include input i
		varsWithIJ = find(all(vars(:,[i j]),2));
		
		varsWithoutI = varsWithIJ - 2^(i-1);
		varsWithoutJ = varsWithIJ - 2^(j-1);
		varsWithoutIJ = varsWithIJ - 2^(i-1) - 2^(j-1);
		
		% Cardinal of sets without current input
		K = sum(vars(varsWithIJ,:),2) - 2;
		
		zetaK = factorial(N - K - 2).*factorial(K)/factorial(N - 1);
		
		InteractionIndex(i,j) = sum(zetaK.*(FMbinary(varsWithIJ+1) - FMbinary(varsWithoutI+1) - FMbinary(varsWithoutJ+1) + FMbinary(varsWithoutIJ+1)));
		end
	end
	InteractionIndex = triu(InteractionIndex) + triu(InteractionIndex,1)';

	InteractionIndex(isnan(InteractionIndex)) = 1;
	
end