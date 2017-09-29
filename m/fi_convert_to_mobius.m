function [FMmobbinary] = fi_convert_to_mobius(N, FMbinary)
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
%		N = number of inputs
%		FMbinary = binary coded fuzzy measure
%	Outputs:
%		FMmobbinary = resultant fuzzy measure in mobius coded format
%
% Program description:
%   make mobius fuzzy measure versus binary fuzzy measure

	FMmobbinary = zeros(size(FMbinary));

	FMbinary = FMbinary(:);
	% number of FM variables
	nVars = length(FMbinary);

	% number of inputs
	%N = ceil(log2(nVars));

	if ((nVars ~= 2^N-1))
	  error('Number of FM variables must be 2^N-1');
    end

	% binary representation of variable Index
	vars = de2bi(1:nVars, N);

	% layer index of variables
	varsLayers = sum(vars, 2);
	for i = 1:nVars
		curVar = vars(i,:);
		inputsNotPresent = find(~curVar);
		
		% subset of current variables
		subsetIdx = find(~any(vars(:,inputsNotPresent),2));
		
		% compute the mobius value
		if (mod(sum(curVar),2) == 0)
			temp = (-1).^(varsLayers(subsetIdx));
		else
			temp = (-1).^(varsLayers(subsetIdx)+1);
		end
		FMmobbinary(i) = sum(FMbinary(subsetIdx).*temp);
	end
	FMmobbinary = FMmobbinary(:)';

end	