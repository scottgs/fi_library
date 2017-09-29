function [FMbinary] = fi_convert_to_standard_FM(N, FMmobius)
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
%		FMmobius = fuzzy measure in mobius form
%	Outputs:
%		FMbinary = resultant fuzzy measure in binary coded format
%
% Program description:
%   convert from mobius to binary coded fuzzy measure format

	nVars = length(FMmobius);

	varsBin = de2bi(1:nVars, N);

	if ((nVars != 2^N-1))
	  error("Number of FM variables must be 2^N-1");
	endif

	FMbinary = zeros(nVars, 1);
	for i=1:nVars
		curVarBin = find(~varsBin(i,:));
		idxVarsBin = find(~any(varsBin(:,curVarBin),2));
		if (~isempty(idxVarsBin))
			FMbinary(i) = sum(FMmobius(idxVarsBin));
		else
			FMbinary(i) = sum(FMmobius);
		end
	end
	FMbinary = reshape(FMbinary, size(FMmobius));
	
end