function [ ShapleyVector ] = fi_shapley( FMbinary )
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
%		FMbinary = binary coded fuzzy measure
%	Outputs:
%		ShapleyVector = shapley vector/values
%
% Program description:
%   Compute the N Shapley values for a given fuzzy measure

	% number of FM variables
	nVars = length(FMbinary);

	% append FM with  value for empty set 
	FMbinary = [ 0 ; FMbinary(:) ];
	% number of inputs
	N = ceil(log2(nVars));

	ShapleyVector = zeros(1, N);

	% binary representation of variable Index
	vars = de2bi(1:nVars, N);

	for i = 1:N
		% set of variables that include input i
		varsWithCurInput = find(vars(:,i));
		varsWithoutCurIput = varsWithCurInput - 2^(i-1);
		% cardinal of sets without current input
		K = sum(vars(varsWithCurInput,:),2) - 1;
		gammaK = factorial(N - K - 1).*factorial(K)/factorial(N);
		ShapleyVector(i) = sum(gammaK.*(FMbinary(varsWithCurInput+1) - FMbinary(varsWithoutCurIput+1)));
	end
	
end