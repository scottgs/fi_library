function [FMbinaryResult] = fi_set_value_for_specific_variable( query, FMbinary , value )
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
%		query = binary index, e.g., for set {1,3} you pass [1 3]
%		FMbinary = binary, e.g., [ x1 x2 x12 x3 x13 x23 x123 ]
%		value = what value to store?
%	Outputs:
%		FMbinaryResult = result of storage
%	Example:
%		[FMbinaryResult] = fi_set_value_for_specific_variable( [1 3], FMbinary , 0.9 )
%
% Program description:
%   Set a specific variable in the FM 

	FMbinaryResult = FMbinary;

	index_terms = cumsum(2.^[query-1]); % get all indices in their order provided
	index_terms = index_terms(end); % take just the last one for complete tuple

	FMbinaryResult(index_terms) = value;	

end