function [res] = fi_choquet_integral_h_and_g_form( inputs , FM )
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
%		inputs = e.g., [ 0.2 0.3 0.4 ]
%		FM = binary, e.g., [ x1 x2 x12 x3 x13 x23 x123 ]
%	Outputs:
%		res = result of ChI
%
% Program description:
%   Compute the discrete (finite X) Choquet integral
%		this is the "h's and g's form/representation"
%	Assumes the FM is in "binary coded" format
	
    [SortVal, SortInd] = sort( inputs , 2, 'descend' ); % sort
	
	SortVal = [SortVal zeros(size(inputs,1),1)]; % append 0 for difference form below
	
	% compute the integral
	i = cumsum(2.^(SortInd-1),2);
    res = sum(FM(i).*(SortVal(:,1:end-1)-SortVal(:,2:end)),2);

end