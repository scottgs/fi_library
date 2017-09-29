function [FiResult] = fi_choquet_integral_k_additive( inputs , FMmobius )
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
%		inputs = DxN, where D is # of instances and N is # of inputs (so can batch process)
%		FMmobius = mobius fuzzy measure
%	Outputs:
%		FiResult = result of k-additive fuzzy integral
%
% Program description:
%   k-additive fuzzy integral

	[M, N] = size(inputs);

	% find the additive order, k
	nVars = length(FMmobius);

	vars = [];
	i = 1;
	while(length(vars) < nVars)
	  varsKLayer = sum(2.^(nchoosek(1:N, i)-1),2);
	  vars = [vars; varsKLayer(:)];
	  i = i+1;
	end

	if (nVars ~= length(vars))
	  error('number of Mobius terms does not with any of the k-additive orders');
	end

	vars = sort(vars,'ascend');
	varsBinary = de2bi(vars, N);

	coeff = cell2mat(arrayfun(@(x) min(inputs(:,logical(varsBinary(x,:))),[],2),1:nVars,'Un',0));

	FiResult = coeff*FMmobius(:);

end