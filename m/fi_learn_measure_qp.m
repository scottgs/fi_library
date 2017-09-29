function [FMbinary] = fi_learn_measure_qp( data_set_with_labels )
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
%		data_set_with_labels = matrix of size (I,N+1), where I is number of instances with N values and N+1th element is "label"
%	Outputs:
%		FMbinary = resultant fuzzy measure in binary coded format
%
% Program description:
%   learn the fuzzy measure and ChI w.r.t. quadratic programming 

	% just run with no regularization
	FMbinary = fi_learn_measure_qp_reg( data_set_with_labels , 0 );
	
end