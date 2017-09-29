function [ FIresult ] = fi_ndfi( inputs , FMbinary , whichintegral , plotit )
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
%		inputs = discrete domain due to ndfi, e.g., [ 0.2 0.3 0.1 ; 0.6 0.1 1 ] is two fuzzy sets on domain X such that set A is 0.2/a 0.3/b 0.1/c and set B is 0.6/a 0.1/b and 1/c
%		FMbinary = binary, e.g., [ x1 x2 x12 x3 x13 x23 x123 ]
%		whichintegral = 0 is Sugeno integral and 1 is Choquet integral based
%		plotit = 1 for if you want to plot the result
%	Outputs:
%		FIresult = result of NDFI
%
% Program description:
%   Compute the discrete (finite X) NDFI for set-valued data
%	Assumes the FM is in "binary coded" format

	% how many inputs and what is the number (N) of discrete elements (D)
	[N D] = size( inputs );

	% compute the FI at each element in D
	FIresult = zeros( 1 , D );
	for i=1:D
		if( whichintegral == 0 ) %sugeno
			FIresult(i) = fi_sugeno_integral_h_and_g_form( inputs(:,i)' , FMbinary );		
		else %choquet
			FIresult(i) = fi_choquet_integral_h_and_g_form( inputs(:,i)' , FMbinary );
		end
	end
	
	% plot it?
	if( plotit == 1 )
		figure; 
		hold on; 
		area( FIresult );
		termset = {};   
		if( whichintegral == 0 ) %sugeno
			termset{1} = sprintf( 'sugeno result', i );
		else
			termset{1} = sprintf( 'choquet result', i );		
		end  		
        for i=1:N
			hold on; plot( inputs(i,:) , '--' , 'Color' , hsv2rgb([i/N 1 1]) , 'LineWidth' , 2 );
			termset{i+1} = sprintf( 'input%d', i );
        end           
		hold on; legend(termset);
	end
	
end