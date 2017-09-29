function [ interval ] = fi_choquet_integral_interval( inputs , FM , plotit )
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
%	[interval] = fi_choquet_integral_interval( inputs , FM )
%	Inputs:
%		inputs = e.g., [ 0.2 0.3 ; 0.1 0.8 ] for two intervals of [0.2,0.3] and [0.1,0.8]
%		FM = binary, e.g., [ x1 x2 x12 x3 x13 x23 x123 ]
%       plotit = plot the result?
%	Outputs:
%		interval = result of ChI (an interval)
%
% Program description:
%   Compute the discrete (finite X) Choquet integral for interval-valued data
%	Assumes the FM is in "binary coded" format

	% result
	interval = zeros(1,2);
	
	% do on left most endpoints
	interval(1) = fi_choquet_integral_h_and_g_form( inputs(:,1)' , FM );
	% do on right most endpoints
	interval(2) = fi_choquet_integral_h_and_g_form( inputs(:,2)' , FM );
    
    % plot it?
    if( plotit == 1 )
        figure;
        legendnames = [];
        for i=1:size(inputs,1)
            hold on;
            plot( [inputs(i,1) inputs(i,2)] , [i i] , 'Color', hsv2rgb([i/size(inputs,1) 1 1]) );
            legendnames{i} = sprintf('Input %d',i);
        end
        hold on;
        plot( [interval(1) interval(2)] , [i+1 i+1] , 'k' );
        legendnames{i+1} = sprintf('ChI result');
        axis( [ 0 1 0 size(inputs,1)+2 ] );
        legend(legendnames);
    end

end