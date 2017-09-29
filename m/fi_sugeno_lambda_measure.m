function [FMbinary, lambda] = fi_sugeno_lambda_measure( densities )
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
%		densities = the N densities (for N inputs)
%	Outputs:
%		FMbinary = the resultant fuzzy measure
%		lambda = the Sugeno lambda value

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first, calc lambda
    %    note: not doing any error checking, you can
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % how many densities?
    N = length( densities );

    % some conditioning of the densities
    densities = densities + eps; 
    
    % hold onto the polynomial
    Coefficents = [];
    
    % compute the first, which is the sum of all the densities
    Coefficents(N) = sum( densities );

    % now compute the rest
    for i=N:-1:2
        
        %Generate the combinations ... get all pairs, 3-tuples, etc., ....
        Combos = nchoosek(1:N,i);
        %How many elements were generated?
        [NumberOf Elems] = size(Combos);
        %Reset the Coefficent value to start
        Coefficents(N-i+1) = 0;
        %Now, for each coefficent, compute its product, summing them all up
        for j=1:NumberOf
            %Tmp used for product calculation
            tmp = 1;
            %Now, for each combo, compute the product
            for k=1:Elems
                tmp = tmp * densities(Combos(j,k));
            end
            %Sum these up
            Coefficents(N-i+1) = Coefficents(N-i+1) + tmp;
        end
        
    end
    
    %Following equation to solve for has to equal 1, so subtract from scalar
    Coefficents(N) = Coefficents(N) - 1;
    
    %Now use matlab's roots fx (solution is the largest real root)
    lambda = max(real(roots(Coefficents)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now, calc the full measure 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get our binary mapping
    nVars = 2^N-1;
    varsBin = de2bi(1:nVars, N);   
    
    % our result
    FMbinary = zeros(1,2^N-1);
    
    % which ones are the densities?
    for i=1:N
        FMbinary( 2^(i-1) ) = densities(i);
    end
    % now do the rest
    for i=2:1:N
        
        % generate the combinations ... get all pairs, 3-tuples, etc., ....
        Combos = nchoosek(1:N,i);
        % how many elements were generated?
        [NumberOf Elems] = size(Combos);
        for j=1:NumberOf
            
            ThisSet = sort(Combos(j,:)); % sort numerically, e.g., [1 3 4]
            ThisSet1 = ThisSet(1); % get the first element, e.g., [1]
            ThisSet2 = ThisSet(2:end); % get the remaining elements, e.g., [3 4]
            
            % fetch their values calc'd so far
            Val1 = cumsum(2.^[([ThisSet1])-1]);
            Term1 = FMbinary(Val1(end));
            Val2 = cumsum(2.^[([ThisSet2])-1]);
            Term2 = FMbinary(Val2(end));
            
            % final result (per variable)
            Val3 = cumsum(2.^[([ThisSet])-1]);
            % e.g., so doing g1 + g34 + lambda * g1 * g34
            FMbinary(Val3(end)) = Term1 + Term2 + lambda * Term1 * Term2;
            
        end
        
    end    
    
end