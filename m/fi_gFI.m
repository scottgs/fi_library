function [ FIresult ] = fi_gFI( inputs , FMbinary , whichintegral , noalphacuts , plotit )
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
%	[ FIresult ] = fi_gFI( inputs , FMbinary , whichintegral , plotit )
%	Inputs:
%		inputs = cell structure
%			inputs{1}.type = 'gauss'
%			inputs{1}.params = [ 0.3 0.1 1 ]  for sigma, mean and height
%			inputs{2}.type = 'trap'
%			inputs{2}.params = [ 0.2 0.4 0.6 0.7 0.5 ] for a,b,c,d and 0.5 is height
%			inputs{3}.type = 'tri'
%			inputs{3}.params = [ 0.5 0.7 0.9 0.4 ] for a,b,c and 0.4 is height
%		FMbinary = binary, e.g., [ x1 x2 x12 x3 x13 x23 x123 ]
%		whichintegral = 0 is Sugeno integral and 1 is Choquet integral based
%    	noalphacuts = how many alpha cuts on [0,1] do you want? 
%		plotit = 1 for if you want to plot the result
%	Outputs:
%		FIresult = result of gFI 
%			e.g., FIresult.alpha - vector of alpha cut values
%			e.g., FIresult.intervals - results at those alpha cuts
%
% Example
%   inputs{1}.type = 'gauss';
%   inputs{1}.params = [ 0.1 0.5 1 ];
%   inputs{2}.type = 'trap';
%   inputs{2}.params = [ 0.1 0.2 0.3 0.4 0.75 ];
%   inputs{3}.type = 'tri';
%   inputs{3}.params = [ 0.7 0.9 1 ];
%   inputs{4}.type = 'gauss';
%   inputs{4}.params = [ 0.2 0.9 1 ];
%   FMbinary = fi_owa( ones(1,4) ./ 4 )';
%   whichintegral = 1;
%   noalphacuts = 100;
%   plotit = 1;
%   [ FIresult ] = fi_gFI( inputs , FMbinary , whichintegral , noalphacuts , plotit );
%
% Program description:
%   Compute the discrete (finite X) gFI for set-valued data
%	Assumes the FM is in "binary coded" format
% 	Note, we did implment subnormal
% 	Note, we did not implement non-convex inputs, easy to do though
%		See D. T. Anderson, T. C. Havens, C. Wagner, J. M. Keller, M. F. Anderson, D. J. Wescott, "Extension of the Fuzzy Integral for General Fuzzy Set-Valued Information," IEEE Transactions on Fuzzy Systems, vol. 22 (6), pp. 1625-1639, 2014
	
	% how many inputs
	N = length(inputs);
	
	% results
	FIresult.alpha = eps:(1/noalphacuts):1;
	FIresult.intervals = zeros(length(FIresult.alpha),2);
	
	% max height
	MaxHeight = 1;
	for i=1:N
		MaxHeight = min( MaxHeight , inputs{i}.params(end) );
	end
	
	% do the gFI
	i=1;
	if(plotit==1)
		figure;
    end
	for a=(FIresult.alpha)
        % get the intervals
        Ints = zeros(N,2);
        for j=1:N
            if( strcmp('gauss',inputs{j}.type) == 1 ) 
                if( a <= inputs{j}.params(end) )
                    n1=1;
                    n2=(-1)*(2*inputs{j}.params(2));
                    n3=inputs{j}.params(2)*inputs{j}.params(2)+log((a*(1/inputs{j}.params(end))))*2*inputs{j}.params(1)*inputs{j}.params(1);
                    sols1=( (-1)*n2 + sqrt(n2*n2 - 4*n1*n3) ) / ( 2*n1 );
                    sols2=( (-1)*n2 - sqrt(n2*n2 - 4*n1*n3) ) / ( 2*n1 );
                    Ints(j,1)=min([sols1 sols2]);
                    Ints(j,2)=max([sols1 sols2]);	
                end
            elseif( strcmp('trap',inputs{j}.type) == 1 ) 
                if( a <= inputs{j}.params(end) )  
                    if(a==inputs{j}.params(end))
                        Ints(j,1) = inputs{j}.params(2);
                        Ints(j,2) = inputs{j}.params(3);
                    else
                        Ints(j,1)=(a*(1/inputs{j}.params(end)))*(inputs{j}.params(2)-inputs{j}.params(1))+inputs{j}.params(1);
                        Ints(j,2)=inputs{j}.params(4)-(a*(1/inputs{j}.params(end)))*(inputs{j}.params(4)-inputs{j}.params(3));
                    end
                end
            elseif( strcmp('tri',inputs{j}.type) == 1 ) 
                if( a <= inputs{j}.params(end) )  
                    if(a==inputs{j}.params(end))
                        Ints(j,1) = inputs{j}.params(2);
                        Ints(j,2) = inputs{j}.params(2);                        
                    else
                        Ints(j,1)=(a*(1/inputs{j}.params(end)))*(inputs{j}.params(2)-inputs{j}.params(1))+inputs{j}.params(1);
                        Ints(j,2)=inputs{j}.params(3)-(a*(1/inputs{j}.params(end)))*(inputs{j}.params(3)-inputs{j}.params(2));
                    end
                end
            end
        end
        % plot them?
        if(plotit==1)
            hold on;
            for j=1:N
                if( a <= inputs{j}.params(end) )                
                    plot( [Ints(j,1) Ints(j,2)] , [a a]+(rand()*0.2*(1/noalphacuts)) , '-' , 'Color', hsv2rgb([j/N 1 1]) );
                end
            end
            hold on; 
            axis( [0 1 0 1] );
        end	         
        % compute the FI in interval fashion
		if( a <= MaxHeight )			
			if(whichintegral==0)
				FIresult.intervals(i,1) = fi_sugeno_integral_h_and_g_form( Ints(:,1)' , FMbinary );
				FIresult.intervals(i,2) = fi_sugeno_integral_h_and_g_form( Ints(:,2)' , FMbinary );	
            else                
				FIresult.intervals(i,1) = fi_choquet_integral_h_and_g_form( Ints(:,1)' , FMbinary );
				FIresult.intervals(i,2) = fi_choquet_integral_h_and_g_form( Ints(:,2)' , FMbinary );			
            end		
            % plot?
            if(plotit==1)
                hold on;
                plot( [FIresult.intervals(i,1) FIresult.intervals(i,2)] , [a a] , '--k' );
                hold on; 
                axis( [0 1 0 1] );
            end	               
        end     
		i=i+1;
    end
     
    % label the inputs
    if(plotit==1)    
        legendnames = [];
        for j=1:N
            legendnames{j} = sprintf('Input %d: %s',j,inputs{j}.type);
        end
        legendnames{N+1} = 'Fuzzy integral';
        hold on; 
        legend(legendnames);
    end
    
end