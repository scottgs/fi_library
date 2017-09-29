function [FMbinary] = fi_learn_measure_qp_reg_matlab( data_set_with_labels , reg )
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
%		reg = regularizer value
%	Outputs:
%		FMbinary = resultant fuzzy measure in binary coded format
%
% Program description:
%   learn the fuzzy measure and ChI w.r.t. quadratic programming and regularization

if (nargin == 1)
reg = 0;
end
data_set = data_set_with_labels(:,1:end-1);
labels = data_set_with_labels(:,end);

NumOfInputs = size(data_set,2);
g = 2^(NumOfInputs)-1; 

% Ceofficients and Constraints for the QP
[C,D,f, variablesInO]=QPmatrices(data_set,labels);
% options = optimset('Display','off');

if (reg<=0)
  FMbinary = quadprog(2*D,f,C,zeros(size(C,1),1),[],[],[zeros(g-1,1); 1],ones(g,1));
else
  regConstr = ones(1,g);
  rhs = 1/reg;
  FMbinary = quadprog(2*D,f,[C; regConstr] ,[zeros(size(C,1),1);rhs],[],[],[zeros(g-1,1); 1],ones(g,1));
end

% Standard QP minimization
end

function [C,D,f, variablesInT] = QPmatrices(T,alpha)
%QPmatrices_v2 Generate the coefficient and constraint matrices for the QP.
%   [C,D,f] = QPmatrices_v2(H,alpha) generates coefficients D and F; and 
%   inequality monotonicity constraints, C for quadprog from tranining data, T and 
%   training label vector, alpha 
%   
%   T: Matrix of training observation (N X M )
%   alpha: Training labels (N X 1)
%
%   Note: the index of a variables in the returned matrices corresponds to 
%   the decimal valued binary position of the FM variable

[no,N]=size(T);

% Form the constraint matrix
nc = 2^N-1;
k = [0:N-2];
nr = sum(factorial(N)./(factorial(N-k-1).*factorial(k)));
C = zeros(nr,nc);

index = 1;
n = [0:N-1];
for i=1:N,
    A = nchoosek(n,i);
    ncom = size(A,1);
    B = cell2mat(cellfun(@(x)setdiff(n,x),num2cell(A,2),'UniformOutput',false));
    cA = sum(2.^A,2);
    for j=1:N-i,
        %[A B(:,j)]
        rc = sub2ind([nr nc],[index:index+ncom-1]',cA); C(rc)=1;
        rc = sub2ind([nr nc],[index:index+ncom-1]',sum(2.^[A B(:,j)],2)); C(rc)=-1;
        index = index + ncom;
    end;
end;

[SortVal, SortInd] = sort( T , 2, 'descend' );

%Append a 0 for the difference calculation below
SortVal = [SortVal zeros(no,1)];

% These are the A vectors (without all the zeros)
Hdiff = SortVal(:,1:end-1)-SortVal(:,2:end);

% Compute the index of the values to fill the D and f matrices
i = cumsum(2.^(SortInd-1),2);

D = zeros(nc,nc);
f = zeros(nc,1);
for j = 1:no,
    D(i(j,:),i(j,:)) = D(i(j,:),i(j,:)) + Hdiff(j,:)'*Hdiff(j,:);
    f(i(j,:)) = f(i(j,:)) -2*alpha(j)*Hdiff(j,:)';
end;
variablesInT = unique(sort(i));
end