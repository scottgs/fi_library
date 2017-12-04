function [FMbinary] = fi_learn_measure_qp_reg_matlab2( data_sets_with_labels , reg )

    if (nargin == 1)
        reg = 0;
    end

    NumOfOutputs = length(data_sets_with_labels);
    data_set = data_sets_with_labels{1}(:,1:end-1);
    labels = data_sets_with_labels{1}(:,end);
    NumOfInputs = size(data_set,2);
    g = 2^(NumOfInputs)-1;
    rollingD = zeros(g,g);
    rollingf = zeros(g,1);
    for i=1:NumOfOutputs

        data_set = data_sets_with_labels{i}(:,1:end-1);
        labels = data_sets_with_labels{i}(:,end);

        NumOfInputs = size(data_set,2);

        % Ceofficients and Constraints for the QP
        [C,D,f,variablesInO]=QPmatrices(data_set,labels);
        % options = optimset('Display','off');

        rollingD = rollingD + D;
        rollingf = rollingf + f;

    end

    if (reg<=0)
      FMbinary = quadprog(2*rollingD,rollingf,C,zeros(size(C,1),1),[],[],[zeros(g-1,1); 1],ones(g,1));
    else
      regConstr = ones(1,g);
      rhs = 1/reg;
      FMbinary = quadprog(2*rollingD,rollingf,[C; regConstr] ,[zeros(size(C,1),1);rhs],[],[],[zeros(g-1,1); 1],ones(g,1));
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