clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load our training data (they call augmented data - they did training only, not training + validation here)

[ data ] = Grant_loadcsvfile( 'fusion_files_batch_1\\caffenet_augmented_fold_A.csv' );
[ data2 ] = Grant_loadcsvfile( 'fusion_files_batch_1\\googlenet_augmented_fold_A.csv' );
[ data3 ] = Grant_loadcsvfile( 'fusion_files_batch_1\\resnet50_augmented_fold_A.csv' );
[ data4 ] = Grant_loadcsvfile( 'fusion_files_batch_1\\resnet101_augmented_fold_A.csv' );
% format is [class conf label 21_neuron_output_values] (so 24 length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load our test data (not seen during training)

[ tdata ] = Grant_loadcsvfile( 'fusion_files_batch_1\\caffenet_testing_fold_A.csv' );
[ tdata2 ] = Grant_loadcsvfile( 'fusion_files_batch_1\\google_testing_fold_A.csv' );
[ tdata3 ] = Grant_loadcsvfile( 'fusion_files_batch_1\\resnet50_testing_fold_A.csv' );
[ tdata4 ] = Grant_loadcsvfile( 'fusion_files_batch_1\\resnet101_testing_fold_A.csv' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Our high-level program parameters

% number of inputs (number of DCNNs)
N = 4;

% how many variables (in the fuzzy measure (FM))
V = 2^N - 1;

% how many ChIs? 
R = 21;

% lets store our different ChIs
Chis = zeros( R + 1 + 4 + 1 + R ,V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train 21 different fusions

% learn the FMs/ChIs
StoredData = [];
for i=1:R
    
    fprintf(1,sprintf('############### LEARN THE CHI: iteration %d ################\n',i));
    
    % extract test data
    d1 = tdata(:,3+i); % for neuron i, extract output for all samples
    d2 = tdata2(:,3+i); % same ...
    d3 = tdata3(:,3+i);
    d4 = tdata4(:,3+i);
    
    % whats our labels?
    L = tdata(:,1); % get class label out (Grant has them 0 indexed)
    F = find( L == (i-1) ); % find which inds are this class
    L = zeros(size(L)); L(F) = 1; % this is our label vector (1==this class, 0==else)

    % make the data to pass to optimization code (output1 output2 output3 output4 label)
    D = [ d1 d2 d3 d4 L ];
    
    % remember these data's
    StoredData{i} = D;
    
    % learn it (second parameter is regularizer value; 0 == no reg)
    Chis(i,:) = fi_learn_measure_qp_reg_matlab( D , 0 )';
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Learn a single FM/ChI

Chis(R+1,:) = fi_learn_measure_qp_reg_matlab2( StoredData , 0 )';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try some fixed ones

% how many fixed operators?
NumFixedOps = 4;

% default seed to test against
Chis(R+2,:) = fi_owa( [ 1 0 0 0 ] )'; % max
Chis(R+3,:) = fi_owa( ones(1,4) ./ 4 )'; % avg
Chis(R+4,:) = fi_owa( [ 0.05 0.45 0.45 0.05 ] )'; % trimed mean
Chis(R+5,:) = fi_owa( [ 0 0 0 1 ] )'; % min

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now, do the SLFM

% learn one SLFM

IndConfMatrix = zeros(R,R,N);
for i=1:size(tdata,1)

    class = tdata(i,1);

    [val1 ind1] = max( tdata(i,4:end) ); 
    [val2 ind2] = max( tdata2(i,4:end) ); 
    [val3 ind3] = max( tdata3(i,4:end) ); 
    [val4 ind4] = max( tdata4(i,4:end) );
    IndConfMatrix( ind1 , class+1 , 1 ) = IndConfMatrix( ind1 , class+1 , 1 ) + 1; %val1;
    IndConfMatrix( ind2 , class+1 , 2 ) = IndConfMatrix( ind2 , class+1 , 2 ) + 1; %val2;
    IndConfMatrix( ind3 , class+1 , 3 ) = IndConfMatrix( ind3 , class+1 , 3 ) + 1; %val3;
    IndConfMatrix( ind4 , class+1 , 4 ) = IndConfMatrix( ind4 , class+1 , 4 ) + 1; %val4;

end

% accuracy number
accs = zeros(1,N);
for i=1:N
    accs(i) = sum(diag(IndConfMatrix(:,:,i))) / size(tdata,1);
end

% now, get the Sugeno lambda FM
[slfm, lam] = fi_sugeno_lambda_measure( accs );    
Chis(R+6,:) = slfm;

% now, do per neuron

for i=1:R % loop over classes
    ClassResults = zeros(N,2); % right and wrong
    for k=1:size(tdata,1) % go through each data point
        class = tdata(k,1)+1; % what is the real underlying class
        [val ind] = max( tdata(k,4:end) ); 
        if( class == ind ) % right!
            ClassResults(1,1) = ClassResults(1,1) + 1;
        else
            ClassResults(1,2) = ClassResults(1,2) + 1;            
        end
        [val ind] = max( tdata2(k,4:end) ); 
        if( class == ind ) % right!
            ClassResults(2,1) = ClassResults(2,1) + 1;
        else
            ClassResults(2,2) = ClassResults(2,2) + 1;            
        end       
        [val ind] = max( tdata3(k,4:end) ); 
        if( class == ind ) % right!
            ClassResults(3,1) = ClassResults(3,1) + 1;
        else
            ClassResults(3,2) = ClassResults(3,2) + 1;            
        end           
        [val ind] = max( tdata4(k,4:end) ); 
        if( class == ind ) % right!
            ClassResults(4,1) = ClassResults(4,1) + 1;
        else
            ClassResults(4,2) = ClassResults(4,2) + 1;            
        end            
    end
    ClassResults = ClassResults(:,1) ./ sum(ClassResults')';
    [slfm, lam] = fi_sugeno_lambda_measure( ClassResults' );  
    Chis(R+6+i,:) = slfm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze what we learned (information theoretic indices - focused on the FM)

% what did they learn that was different?
stds = std(Chis(1:R,:));
figure; stem(stds); title('difference (std) in FM values');
% plot them
figure; 
for i=1:R
    hold on; plot(Chis(i,:),'Color',hsv2rgb([i/R 1 1]));
end
title('FM variable values');
% plot them
figure; 
for i=1:R
    hold on; stem(Chis(i,:),'Color',hsv2rgb([i/R 1 1]));
end
title('FM variable values');

% whats the importance of the variables?
Shaps = zeros(R,N);
for i=1:R
   Shaps(i,:) = fi_shapley( Chis(i,:) );
end
% plot it
stds = std(Shaps);
figure; stem(stds); title('difference (std) in FM Shapley values');
% plot them
figure; 
for i=1:R
    hold on; stem(Shaps(i,:),'Color',hsv2rgb([i/R 1 1]));
end
title('Shapley values');

% what about interaction inds?
Ints = zeros(R,N,N);
for i=1:R
   Ints(i,:,:) = fi_interaction_index( Chis(i,:) );
end
% show them
figure; 
for i=1:R
    subplot(7,3,i); 
    imagesc( squeeze(Ints(i,:,:)) , [-1 1] );
    title(sprintf('Interaction index; class %d',i));
end
% differences
stds = squeeze( std(Ints(:,:,:)) );
imagesc( stds , [-1 1] ); title('[-1,1] difference (std) in interaction inds');
imagesc( stds ); title('[min,max] difference (std) in interaction inds');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% What aggregation operator?

AggDistMat = zeros(R,4);
for i=1:R
    
    L = Chis(i,:); 
    
    % compute the max like
    Ws = [ 0.25 0.5 0.75 1 ] ./ (0.25 + 0.5 + 0.75 + 1);
    D1 = (Ws(1)/2)*(fi_t1( L , 1 , N ) + fi_t4( L , 1 , N ));
    for k=2:N
        D1 = D1 + ( (Ws(k)/3) * ( fi_t1( L , k , N ) + fi_t2( L , k , N ) + fi_t4( L , k , N ) ) );
    end
    
    % compute the min like
    Ws2 = [ 1 0.75 0.5 0.25 ] ./ (0.25 + 0.5 + 0.75 + 1);
    D2 = (Ws2(1)/2)*(fi_t3( L , 1 , N ) + fi_t4( L , 1 , N ));
    for k=2:(N-1)
        D2 = D2 + ( (Ws2(k)/3) * ( fi_t3( L , k , N ) + fi_t2( L , k , N ) + fi_t4( L , k , N ) ) );
    end
    
    % mean like
    D3 = 0;
    for k=1:(N-1)
        vals = fetch_vals_at_layer( k , L , N );
        for j=1:length(vals)
            D3 = D3 + abs( vals(j) - (k/N) );
        end
    end
    D3 = ( 1 / (2^N-2) ) * D3;
    
    % LCOS like
    D4 = 0;
    for k=1:(N-1)
        v = fi_t4( L , k , N );
        D4 = D4 + sqrt( v );
    end
    D4 = ( 1 / (N-1) ) * D4;
    
    AggDistMat(i,1) = D1;
    AggDistMat(i,2) = D2;
    AggDistMat(i,3) = D3;
    AggDistMat(i,4) = D4;    
    
end
fprintf(1,'\n\n distance values; max, min, avg, owa\n');
for i=1:R
    for j=1:4
        fprintf(1,'%f ',AggDistMat(i,j));
    end
    fprintf(1,'\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now, eval on test data

% go through, fuse and make decisions
ConfMatrix = zeros(R,R); % ChI learned per neuron
SugenoConfMatrix = zeros(R,R); % ChI learned per neuron
SingleChIConfMatrix = zeros(R,R); % single ChI
IndConfMatrix = zeros(R,R,4); % ind algorithms
SLFMConfMatrix = zeros(R,R); % SLFM 
PerNeuronSLFMConfMatrix = zeros(R,R); % SLFM per neuron
FixedChIConfMatrix = zeros(R,R,NumFixedOps); % our fixed ones (e.g., max)
for i=1:size(tdata,1)

    % do per neuron SLFM
    if(1)
        % each class/neuron
        NeuronFires = zeros(1,R);
        for j=1:R

            % get our inputs
            v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
            % fuse them
            NeuronFires(j) = fi_choquet_integral_h_and_g_form( v , Chis(R+6+j,:) );

        end

        % take max for the best one
        [val ind] = max( NeuronFires );
        class = tdata(i,1);
        PerNeuronSLFMConfMatrix( ind , class+1 ) = PerNeuronSLFMConfMatrix( ind , class+1 ) + 1;
    end    
    
    % Sugeno - do the learned full FM - one per neuron
    if(1)
        % each class/neuron
        NeuronFires = zeros(1,R);
        for j=1:R

            % get our inputs
            v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
            % fuse them
            NeuronFires(j) = fi_sugeno_integral_h_and_g_form( v , Chis(j,:) );

        end

        % take max for the best one
        [val ind] = max( NeuronFires );
        class = tdata(i,1);
        SugenoConfMatrix( ind , class+1 ) = SugenoConfMatrix( ind , class+1 ) + 1;
    end    
    
    % Choquet - do the learned full FM - one per neuron
    if(1)
        % each class/neuron
        NeuronFires = zeros(1,R);
        for j=1:R

            % get our inputs
            v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
            % fuse them
            NeuronFires(j) = fi_choquet_integral_h_and_g_form( v , Chis(j,:) );

        end

        % take max for the best one
        [val ind] = max( NeuronFires );
        class = tdata(i,1);
        ConfMatrix( ind , class+1 ) = ConfMatrix( ind , class+1 ) + 1;
    end

    % do for our max, min, and seeded ones
    if(1)
        for k=1:NumFixedOps
            NeuronFires = zeros(1,R);
            for j=1:R

                % get our inputs
                v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
                % fuse them
                NeuronFires(j) = fi_choquet_integral_h_and_g_form( v , Chis(R+1+k,:) );

            end    

            % take max for the best one
            [val ind] = max( NeuronFires );
            class = tdata(i,1);
            FixedChIConfMatrix( ind , class+1 , k ) = FixedChIConfMatrix( ind , class+1 , k ) + 1; 
        end
    end      
    
    % also do the single learned ChI
    if(1)
        NeuronFires = zeros(1,R);
        for j=1:R

            % get our inputs
            v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
            % fuse them
            NeuronFires(j) = fi_choquet_integral_h_and_g_form( v , Chis(R+1,:) );

        end    

        % take max for the best one
        [val ind] = max( NeuronFires );
        class = tdata(i,1);
        SingleChIConfMatrix( ind , class+1 ) = SingleChIConfMatrix( ind , class+1 ) + 1; 
    end   
    
    % also do the sugeno lambda FM
    if(1)
        NeuronFires = zeros(1,R);
        for j=1:R

            % get our inputs
            v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
            % fuse them
            NeuronFires(j) = fi_choquet_integral_h_and_g_form( v , slfm );

        end    

        % take max for the best one
        [val ind] = max( NeuronFires );
        class = tdata(i,1);
        SLFMConfMatrix( ind , class+1 ) = SLFMConfMatrix( ind , class+1 ) + 1; 
    end
    
    % what is performance of individuals?
    if(1)
        [val1 ind1] = max( tdata(i,4:end) ); 
        [val2 ind2] = max( tdata2(i,4:end) ); 
        [val3 ind3] = max( tdata3(i,4:end) ); 
        [val4 ind4] = max( tdata4(i,4:end) );
        IndConfMatrix( ind1 , class+1 , 1 ) = IndConfMatrix( ind1 , class+1 , 1 ) + 1;
        IndConfMatrix( ind2 , class+1 , 2 ) = IndConfMatrix( ind2 , class+1 , 2 ) + 1;
        IndConfMatrix( ind3 , class+1 , 3 ) = IndConfMatrix( ind3 , class+1 , 3 ) + 1;
        IndConfMatrix( ind4 , class+1 , 4 ) = IndConfMatrix( ind4 , class+1 , 4 ) + 1;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the confusion matrices

% show the result
figure; imagesc(ConfMatrix); title('confusion matrix (full learned, Choquet)');
figure; imagesc(SLFMConfMatrix); title('confusion matrix (slfm)');
figure; imagesc(SingleChIConfMatrix); title('confusion matrix (single full learned solution)');
for k=1:NumFixedOps
    figure; imagesc(FixedChIConfMatrix(:,:,k)); title(sprintf('confusion matrix (fixed operators): %d one',k));
end
figure; imagesc(PerNeuronSLFMConfMatrix); title('confusion matrix (slfm per neuron)');
figure; imagesc(SugenoConfMatrix); title('confusion matrix (full leanred, Sugeno)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the acc numbers

% accuracy number
acc = sum(diag(ConfMatrix)) / sum(ConfMatrix(:));
accsug = sum(diag(SugenoConfMatrix)) / sum(SugenoConfMatrix(:));
accslfm = sum(diag(SLFMConfMatrix)) / sum(SLFMConfMatrix(:));
accchi = sum(diag(SingleChIConfMatrix)) / sum(SingleChIConfMatrix(:)); 
accs = zeros(1,N);
for i=1:N
    accs(i) = sum(diag(IndConfMatrix(:,:,i))) / sum(sum(IndConfMatrix(:,:,i)));
end
accfixed = zeros(1,NumFixedOps);
for k=1:NumFixedOps
    accfixed(k) = sum(diag(FixedChIConfMatrix(:,:,k))) / sum(sum(FixedChIConfMatrix(:,:,k)));
end
accslfmperneuron = sum(diag(PerNeuronSLFMConfMatrix)) / sum(PerNeuronSLFMConfMatrix(:));
% print it out
figure; stem( [acc accs accchi accfixed accslfmperneuron accsug] ); title('performance'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print out the acc numbers

fprintf(1,'############################################################\n');
fprintf(1,'############################################################\n');
fprintf(1,'Full Chi (per neuron, Choquet) = %f\n',acc);
fprintf(1,'Full Chi (per neuron, Sugeno) = %f\n',accsug);
fprintf(1,'Full Chi (one shared) = %f\n',accchi);
fprintf(1,'Imputed Chi (one SLFM) = %f\n',accslfm);
fprintf(1,'Imputed Chi (SLFM per neuron) = %f\n',accslfmperneuron);
for i=1:N
    fprintf(1,'Individual performer %d = %f\n',i,accs(i));
end
for i=1:N
    fprintf(1,'Fixed method (max,avg,trimmed,min) %d = %f\n',i,accfixed(i));
end
fprintf(1,'############################################################\n');
fprintf(1,'############################################################\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vis one of our lattices?

i = 1; % which neuron
d1 = tdata(:,3+i); 
d2 = tdata2(:,3+i); 
d3 = tdata3(:,3+i);
d4 = tdata4(:,3+i);
D = [ d1 d2 d3 d4 ];
[Sv Si] = sort(D');
SortOrder = Si';
hfig = figure(1);
showFMvals = 0;
showStats = 0;
showPaths = 1;
fi_prettyLattice(SortOrder, Chis(R+1,:), showFMvals, showStats, showPaths);
