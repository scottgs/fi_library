clear all; close all; clc;

WhichSet = 2; % 1 == UCMERCED (4 nets) and 2 == RSD (3 nets)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load our test data (not seen during training)

if( WhichSet == 1 )
    
    DD1{1} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\caffe_testing_fold_A.csv' );
    DD2{1} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\google_testing_fold_A.csv' );
    DD3{1} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet50_testing_fold_A.csv' );
    DD4{1} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet101_testing_fold_A.csv' );

    DD1{2} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\caffe_testing_fold_B.csv' );
    DD2{2} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\google_testing_fold_B.csv' );
    DD3{2} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet50_testing_fold_B.csv' );
    DD4{2} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet101_testing_fold_B.csv' );

    DD1{3} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\caffe_testing_fold_C.csv' );
    DD2{3} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\google_testing_fold_C.csv' );
    DD3{3} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet50_testing_fold_C.csv' );
    DD4{3} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet101_testing_fold_C.csv' );

    DD1{4} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\caffe_testing_fold_D.csv' );
    DD2{4} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\google_testing_fold_D.csv' );
    DD3{4} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet50_testing_fold_D.csv' );
    DD4{4} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet101_testing_fold_D.csv' );

    DD1{5} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\caffe_testing_fold_E.csv' );
    DD2{5} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\google_testing_fold_E.csv' );
    DD3{5} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet50_testing_fold_E.csv' );
    DD4{5} = Grant_loadcsvfile( 'fusion_files\\UCMERCED\\resnet101_testing_fold_E.csv' );

elseif( WhichSet == 2 )
    
    DD1{1} = Grant_loadcsvfile( 'fusion_files\\RSD\\caffeTest_0_accuracy.csv' );
    DD2{1} = Grant_loadcsvfile( 'fusion_files\\RSD\\google_0_accuracy.csv' );
    DD3{1} = Grant_loadcsvfile( 'fusion_files\\RSD\\ResNet50_0_accuracy.csv' );

    DD1{2} = Grant_loadcsvfile( 'fusion_files\\RSD\\caffeTest_1_accuracy.csv' );
    DD2{2} = Grant_loadcsvfile( 'fusion_files\\RSD\\google_1_accuracy.csv' );
    DD3{2} = Grant_loadcsvfile( 'fusion_files\\RSD\\ResNet50_1_accuracy.csv' );

    DD1{3} = Grant_loadcsvfile( 'fusion_files\\RSD\\caffeTest_2_accuracy.csv' );
    DD2{3} = Grant_loadcsvfile( 'fusion_files\\RSD\\google_2_accuracy.csv' );
    DD3{3} = Grant_loadcsvfile( 'fusion_files\\RSD\\ResNet50_2_accuracy.csv' );

    DD1{4} = Grant_loadcsvfile( 'fusion_files\\RSD\\caffeTest_3_accuracy.csv' );
    DD2{4} = Grant_loadcsvfile( 'fusion_files\\RSD\\google_3_accuracy.csv' );
    DD3{4} = Grant_loadcsvfile( 'fusion_files\\RSD\\ResNet50_3_accuracy.csv' );

    DD1{5} = Grant_loadcsvfile( 'fusion_files\\RSD\\caffeTest_4_accuracy.csv' );
    DD2{5} = Grant_loadcsvfile( 'fusion_files\\RSD\\google_4_accuracy.csv' );
    DD3{5} = Grant_loadcsvfile( 'fusion_files\\RSD\\ResNet50_4_accuracy.csv' );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Our high-level program parameters

% number of inputs (number of DCNNs)
if(WhichSet==1)
    N = 4;
elseif(WhichSet==2)
    N = 3;
end

% how many variables (in the fuzzy measure (FM))
V = 2^N - 1;

% how many ChIs? 
if(WhichSet==1)
    R = 21;
elseif(WhichSet==2)
    R = 19;
end

% lets store our different ChIs
Chis = zeros( R + 1 + 4 + 1 + R ,V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process our folds
ACCSUMS = []; % keep track of performance of solutions
ACCNAMES = []; % keep track of their string IDs 
OUTPUTRESULTS = []; % keep track of what we classified everything as; first field is 1==yes 0==mistake
for fold=1:5

    tdata = DD1{fold};
    tdata2 = DD2{fold};
    tdata3 = DD3{fold};
    if( WhichSet==1 )
        tdata4 = DD4{fold};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Train R different fusions

    % learn the FMs/ChIs
    StoredData = [];
    for i=1:R

        fprintf(1,sprintf('############### LEARN THE CHI: iteration %d ################\n',i));

        % extract test data
        d1 = tdata(:,3+i); % for neuron i, extract output for all samples
        d2 = tdata2(:,3+i); % same ...
        d3 = tdata3(:,3+i);
        if(WhichSet==1)
            d4 = tdata4(:,3+i);
        end
        
        % whats our labels?
        L = tdata(:,1); % get class label out (Grant has them 0 indexed)
        F = find( L == (i-1) ); % find which inds are this class
        L = zeros(size(L)); L(F) = 1; % this is our label vector (1==this class, 0==else)

        % make the data to pass to optimization code (output1 output2 output3 output4 label)
        if(WhichSet==1)
            D = [ d1 d2 d3 d4 L ];
        elseif(WhichSet==2)
            D = [ d1 d2 d3 L ];
        end
        
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
    if(WhichSet==1)
        Chis(R+2,:) = fi_owa( [ 1 0 0 0 ] )'; % max
        Chis(R+3,:) = fi_owa( ones(1,4) ./ 4 )'; % avg
        Chis(R+4,:) = fi_owa( [ 0.05 0.45 0.45 0.05 ] )'; % trimed mean
        Chis(R+5,:) = fi_owa( [ 0 0 0 1 ] )'; % min
    elseif(WhichSet==2)
        Chis(R+2,:) = fi_owa( [ 1 0 0 ] )'; % max
        Chis(R+3,:) = fi_owa( ones(1,3) ./ 3 )'; % avg
        Chis(R+4,:) = fi_owa( [ 0.05 0.9 0.05 ] )'; % trimed mean
        Chis(R+5,:) = fi_owa( [ 0 0 1 ] )'; % min        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now, do the SLFM

    % learn one SLFM

    IndConfMatrix = zeros(R,R,N);
    for i=1:size(tdata,1)

        class = tdata(i,1);

        [val1 ind1] = max( tdata(i,4:end) ); 
        [val2 ind2] = max( tdata2(i,4:end) ); 
        [val3 ind3] = max( tdata3(i,4:end) ); 
        if(WhichSet==1)
            [val4 ind4] = max( tdata4(i,4:end) );
        end
        IndConfMatrix( ind1 , class+1 , 1 ) = IndConfMatrix( ind1 , class+1 , 1 ) + 1; %val1;
        IndConfMatrix( ind2 , class+1 , 2 ) = IndConfMatrix( ind2 , class+1 , 2 ) + 1; %val2;
        IndConfMatrix( ind3 , class+1 , 3 ) = IndConfMatrix( ind3 , class+1 , 3 ) + 1; %val3;
        if(WhichSet==1)
            IndConfMatrix( ind4 , class+1 , 4 ) = IndConfMatrix( ind4 , class+1 , 4 ) + 1; %val4;
        end
        
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
            if(WhichSet==1)
                [val ind] = max( tdata4(k,4:end) ); 
                if( class == ind ) % right!
                    ClassResults(4,1) = ClassResults(4,1) + 1;
                else
                    ClassResults(4,2) = ClassResults(4,2) + 1;            
                end    
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
    title(sprintf('Fold %d: FM variable values',fold));
    % plot them
    figure; 
    for i=1:R
        hold on; stem(Chis(i,:),'Color',hsv2rgb([i/R 1 1]));
    end
    title(sprintf('Fold %d: FM variable values',fold));

    % whats the importance of the variables?
    Shaps = zeros(R,N);
    for i=1:R
       Shaps(i,:) = fi_shapley( Chis(i,:) );
    end
    % plot it
    stds = std(Shaps);
    figure; stem(stds); title(sprintf('Fold %d: difference (std) in FM Shapley values',fold));
    % plot them
    figure; 
    for i=1:R
        hold on; stem(Shaps(i,:),'Color',hsv2rgb([i/R 1 1]));
    end
    title(sprintf('Fold %d: Shapley values',fold));

    % what about interaction inds?
    Ints = zeros(R,N,N);
    for i=1:R
       Ints(i,:,:) = fi_interaction_index( Chis(i,:) );
    end
    % show them
    figure; 
    for i=1:R
        if(WhichSet==1)
            subplot(7,3,i); 
        elseif(WhichSet==2)
            subplot(7,3,i);
        end
        imagesc(squeeze(Ints(i,:,:)),[-1 1]); colormap('gray');
        title(sprintf('Fold %d: neuron %d',fold,i));
    end
    % differences
    figure;
    stds = squeeze( std(Ints(:,:,:)) );
    imagesc( stds , [-1 1] ); title(sprintf('Fold %d: [-1,1] difference (std) in interaction inds',fold));
    figure; 
    imagesc( stds ); title(sprintf('Fold %d: [min,max] difference (std) in interaction inds',fold));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% What aggregation operator?

    AggDistMat = zeros(R,4);
    for i=1:R

        L = Chis(i,:); 

        % compute the max like
        if(WhichSet==1)
            Ws = [ 0.25 0.5 0.75 1 ] ./ (0.25 + 0.5 + 0.75 + 1);
        elseif(WhichSet==2)
            Ws = [ 0.3 0.6 1 ] ./ (0.3 + 0.6 + 1);
        end
        D1 = (Ws(1)/2)*(fi_t1( L , 1 , N ) + fi_t4( L , 1 , N ));
        for k=2:N
            D1 = D1 + ( (Ws(k)/3) * ( fi_t1( L , k , N ) + fi_t2( L , k , N ) + fi_t4( L , k , N ) ) );
        end

        % compute the min like
        if(WhichSet==1)
            Ws2 = [ 1 0.75 0.5 0.25 ] ./ (0.25 + 0.5 + 0.75 + 1);
        elseif(WhichSet==2)
            Ws2 = [ 0.3 0.6 1 ] ./ (0.3 + 0.6 + 1);
        end
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
    fprintf(1,'\n\n\n');
    fprintf(1,'Fold %d: distance values; max, min, avg, owa\n',fold);
    for i=1:R
        for j=1:4
            fprintf(1,'%f ',AggDistMat(i,j));
        end
        fprintf(1,'\n');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% cool way to plot distance matrix

    figure;
    title(sprintf('Fold %d: Distance matrix',fold));
    AggDistMatN = AggDistMat ./ max(AggDistMat(:)); % convert to [0,1] to get grayscale to threshold on
    imagesc(AggDistMat,[0,max(AggDistMat(:))]);
    colormap('gray');
    textStrings = num2str(AggDistMat(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:size(AggDistMat,2),1:size(AggDistMat,1));   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    %textColors = repmat(AggDistMat(:) > midValue,1,3);  %# Choose white or black for the
                                                 %#   text color of the strings so
                                                 %#   they can be easily seen over
                                                 %#   the background color
    %textColors = repmat(AggDistMat(:) > 0,1,3);                                              
    %textColors = repmat(AggDistMat(:) > midValue,1,3);
    textColors = zeros((size(AggDistMat,1)*size(AggDistMat,2)),3);
    for i=1:size(textColors,1)
        if( AggDistMatN(i) > 0.5 )
            textColors(i,:) = [ 0 0 0 ];
        else
            textColors(i,:) = [ 1 1 1 ];
        end
    end
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

    %set(gca,'XTick',1:5,...                         %# Change the axes tick marks
    %        'XTickLabel',{'A','B','C','D','E'},...  %#   and tick labels
    %        'YTick',1:5,...
    %        'YTickLabel',{'A','B','C','D','E'},...
    %        'TickLength',[0 0]);
    if(WhichSet==1)
        set(gca,'XTick',1:4,...                         %# Change the axes tick marks
            'XTickLabel',{'Max','Min','Mean','OWA'},...  %#   and tick labels
            'YTick',1:21,...
            'YTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21'},...
            'TickLength',[0 0]);
    elseif(WhichSet==2)
        set(gca,'XTick',1:4,...                         %# Change the axes tick marks
            'XTickLabel',{'Max','Min','Mean','OWA'},...  %#   and tick labels
            'YTick',1:19,...
            'YTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19'},...
            'TickLength',[0 0]);    
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
    for i=1:size(tdata,1) % go across each data point

        % do per neuron SLFM
        if(1)
            % each class/neuron
            NeuronFires = zeros(1,R);
            for j=1:R

                % get our inputs
                if(WhichSet==1)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
                elseif(WhichSet==2)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) ];                 
                end
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
                if(WhichSet==1)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
                elseif(WhichSet==2)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) ];                     
                end
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
                if(WhichSet==1)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
                elseif(WhichSet==2)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) ];                     
                end
                % fuse them
                NeuronFires(j) = fi_choquet_integral_h_and_g_form( v , Chis(j,:) );

            end

            % take max for the best one
            [val ind] = max( NeuronFires );
            class = tdata(i,1);
            ConfMatrix( ind , class+1 ) = ConfMatrix( ind , class+1 ) + 1;
            % record the classification result (ind is what we called it, class+1 is what it was)
            OUTPUTRESULTS{fold,i,1} = ( ind == (class+1) );
            OUTPUTRESULTS{fold,i,2} = ind;
            OUTPUTRESULTS{fold,i,3} = class+1;
        end

        % do for our max, min, and seeded ones
        if(1)
            for k=1:NumFixedOps
                NeuronFires = zeros(1,R);
                for j=1:R

                    % get our inputs
                    if(WhichSet==1)
                        v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
                    elseif(WhichSet==2)
                        v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) ]; 
                    end
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
                if(WhichSet==1)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
                elseif(WhichSet==2)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) ]; 
                end
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
                if(WhichSet==1)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) tdata4(i,3+j) ]; 
                elseif(WhichSet==2)
                    v = [ tdata(i,3+j) tdata2(i,3+j) tdata3(i,3+j) ]; 
                end
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
            if(WhichSet==1)
                [val4 ind4] = max( tdata4(i,4:end) );
            end
            IndConfMatrix( ind1 , class+1 , 1 ) = IndConfMatrix( ind1 , class+1 , 1 ) + 1;
            IndConfMatrix( ind2 , class+1 , 2 ) = IndConfMatrix( ind2 , class+1 , 2 ) + 1;
            IndConfMatrix( ind3 , class+1 , 3 ) = IndConfMatrix( ind3 , class+1 , 3 ) + 1;
            if(WhichSet==1)
                IndConfMatrix( ind4 , class+1 , 4 ) = IndConfMatrix( ind4 , class+1 , 4 ) + 1;
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot the confusion matrices

    % show the result
    figure; imagesc(ConfMatrix); title(sprintf('Fold %d: confusion matrix (full learned, Choquet)',fold));
    figure; imagesc(SLFMConfMatrix); title(sprintf('Fold %d: confusion matrix (slfm)',fold));
    figure; imagesc(SingleChIConfMatrix); title(sprintf('Fold %d: confusion matrix (single full learned solution)',fold));
    for k=1:NumFixedOps
        figure; imagesc(FixedChIConfMatrix(:,:,k)); title(sprintf('Fold %d: confusion matrix (fixed operators): %d one',fold,k));
    end
    figure; imagesc(PerNeuronSLFMConfMatrix); title(sprintf('Fold %d: confusion matrix (slfm per neuron)',fold));
    figure; imagesc(SugenoConfMatrix); title(sprintf('Fold %d: confusion matrix (full leanred, Sugeno)',fold));

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
    fprintf(1,'Fold %d',fold);
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
    for i=1:NumFixedOps
        fprintf(1,'Fixed method (max,avg,trimmed,min) %d = %f\n',i,accfixed(i));
    end
    fprintf(1,'############################################################\n');
    fprintf(1,'############################################################\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Store those acc's
    
    ACCSUMS{fold,1} = acc;                  ACCNAMES{fold,1} = 'Full Chi (per neuron, Choquet)';
    ACCSUMS{fold,2} = accsug;               ACCNAMES{fold,2} = 'Full Chi (per neuron, Sugeno)';
    ACCSUMS{fold,3} = accchi;               ACCNAMES{fold,3} = 'Full Chi (one shared)';
    ACCSUMS{fold,4} = accslfm;              ACCNAMES{fold,4} = 'Imputed Chi (one SLFM)';    
    ACCSUMS{fold,5} = accslfmperneuron;     ACCNAMES{fold,5} = 'Imputed Chi (SLFM per neuron)';    
    for i=1:N
        ACCSUMS{fold,5+i} = accs(i);        ACCNAMES{fold,5+i} = 'Individual performer';
    end
    for i=1:NumFixedOps
        ACCSUMS{fold,5+N+i} = accfixed(i);    ACCNAMES{fold,5+N+i} = 'Fixed method (max,avg,trimmed,min)';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Vis one of our lattices?

    if(1)
        i = 5; % which neuron
        d1 = tdata(:,3+i); 
        d2 = tdata2(:,3+i); 
        d3 = tdata3(:,3+i);
        if(WhichSet==1)
            d4 = tdata4(:,3+i);
            D = [ d1 d2 d3 d4 ];
        elseif(WhichSet==2)
            D = [ d1 d2 d3 ];
        end
        [Sv Si] = sort(D');
        SortOrder = Si';
        hfig = figure;
        showFMvals = 1;
        showStats = 0;
        showPaths = 1;
        fi_prettyLattice(SortOrder, Chis(i,:), showFMvals, showStats, showPaths);
        title(sprintf('Fold %d: Lattice for neuron %d',fold,i));
    end
    
    %keyboard;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% what are the average scores?

fprintf(1,'####################################################\n');
fprintf(1,'####################################################\n');
fprintf(1,'####################################################\n');
fprintf(1,'####################################################\n');
fprintf(1,'Average scores\n');
for i=1:size(ACCSUMS,2)
    v = 0;
    for j=1:5
        v = v + ACCSUMS{j,i};
    end
    fprintf(1,'%s: score %f\n',ACCNAMES{1,i},v/5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% what did we miss?

fprintf(1,'####################################################\n');
fprintf(1,'####################################################\n');
fprintf(1,'####################################################\n');
fprintf(1,'####################################################\n');
DataSetInds = [ size(DD1{1},1) size(DD1{2},1) size(DD1{3},1) size(DD1{4},1) size(DD1{5},1) ];
for i=1:5
    
    whomiss = [];
    whatdidwecallit = [];
    whatshoulditbe = [];
    for j=1:DataSetInds(i)
        if( OUTPUTRESULTS{i,j,1} == 0 )
            whomiss = [ whomiss j ];
            whatdidwecallit = [ whatdidwecallit OUTPUTRESULTS{i,j,2} ];
            whatshoulditbe = [ whatshoulditbe OUTPUTRESULTS{i,j,3} ];
        end
    end

    fprintf(1,'---- MISSED the following: fold %d ----\n',i);
    for j=1:length(whomiss)
        fprintf(1,'\t File index (1 ind) %d, called it (1 ind) %d, should be (1 ind) %d\n',whomiss(j),whatdidwecallit(j),whatshoulditbe(j));
    end
    
end