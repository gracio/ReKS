function [final_newA final_allii final_alljj final_numE final_tau final_numU final_levelLabel, final_ComponentIndices,S] = myPerturbSnip(A,beta0,tauVal, epsilon, plotFigs,valOrPrctile,levelLabel,componentIndices)
% beta0values = 10:10:100
% tauPrctilevalues = 0.001:0.001:0.01
% epsilonvalues = [1/8 1/4 1/6 1/3 1/2]
% tauVal = -0.5:0.05:-0.05
%
% plotFigs: 1 if one wishes to plot figures, 0 if not
% valOrPrctile: 'prctile' to use prctile method to pick tau, 'val' to use a
% fraction of the median to pick tau


k=1
allii=[];
alljj=[];
tau = [];
numE = 0;
numU = 0;
ii=1;
S=[];
tau = NaN;

% keep performing iterations until there are no eligible eigenvectors to
% use, or disconnected components have been found. only consider A's with
% size > 2
while ~isempty(ii) && length(A)>2
    numel(A)
    ii=[];
    % step 1-------------------------------------------------------------
    % A is given as input
    
    % step 2 ------------------------------------------------------------
    D = sum(A,2);
    delta = median(D);
    L = diag(D.^-0.5)*A*diag(D.^-0.5) ;
    
    % step 3 ------------------------------------------------------------
    [U,S,V] = svds(L,min(20,size(L,2)));
    S = diag(S);
    
    
    % step 4 ------------------------------------------------------------
    % calculate beta_k, the half life of the eigenvector u, <eq.7.19>
    for i=1:length(S)
        beta_k(i) = -1*(log(2)/log(S(i)));
    end
    
    % to be used in perturbBeta
    DD = D.^0.5;
    if (size(DD,2) ~= 1)
        DD = diag(DD);
    end
    
    % the part of matrix that originally had affinity value
    binNhbr = full(A) > 0;
    dLogHalfAll = [];
    dLogHalfMinStack = [];
    k2use = setdiff(find(beta_k>beta0*epsilon),find(abs(S'-1)< 0.000001)); % exclude eigenvectors with singular values of 1. take care of floating point precision
    fprintf(['eigenvectors with eigenvalues above the thresholds and <1: ' num2str(k2use) '\n'])
    
    % if one goes ahaed with perturbation
    if ~isempty(k2use) && length(find(abs(S'-1)< 0.000001))==1
        
        
        % take all Uk with half life > beta0* epsilon
        for i = 1:length(k2use)
            k2use(i)
            id = k2use(i);
            % perturbBeta implements <equation 7.28>
            fprintf('performing eigenCut\n')
            dLogHalf{id} = perturbBeta(beta0,S(id),U(:,id),DD) ;
            dLogHalf{id} = dLogHalf{id} .* binNhbr; % only take parts of the matrix that originally had affinity value
            
            % dLogHalfMinStack is keeping track of all the lowest sensitivity values
            % across eigenvectors
            if isempty(dLogHalfMinStack)
                dLogHalfMinStack = dLogHalf{id};
            else
                dLogHalfMinStack = min(dLogHalfMinStack,dLogHalf{id});
            end
            fprintf('finished. now off to do non-maximal supressions\n')
            
            %step 5 non-maximal supression...potential speedup available-----
            minvali = full(min(dLogHalf{id},[],1));
            minvalj = full(min(dLogHalf{id},[],2));
            fprintf('now\n')
            isMinimal{id} = (full(dLogHalf{id}) <= repmat(minvali,size(dLogHalf{id},1),1)) & (dLogHalf{id} <= repmat(minvalj,1,size(dLogHalf{id},2)));
            fprintf('done\n')
            
            dLogHalfAll = [dLogHalfAll ; dLogHalf{id}(find(binNhbr))]; % put all the half life sensitivity values in one single vector for plotting later
        end
        
        % step 6---------------------------------------------------------
        if strmatch(valOrPrctile,'val')
            % define tau relative to delta
            tau(k) = tauVal/delta;
            fprintf(['tau is fractionValue(' num2str(tauVal) ')/delta(' num2str(delta) ').\n'])
        else
            % define tau relative to percentile
            [f,xx]=ecdf(dLogHalfAll);
            % using tau percentile to find a tau value as the cut-off
            tp = find(f<tauVal);
            tau(k) = xx(tp(end));
            fprintf(['tau is prctile(' num2str(tauVal) ')\n'])
        end
        
        % find all edges less than the tau cut-off to snip in the upper
        % triangle for all eigenvectors, make sure they are minimal
        toCut = zeros(size(A));
        for i = 1:length(k2use)
            id = k2use(i);
            toCut_id = (triu(dLogHalf{id},1)<tau(k)) & (triu(dLogHalf{id},1)~=0) & (isMinimal{id});
            %toCut_id = (triu(dLogHalf{id},1)<tau(k)) & (triu(dLogHalf{id},1)~=0) ;
            toCut = toCut | toCut_id;
        end
        [ii,jj]= find(toCut);
        
        % set edges snipped to 0 (snip) for both upper and lower triangle
        % indices
        newA = A;
        for i=1:length(ii)
            newA(ii(i),jj(i))=0;
            newA(jj(i),ii(i))=0;
            % transfer weights of the edges that are snip to diagonal of
            % the affinity matrix
            newA(jj(i),jj(i)) = newA(jj(i),jj(i)) + A(ii(i),jj(i));
            newA(jj(i),jj(i)) = newA(ii(i),ii(i)) + A(ii(i),jj(i));
        end
        
        % plotting functionalities -------------------------------------
        if plotFigs
            
            % parse the levelLabel here
            C=textscan(levelLabel,'%f','Delimiter','-');
            C=C{1};
            figNum = 0;
            for i=1:length(C)
                figNum = figNum*100+C(i);
            end
            figNum = figNum*100;
            
            %fig = gcf;
            fig = figure(figNum+k)
            set(fig, 'Position', [100+k*50 1000 1600 1000]);
            
            % make the half life sensitivities at least dLogHalfMin for display
            dLogHalfMin    = -2;
            dLogHalfAll = max(dLogHalfAll, dLogHalfMin);
            
            % plot cdf of sensitivity
            subplot(2,3,1)
            cdfplot(dLogHalfAll)
            title('CDF plot of half life sensitivity')
            xlabel('sensitivity')
            
            
            % plot affinity heatmap
            subplot(2,3,2)
            imagesc(A)
            colorbar
            title([levelLabel 'Affinity heatmap: before' ])
            fprintf(['non zero elements: ' num2str(nnz(A)) '\n'])
            
            
            subplot(2,3,3)
            imagesc(newA)
            colorbar
            title('Affinity heatmap: after')
            fprintf(['non zero elements: ' num2str(nnz(newA)) '\n'])
            
            
            % plot sensitivity frequency and tau
            subplot(2,3,4)
            totBins = 101; % 101 bins
            [ht,bn] = hist(dLogHalfAll,totBins); % obtain frequency at 101 bins
            ht = ht/sum(ht)/abs(bn(1)-bn(2)); % normalize the frequency
            [mm,xx] = max(log10(ht+1));
            grid on;
            plot(bn,log10(ht+1),'k-','linewidth',2);
            hold on; grid on;
            %%set(gca,'fontsize',20);
            xlabel('dlog(\beta + \beta_0)/d\alpha_{i,j}')
            ylabel('Log Frequency');
            title({['iteration ' num2str(k)]; [' removed ' num2str(length(ii)) ' elements']})
            ax = axis;
            %%%set(gca,'xtick',-0.2:0.05:0.2);
            plot([tau(k) tau(k)],[ax(3) ax(4)],'k--','linewidth',2)
            %%%text(tr-0.001,2,'\tau','fontsize',50,'fontweight','bold','color',[0 0 0])
            hold off
            
            %plot minimal sensitivity across eigenvectors
            subplot(2,3,5)
            imagesc(dLogHalfMinStack)
            colorbar
            title('sensitivity heatmap')
            
            
            %plot edges snipped
            subplot(2,3,6)
            toCut = toCut | toCut';
            imagesc(toCut*100)
            title('edges removed')
            
            %close(fig)
        end
        
        % display edges snipped ------------------------------------------
        % display number of edges snipped
        numE(k) = length(ii);
        fprintf(['num edges snipped: ' num2str(numE(k)) '\n'])
        % actual tau used
        fprintf(['tau: ' num2str(tau(k)) '\n'])
        % number of eigen modes used
        numU(k) = length(k2use);
        
        fprintf(['num eigen modes used: ' num2str(numU(k)) '\n'])
        fprintf(['ii: ' num2str(reshape(ii,1,length(ii))) '\n'])
        fprintf(['jj: ' num2str(reshape(jj,1,length(jj))) '\n'])
        
        % advance to next iteration, reassign A
        k=k+1
        allii=[allii;ii];
        alljj=[alljj;jj];
        A = newA;
        
    end %end of if-statement: there exists eigenvectors with high eigenvalues && there is only one eigenvalue==1
    
    
end % end of while loop for k


% prepare to return these guys unless changed further
% this A is the same size as the input - whatever the input A is
final_newA = A;
% relative to component indices, translate back to component
final_allii = componentIndices(allii);
final_alljj = componentIndices(alljj);
final_tau = tau;
final_numE = numE;
final_numU = numU;
final_levelLabel = levelLabel;
final_ComponentIndices = {componentIndices};

% % if disconnected components are found
% if ( length( find( abs(S'-1)< 0.000001 ) ) >1 )
%     % display connected components at this point
%     numConnected = length(find(abs(S'-1)< 0.000001));
%     grouping = kmeans(U(:,1:numConnected),numConnected);
%     for i=1:numConnected
%         components{i} = find(grouping==i); % potentially a problem
%     end
%     fprintf(['Connected Components: ' num2str(numConnected) '\n'])
%     fprintf(['number of elements: ' num2str(length(components)) '\n']) %
%     
%     % recurse down
%     for i=1:numConnected
%         
%         % A and level are altered
%         [sub_newA sub_allii sub_alljj sub_numE sub_tau sub_numU sub_levelLabel sub_ComponentIndices] = myPerturbSnip(A(components{i},components{i}),beta0,tauVal, epsilon, plotFigs,valOrPrctile,[levelLabel '-' num2str(i)], componentIndices(components{i}));
%         
%         % update output
%         % reorder indices
%         final_newA(components{i},components{i}) = sub_newA;
%         final_allii= [final_allii sub_allii];
%         final_alljj= [final_alljj sub_alljj];
%         
%         % append
%         final_tau = [final_tau sub_tau];
%         final_numE = [final_numE sub_numE];
%         final_numU = [final_numU sub_numU];
%         
%         if iscellstr(sub_levelLabel)
%             final_levelLabel = [final_levelLabel sub_levelLabel];
%         else
%         final_levelLabel = [final_levelLabel {sub_levelLabel}];
%         end
%         final_ComponentIndices = [final_ComponentIndices sub_ComponentIndices];
%     end % go through every component
%     
%     
% end





