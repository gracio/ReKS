function [discComp,U,S] = clusterKmeans(A,A2,dm,U,S)

%function [discComp,U,S] = clusterKmeans(Pts,A,dm,U,S)

% grace comment: A(required) is the affinity matrix(made sparse), dm(optional) is number of
% clusters considered, U(optional) is eigenvectors in columns, S(optinoal) is eigenvalues.

FALSE = (0 == 1);
TRUE = ~FALSE;

Pts = A;

% grace modification: do SVD when both A,dm are supplied, or only A is
% supplied
if (nargin <= 3) % or 3
    
    %A = A/median((sum(A,1)));
    if ~isempty(A2)
        D = sum(A2,1)';
        fprintf('use A2 \n');
    else
        D = sum(A, 1)';
    end
    % Normalize column sum to one.
    sqrtD = D .^ 0.5;
    Q = (sqrtD .^ -1) * ones(1, length(D));
    Mcut = Q .* A .* Q';         % M = D^-0.5 Markov D^0.5
    %fprintf(['ishermitianMcut: ' num2str(ishermitian(Mcut)) '\n'])
    Mcut(isnan(Mcut)) = 0;
    % grace modification: change from svd to svds, examine top 20
    % eigenvectors, or whatever the group size is if less than 20
    if (nargin < 3) % if only A is supplied
        [U S V] = svds(Mcut,min(20,size(A,1)));
    else % if dm is also supplied
        [U S V] = svds(Mcut,dm);
    end
    S = diag(S);
    %[U,S] = eigs(Mcut,20,'LM');
    clear Q
    clear V
    
    % grace modification: handle individual components
    % if there are individual components i.e. more than one S==1
    if length(find(abs(S'-1)< 0.000001))>1
        fprintf('connected components detected.\n')
        S
        numConnected = length(find(abs(S'-1)< 0.000001));
        [grouping nclasses] = gingcca(A,0);
        uqdisc = unique(grouping);
        for i=1:length(uqdisc)
            discComp(grouping==uqdisc(i),i) = 1;
        end
        
        return
    end
    
end

if (nargin < 3) | (nargin > 3)
    % grace modification: determine dm using eigengap only if dm is not
    % supplied
    [maxy,maxi] = max((S(1:length(S)-1)-S(2:length(S))) / S(1));
    dm = maxi + 1;
end

% grace modification: moved from begining.
% variable dm gets overwritten at the end
dmo = dm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Projection Clustering Using K-means %%%%%%%%%%%%%%%%%%%%%
%%%%  See Ng, Jordan, and Wiess           %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Params for selecting initial K-means centers
minRes2 = 0.25;
MARKOV = FALSE; % MARKOV = TRUE; minRes2 = 0.5; nIts = 64;


%%Rerun: rand('state', saveRandState);

n = size(Pts);
nIts  = 2;
%% Set up data for center selection:
if MARKOV
    c = U .* (ones(n(1), 1) * (S .^ nIts)');
    E = c;
    E = ((sum(E.*E, 2).^-0.5) * ones(1, size(E,2))) .* E;
else
    E = U(:, 1:dm);
    E = ((sum(E.*E, 2).^-0.5) * ones(1, dm)) .* E;
end

%% D0 center selection:
%% Rerun: rand('state', saveRandState);
saveRandState = rand('state');

%% Select initial centers to be nearly orthogonal.
j = round(0.5 + n(1) * rand(1));
j = max(j, 1); j = min(j,n(1));
c = E(j,:);
res = E;
k = 2;
while MARKOV | k <= dm
    res = res - (res * (c(k-1,:)')) * c(k-1,:);
    nrmRes2 = sum(res .* res, 2);
    samp = cumsum(nrmRes2 .* (nrmRes2>minRes2));
    if samp(n(1)) == 0 | isnan(samp(n(1))) | isinf(samp(n(1)))
        break;
    end
    samp = samp/samp(n(1));
    r = rand(1);
    idx = find(samp>=r);
    if any(idx)
        c = [c ; E(idx(1), :)];
        k = k+1;
    else
        error('Random draw fanned!??!');
    end
end
k = k-1;
if k < dm & ~MARKOV
    fprintf(2,' Got only %d basis elements\n', k);
    dm = k;
else
    if (k ~= dmo)
        fprintf(2,' Got %d basis elements\n', k);
    end
    dm = k;
end

%% Call kmeans
options = foptions;
[centers options post errlog] = my_kmeans(c, E, options);
discComp = post>0;
