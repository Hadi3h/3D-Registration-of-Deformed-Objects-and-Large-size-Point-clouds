% ICP algorithm: adapted from  Per Bergström's implementation 
% ICP algorithm takes as input TR (rotation matrix), TT (translation matrix), Model
% and Data points. 
%The output of the algorithm is the transformation matrices (TR and TT) and
%res is the error of registration. Data is a matrix of the aligned data points
function [TR,TT, res, Data] = icp(TR, TT,Data, Model)
global ICP_iter;
global ICP_thres;


if nargin<2
    
    error('To few input arguments');
end

if or(isempty(Model),isempty(Data))
    error('Something is wrong with the Model points and Data points');
end


if (size(Model,2)<size(Model,1))
    mTranspose=true;
    dim=size(Model,2);
    M=size(Model,1);
else
    mTranspose=false;
    dim=size(Model,1);
    M=size(Model,2);
end

if (size(Data,2)<size(Data,1))
    Data=Data';
end

if dim~=size(Data,1)
    error('The dimension of the Model points and Data points must be equal');
end

N=size(Data,2);

if isempty(TR)
    % Initiate Rotation
    fprintf('initial Rot -> Identity matrix\n')
    TR=eye(dim);
end
if isempty(TT)
    % Initiate Translation
    fprintf('init Tran -> zero matrix')
    TT=zeros(dim,1);
end

% Create closest point search structure


if mTranspose
    DT=delaunayTriangulation(Model);
else
    DT=delaunayTriangulation(Model');
end

% Start the ICP algorithm

res=9e99;

for iter=1:ICP_iter
    
    oldres=res;
    
    % Find closest Model points to Data points
        [vi,resid] = nearestNeighbor(DT,Data');
    % Find transformation
   res=sum(resid.^2);
            
            med=mean(Data,2);
            if mTranspose
                mem=mean(Model(vi,:),1);
                C=Data*Model(vi,:)-(N*med)*mem;
                [U,~,V]=svd(C);
                Ri=V*U';
                if det(Ri)<0
                    V(:,end)=-V(:,end);
                    Ri=V*U';
                end
                Ti=mem'-Ri*med;
            else
                mem=mean(Model(:,vi),2);
                C=Data*Model(:,vi)'-(N*med)*mem';
                [U,~,V]=svd(C);
                Ri=V*U';
                if det(Ri)<0
                    V(:,end)=-V(:,end);
                    Ri=V*U';
                end
                Ti=mem-Ri*med;
            end
    Data=Ri*Data;                       % Apply transformation
    for i=1:dim
        Data(i,:)=Data(i,:)+Ti(i);      %
    end
    
    TR=Ri*TR;                           % Update transformation
    TT=Ri*TT+Ti;                        %
        if abs(oldres-res) < ICP_thres
            break
        end
    
    
end

