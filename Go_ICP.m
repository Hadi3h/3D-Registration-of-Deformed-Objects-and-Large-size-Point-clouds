% Go-ICP implementation in Matlab

function [] = main ()
global MaxRotLevel;
global trim_fraction;
global do_trim;
global Data;
global Model;
global DT_Model;
global Opt_error;

%initial optimal error is a large number
Opt_error = 1e5;
MaxRotLevel = 20;
trim_fraction = 0;
do_trim = 1;

global ICP_iter;
global ICP_thres;
ICP_thres =  1e-5;
ICP_iter = 10000;

%Reading data points and model points
Model =read_plyFile ('CAD_Model.ply');
Data = read_plyFile ('CAD_Data.ply');
DT_Model=delaunayTriangulation(Model);

[OptRotMatrix, OptTransMatrix, Rot_Uncer_radii, Qr] = initialization ();

[ OptRotMatrix, OptTranMatrix]= outBnB(OptRotMatrix, OptTransMatrix, Rot_Uncer_radii, Qr);

fprintf('Optimal Error = %d\n', Opt_error);
disp (OptRotMatrix);
Translated_Data = OptRotMatrix*Data'+repmat(OptTranMatrix,1,size(Data,1));
Translated_Data = Translated_Data';
    
[TR_icp,TT_icp, Error, ICP_data] = icp([], [], Data, Model);
ICP_data = ICP_data';
icp_Data = TR_icp*Data'+repmat(TT_icp,1,size(Data,1));
icp_Data = icp_Data';
    

subplot(3,1,1)
scatter3(Data(:,1),Data(:,2),Data(:,3));
hold on 
scatter3(Model(:,1),Model(:,2),Model(:,3));
legend('Data','Model');

subplot(3,1,2)

scatter3(icp_Data(:,1),icp_Data(:,2),icp_Data(:,3),'*');
hold on
scatter3(Model(:,1),Model(:,2),Model(:,3));
legend('ICP','Model');

subplot(3,1,3)
scatter3(Translated_Data(:,1), Translated_Data(:,2), Translated_Data(:,3),50,'s');
hold on 
scatter3(Model(:,1),Model(:,2),Model(:,3));
legend('Go_ICP','Model');
savefig('Final_regist');

end


% Outer BnB is defined to find the optimal rotation and optimal translation
% matrices!
function [ OptRotMatrix, OptTransMatrix ]= outBnB( OptRotMatrix, OptTransMatrix, Rot_Uncer_radii, Qr)
global SSEThresh;
global Data;
global Opt_error;
global Model;

Nd = size(Data,1);
    [TR,TT,  Error, data_icp] = icp (OptRotMatrix, OptTransMatrix, Data, Model);
    if(Error < Opt_error)
        OptRotMatrix = TR;
        OptTransMatrix = TT;
        fprintf('Error by ICP is better than the current optimal error; Current optimal error =%d, ICP error =%d\n', Opt_error, Error);
		Opt_error = Error;
    end 

   iteration = 0;
    while(1)
       
        if isempty(Qr)
            fprintf('%s','the rotation queue is empty now!\n');
        end
        
        ParentRot = Qr(1);
        Qr(1)= [];
        
        if Opt_error-ParentRot.lb < SSEThresh
            fprintf('Termination: optimal error = %d with pre-defined epsillon%d\n', Opt_error, SSEThresh);
            break;
        end
        iteration =iteration+ 1;
        fprintf('Outer Iteration: %i, So-far Optimal Error = %d \n', iteration, Opt_error);
            % Define subcubes
            sbcube_wdth = ParentRot.w/2;
            Level_subcub = ParentRot.l+1;
            for i=1:8 % 2^3=8
%                 fprintf('\t******Cube %i*****\n', i)
                %NodeRot =struct( 'stx', {},'sty', {},'stz', {}, 'w', {},'l', {}, 'ub', {}, 'lb', {});
                center_subcub = zeros(3,1);
                if i==2 || i==4  || i==6 || i==8
                    strtx = ParentRot.stx + sbcube_wdth;
                elseif i==3 || i ==4 ||i ==7 || i==8
                    strty = ParentRot.sty+ sbcube_wdth;
                elseif i==5 || i==6 || i==7 || i==8
                    strtz = ParentRot.stz+ sbcube_wdth;
                elseif i==1
                    strtx = ParentRot.stx; strty = ParentRot.sty ; strtz = ParentRot.stz;
                end
                center_subcub(1) = strtx + sbcube_wdth/2;
                center_subcub(2) = strty + sbcube_wdth/2;
                center_subcub(3) = strtz + sbcube_wdth/2;
                norm_center = norm (center_subcub);
                % Skip subcube if it is completely outside the rotation PI-ball 
                if norm_center-(sqrt(3)*sbcube_wdth/2) > pi
                    continue;
                end
                % Equation 3: for converting from axis-angle into a rotation
                % matrix
                Rot = zeros(3,3);
                if(norm_center > 0)
                    center_subcub = center_subcub/norm_center;
                	
                    ct = cos(norm_center);
                    ct2 = 1 - ct;
                    st = sin(norm_center);
                    tmp121 = center_subcub(1)*center_subcub(2)*ct2; tmp122 = center_subcub(3)*st;
                    tmp131 = center_subcub(1)*center_subcub(3)*ct2; tmp132 = center_subcub(2)*st;
                    tmp231 = center_subcub(2)*center_subcub(3)*ct2; tmp232 = center_subcub(1)*st;

                    Rot(1,1) = ct + center_subcub(1)*center_subcub(1)*ct2;		Rot(1,2) = tmp121 - tmp122;		Rot(1,3) = tmp131 + tmp132;
                    Rot(2,1) = tmp121 + tmp122;		Rot(2,2) = ct + center_subcub(2)*center_subcub(2)*ct2;		Rot(2,3) = tmp231 - tmp232;
                    Rot(3,1) = tmp131 - tmp132;		Rot(3,2) = tmp231 + tmp232;		Rot(3,3) = ct + center_subcub(3)*center_subcub(3)*ct2;
                
                    %Rotating all Data points with the current subcube's Rotation matrix
                    %(Rot). Data_temp keeps the rotated data points (Nd*3)
                    RData_temp = (Rot*Data')';
                
                else
                    RData_temp = Data;
                end
                cmptLB_flag = 0;%computing Rotation upper bound. 
%                 Setting uncertainty radius to zero in innerBnB to compute Rotation upper bound
%                 fprintf('UB computation\n');
                [ub_Rot, optNodeTrans] = innerBnB(RData_temp, zeros(Nd,1),cmptLB_flag ); % Upper bound of rotation
                
%                 fprintf('Upper bound = %d\n', ub_Rot);
                
                if(ub_Rot < Opt_error)
			
                    % Update optimal error and rotation/translation nodes
                    Opt_error = ub_Rot;
                    OptRotMatrix = Rot;
                    OptTransMatrix (1,1) = optNodeTrans(1).stx+ (optNodeTrans(1).w)/2;
                    OptTransMatrix (2,1) = optNodeTrans(1).sty+ (optNodeTrans(1).w)/2;
                    OptTransMatrix (3,1) = optNodeTrans(1).stz+ (optNodeTrans(1).w)/2;

                	% Run ICP with from OptRotMatrix, OptTransMatrix
                    disp (OptRotMatrix);
                    disp(OptTransMatrix);
                    [TR_icp,TT_icp, Error, Icp_data] = icp(OptRotMatrix, OptTransMatrix, Data, Model);
                    fprintf('\t############error rate with ICP %d#####################\n', Error);
%                     Error = mean(resid.^2);
                    if(Error < Opt_error)
                        fprintf('\t----> optimal Error updated\n')
                        Opt_error = Error;
                        OptRotMatrix = TR_icp;
                        OptTransMatrix = TT_icp;
                        if isequal(OptRotMatrix, TR_icp)
                            fprintf('\t OptRot is equal to TR_icp\n');
                            
                        end

                    end
                    % remove the cubes from Qr that have lower bound larger
                    % than the optimal error
                    New_Qr =struct( 'stx', {},'sty', {},'stz', {}, 'w', {},'l', {}, 'ub', {}, 'lb', {});
                    BF_Qrsize = length(Qr);
                    while(~isempty(Qr))
                        node = Qr(1);
                        Qr(1)=[];
                        if(node.lb < Opt_error)
                            end_New_Qr = length(New_Qr);
                            New_Qr(end_New_Qr+1)=node;
                        else
                            break;
                        end
                    end
                    fprintf('\t reduced Qr length by %d\n',BF_Qrsize-length(New_Qr));
                    Qr = New_Qr;
                end
                cmptLB_flag = 1;
%                 fprintf('LB computation in %i \n',Level_subcub);
                [lb_Rot,] = innerBnB(RData_temp, Rot_Uncer_radii(Level_subcub,:)',cmptLB_flag );  %Translation Node
               

			    %If the best error so far is less than the lower bound, remove the rotation subcube from the queue
                if(lb_Rot >= Opt_error)
                    %Discard the rotation node
                    continue;
                end

                %Update node and put it in queue
                End_Qr = length (Qr);
                Qr(End_Qr+1).stx = strtx;
                Qr(End_Qr+1).sty = strty;
                Qr(End_Qr+1).stz = strtz;
                Qr(End_Qr+1).w = sbcube_wdth;
                Qr(End_Qr+1).l = Level_subcub;
                Qr(End_Qr+1).lb = lb_Rot;
                Qr(End_Qr+1).ub = ub_Rot;
                Qr = nestedSortStruct(Qr, 'lb');
            
            end
                     
                 
    end% until Qr is empty or obtaining the optimal rotation and translation according to |lb_Rot-EOptimal|<epsillon
end



function [optErrorT, nodeTransOut] = innerBnB(RData, RotUncerradi_L, cmptLB_flag)
    global SSEThresh;
    global inlierNum;
    global do_trim;
    global DT_Model
    global Opt_error;

    Qt_init = struct(...
        'stx', {-1}, ...
        'sty', {-1}, ...
        'stz', {-1}, ...
        'w', {2}, ...
        'ub', {1000}, ...
        'lb', {0});
    
  
    nodeTransOut =struct( 'stx', {},'sty', {},'stz', {}, 'w', {});
    % Set optimal translation error to overall so-far optimal error
	% Investigating translation nodes that are sub-optimal overall is redundant
	optErrorT = Opt_error;
	Qt = Qt_init ;
    iteration=0;
	while(1)
        
		if(isempty(Qt))
%             fprintf('iteration %i empty Qt\n', iteration);
			break;
        end
        %Extract the first element from Qt
		nodeTransParent = Qt(1);
        %Remove the first element from Qt
        Qt(1)= [];

		if(optErrorT - nodeTransParent.lb < SSEThresh)
%             fprintf('Convergence inner BnB; iteration %i \t |E*-E_lb| %d\t SSEThreshold %d\n ',iteration, optErrorT-nodeTransParent.lb, SSEThresh);
			break;
        end
        iteration = iteration+1;
        % Create sub cubes from nodeTransParent 
		sbcube_wdth = nodeTransParent.w/2;
		gamma_t = sqrt(3)*sbcube_wdth/2; % Eq.8

		for i=1:8
		
                centerTrans = zeros(3,1);
                if i==2 || i==4  || i==6 || i==8
                    strtx = nodeTransParent.stx + sbcube_wdth;
                elseif i==3 || i ==4 ||i ==7 || i==8
                    strty = nodeTransParent.sty+ sbcube_wdth;
                elseif i==5 || i==6 || i==7 || i==8
                    strtz = nodeTransParent.stz+ sbcube_wdth;
                elseif i==1
                    strtx = nodeTransParent.stx; strty = nodeTransParent.sty ; strtz = nodeTransParent.stz;
                end
                % center of subcube
                centerTrans(1) = strtx + sbcube_wdth/2;
                centerTrans(2) = strty + sbcube_wdth/2;
                centerTrans(3) = strtz + sbcube_wdth/2;
                                      						
			% For each data point, calculate the distance to it's closest point in the model cloud
			% Find distance between transformed point and closest point in model set ||R_r0 * x + t0 - y||
			% RData is the data points rotated by R0
            
            TranslPointData = RData+ repmat(centerTrans',size(RData,1),1);
            DisMatrix = minDistance(TranslPointData, DT_Model);
            if(cmptLB_flag == 1)
					DisMatrix =DisMatrix - RotUncerradi_L; % Eq.26
            end
                DisMatrix = max (DisMatrix, 0); 
			if(do_trim==1) 
                DisMatrix = sort(DisMatrix);
            end
			% For each data point, find the upper and lower bounds
            % (Eq. 26)
            ub = sum(DisMatrix(1:inlierNum).^2);
%             fprintf('Qt: UB %d\n', ub);
			% If upper bound is better than best, update optErrorT and optTransOut (t*)
            if(ub < optErrorT)
                	optErrorT = ub;
                
%                 if(cmptLB_flag==0)

                    if ~isempty(nodeTransOut)
                        nodeTransOut = [];
                    end
					nodeTransOut(1).stx = strtx;
                    nodeTransOut(1).sty = strty;
                    nodeTransOut(1).stz = strtz;
                    nodeTransOut(1).w = sbcube_wdth;
%                     nodeTransOut(1).ub = ub;
                   % nodeTransOut(1).lb = lb;
%                 end
            end
            %compute lower bound for C_t with r0, t0 (line 12 Algorithm 2)
            lb = sum(max((DisMatrix (1:inlierNum)-gamma_t),0).^2); % eq 27
            
			% Remove subcube from queue if lb is bigger than optErrorT
			if(lb < optErrorT)
                %Update node and put it in queue
%                 fprintf('Update Qt: length %i\n', length(Qt));
                End_Qt = length (Qt);
                Qt(End_Qt+1).stx = strtx;
                Qt(End_Qt+1).sty = strty;
                Qt(End_Qt+1).stz = strtz;
                Qt(End_Qt+1).w = sbcube_wdth;
                Qt(End_Qt+1).lb = lb;
                Qt(End_Qt+1).ub = ub;
                Qt = nestedSortStruct(Qt, 'lb');
                %fprintf('length Qt %i\n', length(Qt));
            end
        end
    end % end while (1) 
    end


function loadeddata = read_plyFile (filename)
ptCloud = pcread(filename);
loadeddata = double (ptCloud.Location);
end

function [ P ] = read_txtFile( filename )

file = fopen(filename, 'r');
N = fscanf(file, '%d', 1);
P = fscanf(file, '%f%f%f', [3,N]);
fclose(file);

end


