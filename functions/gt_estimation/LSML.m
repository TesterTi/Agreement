function [GT_LSP1, intermediate] = LSML(annotations)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% An implementation of LSML:
%   
%   X. Li, B. Aldridge, R. Fisher, J. Rees, Estimating the Ground Truth
%      from Multiple Individual Segmentations Incorporating Prior Pattern
%      Analysis with Applications to Skin Lesion Segmentation, IEEE 
%      International Symposium on Biomedical Imaging: From Nano to Macro,
%      1438-1441, 2011
%
% Provided by X. Li.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


display = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the voting distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annotations = double(annotations);
region = sum(annotations(:,:,1:size(annotations,3)),3)/size(annotations,3);
% imagesc(region), colormap gray


%figure, imagesc(region), colormap(gray)
%%imshow(((region)),[])
%hold on
%%contour(synMat,'color',[1 0 0])
%contour(GT_compact,'color',[0 1 0])
%%contour(GT_detailed,'color',[0 0 1])
%hold off
%legend('GT', 'GT(compact)','GT(detailed)')
%drawnow




%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the initial GT
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INITIAL IMAGE ORIGINALLY USED %%%
% s1 = double(region>=(size(SEG,3)/2)/size(SEG,3));
% SE = strel('disk', 1, 8);
% initseg = imerode(s1,SE);

%%% INITIAL IMAGE USING 50% AGREEMENT %%%
initseg = region >= 0.5;

if display
    figure, imagesc(initseg), colormap(gray)
    drawnow
end




%%%%%%%
% TEST 
%%%%%%%
prior = ones(1,size(annotations,3));
[GT_LSP1, intermediate] = GToptimization(initseg, double(annotations), 1000, 0.0, display, 'pro', prior, [], region);
%[GT_LSP2] = GToptimization(initseg, double(SEG), 400, 0.0, true, 'oth', prior, synMat, region);





%%%%%%%%%%%%%%%%%
% display result
%%%%%%%%%%%%%%%%%

if display
%figure, imagesc(region), colormap(gray)
%title('region')
figure, imagesc(initseg), colormap(gray)
title('initial')
figure, imagesc(GT_LSP1), colormap(gray)
title('LSMLP')
%figure, imagesc(GT_LSP2), colormap(gray)
%title('result 2')

%hold off
end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [seg, intermediate_results] = GToptimization(InitSeg, ManualSegMat, max_its, alpha, display, approach, prior, synMat, region)
    
    % InitSeg - initial grount truth, might use one of the manual result or
    % the ground truth based on voting policy
    % ManualSegMat - manual segmentation results matrix
  
    
    wanted = [100:100:1000];
    ind = 1;
    intermediate_results = zeros(size(InitSeg,1), size(InitSeg,2), numel(wanted));
    
    %-- default value for parameter alpha is .1
    if(~exist('alpha','var')) 
        alpha = 0.1; 
    end
    
    %-- default behavior is to display intermediate outputs
    %if(~exist('display','var'))
    %  display = true;
    %end
    
    %garma = 1;
    
    %-- ensures image is 2D double matrix
    if numel(size(InitSeg)) ~= 2
        exit
    end
    
    
    %-- Create a signed distance map (SDF) from mask
    phi = mask2phi(double(InitSeg)); % in < 0 ; out > 0
    seg = InitSeg;
    
    
    %--main loop
    for its = 1:max_its   % Note: no automatic convergence test
        
        idx = find(phi <= 1.2 & phi >= -1.2);  %get the curve's narrow band

        if strcmp(approach,'pro')
            [STF]= BayesSTF(seg,ManualSegMat,idx, prior);
        else
            STF = DiffSTF(seg, ManualSegMat,idx);
        end

        NAN = sum(isnan(STF));

        if NAN > 0 || length(idx) ~= length(STF)
            break;
        else     
            curvature = get_curvature(phi,idx);  % force from curvature penalty
            
            %dphidt = alpha*curvature + STF; 
            dphidt = STF./max(abs(STF)) + alpha*curvature./max(abs(curvature)); 
            
            %-- maintain the CFL condition
            dt = .5/(max((dphidt))+eps);
            
            %-- evolve the curve
            phi(idx) = phi(idx) + dt*dphidt;
            
            %-- Keep SDF smooth
            phi = sussman(phi, .5);
            
        end

        %-- make mask from SDF
        seg = phi<=0; %-- Get mask from levelset
        
        if display
            imshow(region,[])
            hold on
            %contour(synMat,'r')    
            %contour(seg,'g')

            %     contour(ManualSegMat(:,:,1),'color',[0 0 1])
            %     contour(ManualSegMat(:,:,1),'color',[0 0 0.8])
            %     contour(ManualSegMat(:,:,1),'color',[0 0 0.6])

            %     contour(ManualSegMat(:,:,2),'color',[0 0 1])
            %     contour(ManualSegMat(:,:,2),'color',[0 0 0.8])
            %     contour(ManualSegMat(:,:,2),'color',[0 0 0.6])
            %     
            %     contour(ManualSegMat(:,:,3),'color',[0 0 1])
            %     contour(ManualSegMat(:,:,3),'color',[0 0 0.8])
            %     contour(ManualSegMat(:,:,3),'color',[0 0 0.6])
            %     
            %     contour(ManualSegMat(:,:,4),'color',[0 0 1])
            %     contour(ManualSegMat(:,:,4),'color',[0 0 0.8])
            %     contour(ManualSegMat(:,:,4),'color',[0 0 0.6])
            %     
            %     contour(ManualSegMat(:,:,5),'color',[0 0 1])
            %     contour(ManualSegMat(:,:,5),'color',[0 0 0.8])
            %     contour(ManualSegMat(:,:,5),'color',[0 0 0.6])
            %     
            %     contour(ManualSegMat(:,:,6),'color',[0 0 1])
            %     contour(ManualSegMat(:,:,6),'color',[0 0 0.8])
            %     contour(ManualSegMat(:,:,6),'color',[0 0 0.6])
            %     
            %     contour(ManualSegMat(:,:,7),'color',[0 0 1])
            %     contour(ManualSegMat(:,:,7),'color',[0 0 0.8])
            %     contour(ManualSegMat(:,:,7),'color',[0 0 0.6])
            %     
            %     contour(ManualSegMat(:,:,8),'color',[0 0 1])
            %     contour(ManualSegMat(:,:,8),'color',[0 0 0.8])
            %     contour(ManualSegMat(:,:,8),'color',[0 0 0.6])

            hold off
            title(num2str(its))
            drawnow
        end
        
        % save some intermediate results
        if any(wanted == its)
            
            intermediate_results(:,:,ind) = phi<=0;
            
            ind = ind + 1;
            
        end
        
    end
    
    %-- make mask from SDF
    seg = phi<=0; %-- Get mask from levelset
%   seg = imfill(seg,'holes');
  
end


%-- converts a mask to a SDF
function phi = mask2phi(init_a)

    phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
    
end


%-- compute curvature along SDF
function curvature = get_curvature(phi,idx)
    [dimy, dimx] = size(phi);        
    [y x] = ind2sub([dimy,dimx],idx);  % get subscripts

    %-- get subscripts of neighbors
    ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

    %-- bounds checking  
    ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
    yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;    

    %-- get indexes for 8 neighbors
    idup = sub2ind(size(phi),yp1,x);    
    iddn = sub2ind(size(phi),ym1,x);
    idlt = sub2ind(size(phi),y,xm1);
    idrt = sub2ind(size(phi),y,xp1);
    idul = sub2ind(size(phi),yp1,xm1);
    idur = sub2ind(size(phi),yp1,xp1);
    iddl = sub2ind(size(phi),ym1,xm1);
    iddr = sub2ind(size(phi),ym1,xp1);

    %-- get central derivatives of SDF at x,y
    phi_x  = -phi(idlt)+phi(idrt);
    phi_y  = -phi(iddn)+phi(idup);
    phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
    phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
    phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
             +0.25*phi(iddr)+0.25*phi(idul);
    phi_x2 = phi_x.^2;
    phi_y2 = phi_y.^2;

    %-- compute curvature (Kappa)
    curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
              (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        
end



      
%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
    % forward/backward differences
    a = D - shiftR(D); % backward
    b = shiftL(D) - D; % forward
    c = D - shiftD(D); % backward
    d = shiftU(D) - D; % forward

    a_p = a;  a_n = a; % a+ and a-
    b_p = b;  b_n = b;
    c_p = c;  c_n = c;
    d_p = d;  d_n = d;

    a_p(a < 0) = 0;
    a_n(a > 0) = 0;
    b_p(b < 0) = 0;
    b_n(b > 0) = 0;
    c_p(c < 0) = 0;
    c_n(c > 0) = 0;
    d_p(d < 0) = 0;
    d_n(d > 0) = 0;

    dD = zeros(size(D));
    D_neg_ind = find(D < 0);
    D_pos_ind = find(D > 0);
    dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
    dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;

    D = D - dt .* sussman_sign(D) .* dD;
end


%-- whole matrix derivatives
function shift = shiftD(M)
    shift = shiftR(M')';
end


function shift = shiftL(M)
    shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];
end


function shift = shiftR(M)
    shift = [ M(:,1) M(:,1:size(M,2)-1) ];
end



function shift = shiftU(M)
    shift = shiftL(M')';
end
  
  
  
function S = sussman_sign(D)
    S = D ./ sqrt(D.^2 + 1);
end
  
  

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 2 - difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function STF = DiffSTF(InitSeg, proMat,idx)

    proMat = proMat*8;

    J = [reshape(proMat(:,:),length(proMat(:)),1)];

    J_idx = J(idx,:);

    STF = 4-J_idx;
    
end





%%%%%%%%%
% STAPLE
%%%%%%%%%
function [seg, precision, specificity] = STAPLE_CODE(SEG)

    SEG = (SEG);
    [M,N,number] = size(SEG);
    precision = [0.6 0.6 0.6 0.6 0.99 0.99 0.99 0.99];
    specificity = [0.9 0.9 0.9 0.9 0.6 0.6 0.6 0.6];
    Problesion = ones([M,N]);
    Probskin = ones([M,N]);
    for j = 1:number
        Problesion = Problesion.*(sign(SEG(:,:,j)).*precision(j)+sign(1-SEG(:,:,j)).*(1-precision(j)));
        Probskin = Probskin.*(sign(1-SEG(:,:,j)).*specificity(j)+sign(SEG(:,:,j)).*(1-specificity(j)));
    end
    seg = double(Problesion>=Probskin);
    % seg = SEG(:,:,1);
    XOR = 1;
    fT1 = sum(SEG,3)/number; % prior probablity of being lesion
    fT2 = 1-fT1;% prior probablity of being skin
    k = 1; 
    S = 100;
    oldS = 110;
    while abs(S-oldS)>eps 

        oldS = S;
        % expectation step
        a = ones(size(seg));
        b = ones(size(seg));
        for j = 1:number
            a = a.*(sign(SEG(:,:,j)).*precision(j)+sign(1-SEG(:,:,j)).*(1-precision(j)));
            b = b.*(sign(1-SEG(:,:,j)).*specificity(j)+sign(SEG(:,:,j)).*(1-specificity(j)));
        end
        a = a.*fT1;
        b = b.*fT2;

        W = a./(a+b+eps);


        % maximization step
        for j = 1:number
            precision(j) = sum(sum(sign(SEG(:,:,j)).*W))/sum(sum(W));
            specificity(j) = sum(sum(sign(1-SEG(:,:,j)).*(1-W)))/sum(sum(1-W));
        end

        S = sum(sum(W));

        k = k+1;

    end
    Problesion = ones(size(seg));
    Probskin = ones(size(seg));
    for j = 1:number
        Problesion = Problesion.*(sign(SEG(:,:,j)).*precision(j)+sign(1-SEG(:,:,j)).*(1-precision(j)));
        Probskin = Probskin.*(sign(1-SEG(:,:,j)).*specificity(j)+sign(SEG(:,:,j)).*(1-specificity(j)));
    end
    seg = double(Problesion>=Probskin);
    seg = imfill(seg,'holes');

end

 







%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 1 - probability 
%%%%%%%%%%%%%%%%%%%%%%%%%

function [STF]= BayesSTF(InitSeg, segMat,idx, prior)
    
    % precision - precison parameter matrix
    % specificity - specification parameter matrix
    segMat = double(segMat);
    [M,N,D] = size(segMat);
    for i = 1:D
        [precision(i), specificity(i)]=PrecANDSpec(segMat(:,:,i),double(InitSeg),0);
    end

    Problesion = ones(size(InitSeg));
    Probskin = ones(size(InitSeg));
    for i = 1:D
        Problesion = Problesion.*(sign(segMat(:,:,i)).*precision(i)+sign(1-segMat(:,:,i)).*(1-precision(i))).^prior(i);
        Probskin = Probskin.*(sign(1-segMat(:,:,i)).*specificity(i)+sign(segMat(:,:,i)).*(1-specificity(i))).^prior(i);
    end 

    J2 = [reshape(Probskin,length(Probskin(:)),1)];
    J1 = [reshape(Problesion,length(Problesion(:)),1)];

    J1_idx = J1(idx,:);
    J2_idx = J2(idx,:);

    STF = log((J2_idx)+eps)-log(J1_idx+eps);

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [precision, specificity] = PrecANDSpec(seg,GT,display)
    
    % input - seg1 (STA)
    %         seg2 (extracted region)

    % calculate precision of two segmentations (reflects the accuracy of the
    % extracted area)
    % seg = mat2gray(seg);
    % GT = mat2gray(GT);
    
    GT = GT(:);
    seg = seg(:);
    combine = seg+GT;
    
    TP = sum(combine==2);
    FP = sum(seg == 1 & GT == 0);
    %FP = length(find(seg==1 & GT == 0));
    TN = sum(combine==0);
    FN = sum(seg == 0 & GT == 1);
    %FN = length(find(seg==0 & GT == 1));

    precision = TP/(TP+FN);
    specificity = TN/(TN+FP);


    if display
        seg = ones(size(seg1));
        seg(find(combine==2)) = 2; % overlap
        seg(find(seg1==1&combine==1))=3; % sta
        seg(find(seg2==1&combine==1))=4; % computer

        colmap=[0.5 0.5 0.5;1 0 0; 0 1 0; 0 0 1];
        Img=ind2rgb(seg,colmap);
        imshow(Img)
    end

end