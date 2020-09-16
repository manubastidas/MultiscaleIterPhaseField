%*********************************************************************
% Heuristic macro-scale adaptivity
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [NewActive,NewInactive,newdistances,simplycopy,Micro_dof,stopindic] = ...
    Macro_adaptivity(Sol_Phasefield,Microgeo,u,distances,Active,Macrogeo,...
    Porosity,stopindic)

global Time macro_lambda
global macro_tolr Stop_Poros
% global fine_N

fine_N = 50;
%% Calculate VectSimilarities
dx = 1/fine_N;
x = -0.5+dx/2:dx:0.5-dx/2;
[ym,xm] = meshgrid(x,x);

xx = xm-dx/2;
xx = [xx xx(:,end)];
xx = [xx;0.5*ones(1,size(xx,2))];
yy = ym-dx/2;
yy = [yy;yy(end,:)];
yy = [yy 0.5*ones(size(yy,1),1)];

% [xx,yy]             = meshgrid(0:1/fine_N:1,0:1/fine_N:1);
Finegeo            = struct();
Finegeo.coordinate = [xx(:),yy(:)];
Finegeo.element    = delaunay(xx,yy);
Finegeo.TR         = triangulation(Finegeo.element,Finegeo.coordinate);
% Macrogeo features
[Finegeo]          = myedge(Finegeo);
ref_phase = cell(Macrogeo.nElement,1);
for i=1:Macrogeo.nElement
    [ref_phase{i},~] = Projection_meshes1(Sol_Phasefield{i},Microgeo{i},Finegeo);
end

D = zeros(Macrogeo.nElement);
for i=1:Macrogeo.nElement
    for j=i:Macrogeo.nElement
        D(i,j) = sum(abs(ref_phase{i}.Pres-ref_phase{j}.Pres))*Finegeo.area(i);
        D(j,i) = D(i,j);
    end
end
D = D+ squareform(pdist(u,'cityblock'));
newdistances = exp(-macro_lambda*Time.dt)*distances+Time.dt*(D);
tol_r = macro_tolr*max(max(newdistances));
% tol_e =  0.01*tol_r;
tol_c =  0.2*tol_r;
%
% norm_VectSimilarities = sqrt(sum(u.^2,2));
% [~, maxv] = max(u);
% [~, minv] = min(u);
NewActive = Active; %unique([Active;maxv;minv]);

%% Deactivating
% Dist_a2act = pdist(VectSimilarities(NewActive,:),'Euclidean');
Dist_a2act = newdistances(NewActive,NewActive)+100*eye(length(NewActive),length(NewActive));
% to_deact = setdiff(find(sum(Dist_a2act<tol_c,2)),[maxv;minv]);

to_deact = NewActive(find(sum(Dist_a2act<tol_c,2)));
while ~isempty(to_deact)
    pos = randi(length(to_deact));
    del = to_deact(pos);
    NewActive(NewActive==del) = [];
    Dist_a2act(pos,:) = [];
    Dist_a2act(:,pos) = [];
    %     to_deact = setdiff(find(sum(Dist_a2act<tol_c,2))); %,[maxv;minv]);
    to_deact = NewActive(find(sum(Dist_a2act<tol_c,2)));
end
% NewActive = unique([NewActive;maxv;minv]);

%% Activate - Refinement
NewInactive = setdiff(1:Macrogeo.nElement,NewActive);
Dist_I2act = newdistances(NewInactive,NewActive);
% Dist_I2Iact = newdistances(InActive,InActive);
% Dist_I2act = pdist2(VectSimilarities(InActive,:),VectSimilarities(NewActive,:),'Euclidean');
% Dist_I2Iact = pdist(VectSimilarities(InActive,:),'Euclidean');
% to_act = setdiff(find(sum(Dist_I2act>tol_r,2))); %,[maxv;minv]);
to_act = NewInactive(find(sum(Dist_I2act>tol_r,2)));
while ~isempty(to_act)
    pos = randi(length(to_act));
    add = to_act(pos);
    NewActive = [NewActive;add];
    
    %RELATED INCLUDES ITSELF
    %     related = setdiff(find(newdistances(add,:)<tol_r)); %,[maxv;minv]);
    related = find(newdistances(add,:)<tol_r);
    to_act = setdiff(to_act,related);
    %     Dist_I2Iact(related,:)=[];
    %     Dist_I2Iact(:,related)=[];
end
NewInactive = setdiff(1:Macrogeo.nElement,NewActive);

%% Copy method
Dist_I2act = newdistances(NewInactive,NewActive);
simplycopy = zeros(length(NewInactive),length(NewActive));
for aa = 1:length(NewInactive)
    % nod - number of element
    % aa - number in the list
    %     nod = InActive(aa);
    % min_neig - number in the list
    [~,min_neig]= min(Dist_I2act(aa,:));
    simplycopy(aa,min_neig) = 1;
end

% I N F O R M !!!!
Micro_dof = 0;
% Stop if radii<criteria
for j = 1:length(NewActive)
    others = [NewActive(j); NewInactive(simplycopy(:,j)==1)'];
    Micro_dof = Micro_dof+ Microgeo{NewActive(j)}.nElement;
    if Porosity(NewActive(j)) >= Stop_Poros
        stopindic(others,1) = 1;
    end
end

fprintf('\n Active nodes =  %i \n',length(NewActive))

% Z = zeros(size(Macrogeo.coordinate,1),1);
% figure
% trisurf(Macrogeo.element,Macrogeo.coordinate(:,1),...
%     Macrogeo.coordinate(:,2),Z,...
%     'FaceColor','none')
% hold on
% % colorplot = color(total_actives,:);
% Z = zeros(size(Macrogeo.element,1),1);
% Z(NewActive,1) = 1;
% scatter3(Macrogeo.bari(:,1),...
%     Macrogeo.bari(:,2),Z,80,Z,'filled','MarkerEdgeColor','k')
%
% view([0,90])
% axis tight equal
% grid off
% box on
% pause(5)
% % disp(length(NewActive))

