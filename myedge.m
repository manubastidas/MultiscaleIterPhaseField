%*********************************************************************
% Auxiliar code 
% Read geometry
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium


%function [nodes2element,inneredge]=edge(element,coordinate);
function [geometry]=edge(geometry)

element = geometry.element;
coordinate = geometry.coordinate; 
% nElement -> number of element at each mesh
geometry.nElement  = size(element,1);
% nnodes -> number of nodes at each mesh
geometry.nnodes    = size(coordinate,1);

nodes2element=sparse(geometry.nnodes,geometry.nnodes);
bari = zeros(geometry.nElement,2);

for j=1:geometry.nElement
    nodes2element(element(j,:),element(j,[2 3 1]))=...
        nodes2element(element(j,:),element(j,[2 3 1]))+j*eye(3,3);
    bari(j,:) = sum(coordinate(element(j,:),:))/3;
    area(j)=polyarea(coordinate(element(j,:),1),coordinate(element(j,:),2));
end
minsize = min(sqrt(sum((coordinate(element(:,[1 2 3]),:)-...
    coordinate(element(:,[2 3 1]),:)).^2,2)));
maxsize = max(sqrt(sum((coordinate(element(:,[1 2 3]),:)-...
    coordinate(element(:,[2 3 1]),:)).^2,2)));

B=nodes2element+nodes2element';
[i,j]=find(triu(B));
nodes2edge=sparse(i,j,1:size(i,1),geometry.nnodes,geometry.nnodes);
nodes2edge=nodes2edge+nodes2edge';
noedges=size(i,1);

% to generate element of edge
edge2element=zeros(size(i,1),4);
for m = 1:geometry.nElement
    for k = 1:3
        initial_edge = element(m,k);
        end_edge = element(m,rem(k,3)+1);
        %     p = nodes2edge(element(m,k),element(m,rem(k,3)+1)); % update on 13/2/3
        p = nodes2edge(initial_edge,end_edge);
        if edge2element(p,1)==0
            edge2element(p,:)=[initial_edge end_edge nodes2element(initial_edge,end_edge) ...
                nodes2element(end_edge,initial_edge)];
            %       edge2element(p,:)=[element(m,k) element(m,rem(k,3)+1)  nodes2element(element(m,k),element(m,rem(k,3)+1)) ...
            %                           nodes2element(element(m,rem(k,3)+1),element(m,k))];
            
        end
    end
end

% To produce the interior edges
%[l,m]=find(edge2element(:,4));
interioredge=edge2element(find(edge2element(:,4)),:);

%[r,s]=find(edge2element(:,4)==0);
% exterioredge=edge2element(find(edge2element(:,4)==0),[1,2,3]);

geometry.nodes2element = nodes2element;
geometry.nodes2edge= nodes2edge;
geometry.noedges= noedges;
geometry.edge2element= edge2element;
geometry.interioredge= interioredge;
geometry.bari= bari;
geometry.area= area;
geometry.minsize= minsize;
geometry.maxsize= maxsize;
kkk = boundary(coordinate(:,1),coordinate(:,2));
geometry.Gamma = [kkk(1:end-1) kkk(2:end)];

