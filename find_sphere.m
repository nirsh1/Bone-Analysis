function [C0,R]=find_sphere(V,F,C)
%
% Dr. Nir Shvalb, Ariel University, Nirsh@ariel.ac.il
% Bone toolbox July 2016
% 
% usage [C0,R] = find_sphere(V,F,C);
% V, F, C are from an STL file


%%  MASSAGE
    hh = msgbox('select 5 points to fit a sphere');
    child = get(hh,'Children');
    delete(child(2))
    pause(1)
    close(hh)
    clf
    p=plot_bone(V,F,C,0);
    %% select 5 points
    x=0;y=0;z=0; 
    for pts=1:5
            pause
            [P VV VI] = select3d(p);
            hold on;
            plot3(VV(1),VV(2),VV(3),'+k','MarkerSize',20)
            X(pts,1)=P(1);
            X(pts,2)=P(2);
            X(pts,3)=P(3);
    end
    
[Base_X,Base_Y,Base_Z] = sphere(20);
[Center_LSE,Radius_LSE] = sphereFit(X);
 surf(Radius_LSE*Base_X+Center_LSE(1),...
     Radius_LSE*Base_Y+Center_LSE(2),...
     Radius_LSE*Base_Z+Center_LSE(3),'faceAlpha',0.3,'Facecolor','m')

C0=Center_LSE;
    R=Radius_LSE;

end

function [Center,Radius] = sphereFit(X)
% this fits a sphere to a collection of data using a closed form for the
% solution (opposed to using an array the size of the data set). 
% Minimizes Sum((x-xc)^2+(y-yc)^2+(z-zc)^2-r^2)^2
% x,y,z are the data, xc,yc,zc are the sphere's center, and r is the radius

% Assumes that points are not in a singular configuration, real numbers, ...
% if you have coplanar data, use a circle fit with svd for determining the
% plane, recommended Circle Fit (Pratt method), by Nikolai Chernov
% http://www.mathworks.com/matlabcentral/fileexchange/22643

% Input:
% X: n x 3 matrix of cartesian data
% Outputs:
% Center: Center of sphere 
% Radius: Radius of sphere
% Author:
% Alan Jennings, University of Dayton

A=[mean(X(:,1).*(X(:,1)-mean(X(:,1)))), ...
    2*mean(X(:,1).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,1).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    mean(X(:,2).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,2).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    0, ...
    mean(X(:,3).*(X(:,3)-mean(X(:,3))))];
A=A+A.';
B=[mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,1)-mean(X(:,1))));...
    mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,2)-mean(X(:,2))));...
    mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,3)-mean(X(:,3))))];
Center=(A\B).';
Radius=sqrt(mean(sum([X(:,1)-Center(1),X(:,2)-Center(2),X(:,3)-Center(3)].^2,2)));
end