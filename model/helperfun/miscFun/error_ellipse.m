% Purpose:        Get coordinates of error ellipse with confidence interval.
%                 Error ellipse will be oriented parallel to the data.
%
% Adapated from:  https://gist.github.com/Piyush3dB/bf2c83a8eb7344798644
% Adapted by:     Michael Jigo

function [x y] = error_ellipse(avg,covariance,ci)

[eigenvec, eigenval] = eig(covariance);


% Get the index of the largest eigenvector
[max_evc_ind_c, r] = find(eigenval == max(max(eigenval)));
max_evc = eigenvec(:, max_evc_ind_c);

% Get the largest eigenvalue
max_evl = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(max_evc_ind_c == 1)
    min_evl = max(eigenval(:,2));
    min_evc = eigenvec(:,2);
else
    min_evl = max(eigenval(:,1));
    min_evc = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(max_evc(2), max_evc(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end


% Get the 95% confidence interval error ellipse
%chisquare_val = 2.4477;
chisquare_val = sqrt(chi2inv(ci,2));
theta_grid = linspace(0,2*pi);
phi = angle;
X0  = avg(1);
Y0  = avg(2);
a   = chisquare_val*sqrt(max_evl);
b   = chisquare_val*sqrt(min_evl);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Get coordinates of ellipse
x = r_ellipse(:,1) + X0;
y = r_ellipse(:,2) + Y0;
return

figure;
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-')
hold on;

% Plot the eigenvectors
figure;
quiver(X0, Y0, max_evc(1)*sqrt(max_evl), max_evc(2)*sqrt(max_evl), '-m', 'LineWidth',2);
hold on;
quiver(X0, Y0, min_evc(1)*sqrt(min_evl), min_evc(2)*sqrt(min_evl), '-g', 'LineWidth',2);

% Set the axis labels
hXLabel = xlabel('x');
hYLabel = ylabel('y');

