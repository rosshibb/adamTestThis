
% Generate a basic velocity field with several adjacent vortices
[X,Y] = meshgrid(-5:pi/8:5,-5:pi/8:5);
U = sin(Y);
V = cos(X);

% Compute vorticity (can also use the fonction curl)
[dxU, dyU] = gradient(U);
[dxV, dyV] = gradient(V);
vort = dxV - dyU;

% Increase the resolution of the vorticity field by interpolation
smoothness = 0.01; % parameter controlling the resolution of the vorticity field (smaller values = finer mesh)
[Xq,Yq] = meshgrid(min(X(:)):smoothness:max(X(:)),min(Y(:)):smoothness:max(Y(:))); % makes Y and Y data grid with higher resolution
Vq = interp2(X,Y,vort,Xq,Yq); % interpolates vorticity (vort) data from original XY grid to higher resolution grid

% Plot velocity and vorticity field
hLines = 201; % number of division to plot vorticity colormap
vorticity_colormap = cmocean('balance');  % define colormap

figure; hold on
contourf(Xq,Yq,Vq,hLines,'edgecolor','none');   % Plot vorticity fields
colorbar
colormap(vorticity_colormap)
quiver(X, Y, U, V,'k');

% formatFigure

% Increase the resolution of the vorticity field by interpolation
smoothness = 0.1; % parameter controlling the resolution of the vorticity field (smaller values = finer mesh)
[Xq,Yq] = meshgrid(min(X(:)):smoothness:max(X(:)),min(Y(:)):smoothness:max(Y(:))); % makes Y and Y data grid with higher resolution
Vq = interp2(X,Y,vort,Xq,Yq); % interpolates vorticity (vort) data from original XY grid to higher resolution grid

% Plot velocity and vorticity field
hLines = 201; % number of division to plot vorticity colormap
vorticity_colormap = cmocean('balance');  % define colormap

figure; hold on
contourf(Xq,Yq,Vq,hLines,'edgecolor','none');   % Plot vorticity fields
colorbar
colormap(vorticity_colormap)
quiver(X, Y, U, V,'k');
% 
formatFigure
%% Find center of vortex in selected area

% Calculate the velocity gradient tensor components
D11 = dxU; % du/dx
D12 = dyU; % du/dy
D21 = dxV; % dv/dx
D22 = dyV; % dv/dy

% Symmetric part of the velocity gradient tensor (Rate-of-strain tensor)
S11 = (D11 + D11) / 2;
S12 = (D12 + D21) / 2;
S21 = S12;
S22 = (D22 + D22) / 2;

% Antisymmetric part of the velocity gradient tensor (Vorticity tensor)
Omega11 = 0; % By definition, since it's antisymmetric and 2D
Omega12 = (D12 - D21) / 2;
Omega21 = -Omega12;
Omega22 = 0; % By definition

% Compute the invariants of the velocity gradient tensor
P = -(D11 + D22);
Q = 0.5 * ((Omega12 - Omega21).^2 - (S11.^2 + 2.*S12.*S21 + S22.^2));
R = -(D11.*D22 - D12.*D21); %a determinant, but calculated differently because matrices not square
Qcrit = interp2(X,Y,Q,Xq,Yq); % interpolates vorticity (vort) data from original XY grid to higher resolution grid

%%

% Visualization
figure;
contourf(Xq, Yq, Qcrit, 100, 'LineColor', 'none');
colorbar;
title('Q-criterion for Vortex Detection');
xlabel('X');
ylabel('Y');

%%
% Assuming Qcrit and the velocity field (U, V) preparation is done

% Define the range of nr levels to scan
nr = max(Qcrit, [], 'all');
minPercentage = 0.01;
maxPercentage = 0.99;
numSteps = 10; % Number of steps between min and max percentage
levels = linspace(minPercentage * nr, maxPercentage * nr, numSteps);
levelsPlot = linspace(minPercentage, maxPercentage, numSteps);
circulations = zeros(size(levels));

% Plot the Q-criterion for user to select a vortex
figure; contourf(Xq, Yq, Qcrit, 100, 'LineColor', 'none');
colorbar;
title('Click near the vortex of interest');
[x_click, y_click] = ginput(1); % User clicks near a vortex

for i = 1:length(levels)
    currentLevel = levels(i);
    % Update binaryImage for current level
    binaryImage = Qcrit > currentLevel;
    
    figure
    imshow(binaryImage)
    drawnow;

    % Find boundaries at the current level
    boundaries = bwboundaries(binaryImage);
    
    % Initialize minimum distance to a large number
    minDist = inf;
    selectedBoundary = [];
    
    for k = 1:length(boundaries)
        % For each boundary, find the one closest to the click
        boundary = boundaries{k};
        boundaryX = Xq(sub2ind(size(Xq), boundary(:,1), boundary(:,2)));
        boundaryY = Yq(sub2ind(size(Yq), boundary(:,1), boundary(:,2)));
        distances = sqrt((boundaryX - x_click).^2 + (boundaryY - y_click).^2);
        [minDistK, ~] = min(distances);
        
        if minDistK < minDist
            minDist = minDistK;
            selectedBoundary = boundary;
        end
    end
    
    % Check if a boundary was selected
    if ~isempty(selectedBoundary)
        % Convert boundary indices to coordinates
        boundaryX = Xq(sub2ind(size(Xq), selectedBoundary(:,1), selectedBoundary(:,2)));
        boundaryY = Yq(sub2ind(size(Yq), selectedBoundary(:,1), selectedBoundary(:,2)));
        
        % Approximate circulation calculation for the selected boundary
        circulation = 0;
        for j = 1:length(boundaryX) - 1
            midX = (boundaryX(j) + boundaryX(j+1)) / 2;
            midY = (boundaryY(j) + boundaryY(j+1)) / 2;
            
            velocityX = interp2(X, Y, U, midX, midY, 'linear');
            velocityY = interp2(X, Y, V, midX, midY, 'linear');
            
            dxcirc = boundaryX(j+1) - boundaryX(j);
            dycirc = boundaryY(j+1) - boundaryY(j);
            
            circulation = circulation + velocityX * dxcirc + velocityY * dycirc;

        end
        
        % Close the loop for the circulation calculation
        midX = (boundaryX(end) + boundaryX(1)) / 2;
        midY = (boundaryY(end) + boundaryY(1)) / 2;
        velocityX = interp2(X, Y, U, midX, midY, 'linear');
        velocityY = interp2(X, Y, V, midX, midY, 'linear');
        dx = boundaryX(1) - boundaryX(end);
        dy = boundaryY(1) - boundaryY(end);
        circulation = circulation + velocityX * dx + velocityY * dy;

        circulations(i) = circulation;
    else
        circulations(i) = NaN; % Mark as NaN if no boundary is selected
    end
end

% Plotting the circulation values
figure;
plot(levelsPlot, circulations, '-o');
xlabel('Level');
ylabel('Circulation');
title('Circulation across different levels');

