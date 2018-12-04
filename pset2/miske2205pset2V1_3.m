% Jacob Miske
% 22.05 Pset 2
% Fall 2018

%All neutrons are born at 0,0,0. Medium is homogeneous. Only S and A.
%Cross sections independent of energy (monoenergetics)

%Clear variables
clc; clear all
%start time
tic
%Set problem constants
rng(1) %Select a random seed
neutrons = 1000000; %1e6; %number of neutrons to simulate
%Given cross section constants
waterSigmaA=0.030; waterSigmaS=0.270; waterSt=waterSigmaA+waterSigmaS;
graphiteSigmaA=0.0025; graphiteSigmaS=0.050; graphiteSt=graphiteSigmaA+graphiteSigmaS;
D2OSigmaA=0.0025;D2OSigmaS=0.025; D2OSt=D2OSigmaA+D2OSigmaS;
%Setting boundaries of computation
plottedN=100; scatLimit=200;
%Create placeholders 


%% Run for water

for i=1:neutrons
    %Generate Angles of initial direction
    phi=2*pi*rand(); cosNz=(1-2*rand()); sinNz=sqrt(1-cosNz^2);
    cosNx=sinNz*cos(phi); cosNy=sinNz*sin(phi);
    distanceX=0.0; distanceY=0.0; distanceZ=0.0; %Start at zero
    
    for j=2:1:scatLimit
        %For each neutron up until scatter limit, calculate some random
        %distance until next interaction with medium. Define distance in x,
        %y, and z by random incidence angles.
        distanceTravel=-log(rand())/waterSt;
        distanceX=distanceX+cosNx*distanceTravel; 
        distanceY=distanceY+cosNy*distanceTravel; 
        distanceZ=distanceZ+cosNz*distanceTravel;
        if i<=neutrons
            x(j,i)=distanceX;y(j,i)=distanceY;z(j,i)=distanceZ;
        end
        if rand() < real(waterSigmaS/waterSt)
            phi=2*pi*rand(); cosNz=(1-2*rand()); sinNz=sqrt(1-cosNz^2);
            cosNx=sinNz*cos(phi); cosNy=sinNz*sin(phi);
        else
            break;
        end
    end
end

% % Plot location of water points in 3D
figure(1) %First fig
for i = 1:(plottedN-2) %run through all events (s and eventually a)
    plot3(x(:,i),y(:,i),z(:,i))
    hold on; grid on
end
%Create x at end of neutron travel
for i =(plottedN-1)
    plot3(x(:,i), y(:,i), z(:,i), '-bX')
end
title('100 Neutron Paths in Light Water 3D')
xlabel('X (cm)');ylabel('Y (cm)');zlabel('Z (cm)'); %cross section units are (cm^-1)
saveas(gcf,'Light Water Neutron Path 3D Projection.pdf')

%% Plot location of water points in 2D
figure(4);
for n=1:plottedN;
    plot(x(:,n),y(:,n)); hold on; grid on
end
xlabel('X')
ylabel('Y')
title('100 Light Water Neutron Paths on 2D Projection')
saveas(gcf,'Light Water Neutron Path on 2D Projection.pdf')


%% %% Table Calculations for water
%Max number of collisions
waterMaxCol = size(x,1);
waterX = x; waterY = y; waterZ = z;

%Mean and Max # of collisions to absorption
%Assume neutron doesn't ever return to 0,0,0
collisionCounter = 0;
collisionsPerN = [];
%Relate each point in trace to previous point
distanceBetweenPoints = 0; 
distanceToOrigin = 0; originPoint=[0,0,0]
distanceBetweenN = [];
distanceFromOriginTracker = [];
% Reminds which point is under watch
currentPoint=[]; previousPoint=[]; twoPoints=[];

%Run through waterX, ...Y, ...Z
for j = 1:size(waterX,2) %1:100 number of saved neutrons
    collisionCounter=0; %reset iterator
    
    for i = 2:size(waterX,1) %2:37 (37 happens to be max(col.) -1 (origin))
        if waterX(i,j) ~= 0
            collisionCounter=collisionCounter+1;
        end
        %Define point at each spot in arrays
        currentPoint=[waterX(i,j),waterY(i,j),waterZ(i,j)]; previousPoint=[waterX(i-1,j),waterY(i-1,j),waterZ(i-1,j)];
        twoPoints=[currentPoint;previousPoint];
        %Use pdist function to figure distance along tracks
        distanceBetweenPoints = pdist(twoPoints);
        distanceToOrigin = pdist([currentPoint;originPoint]);
        %Fill arrays    
        distanceBetweenN(i-1,j) = distanceBetweenPoints;
        distanceFromOriginTracker(i-1,j) =distanceToOrigin;
        
    end
    %Save collisionCounter
    collisionsPerN(j) = (collisionCounter);
end
%Using written lists
waterMeanCollisions = mean(collisionsPerN);
waterMaxCollisions = max(collisionsPerN);

%% Calculating mean distance for the individual points, and for each
%neutron's path
%Runs a find function for non-zeros, only accounts them for each column
[~,ii,v] = find(distanceBetweenN);
waterMeanPntPntDistancePerN = accumarray(ii,v,[],@mean);
waterMeanPntPntDistance = mean(waterMeanPntPntDistancePerN);

%%  Derive 1/6 "mean square crow flight distance"
%Take list of distances from origin, square them, and find the mean (*1/6)
waterMeanSquareCrowFlight = max(distanceFromOriginTracker);
waterMeanSquareCrowFlight = waterMeanSquareCrowFlight*1/6;

%% Mean distance travelled between birth and absorption
%Runs a find function for non-zeros, only accounts them for each column
[~,ii,v] = find(distanceFromOriginTracker);
waterMeanPntOriginDistancePerN = accumarray(ii,v,[],@mean);
waterMeanPntOriginDistance = mean(waterMeanPntPntDistancePerN);


%% Maximum distance of any neutron during flight (cm)
waterMaxFromOriginDistancePerN = max(distanceFromOriginTracker);
waterMaxFromOriginDistance = max(waterMaxFromOriginDistancePerN);



%% Run for graphite

for i=1:neutrons
    %Generate Angles of initial direction
    phi=2*pi*rand(); cosNz=(1-2*rand()); sinNz=sqrt(1-cosNz^2);
    cosNx=sinNz*cos(phi); cosNy=sinNz*sin(phi);
    distanceX=0.0; distanceY=0.0; distanceZ=0.0; %Start at zero
    
    for j=2:1:scatLimit
        %For each neutron up until scatter limit, calculate some random
        %distance until next interaction with medium. Define distance in x,
        %y, and z by random incidence angles.
        distanceTravel=-log(rand())/waterSt;
        distanceX=distanceX+cosNx*distanceTravel; 
        distanceY=distanceY+cosNy*distanceTravel; 
        distanceZ=distanceZ+cosNz*distanceTravel;
        if i<=plottedN
            x(j,i)=distanceX;y(j,i)=distanceY;z(j,i)=distanceZ;
        end
        if rand() < real(graphiteSigmaS/graphiteSt)
            phi=2*pi*rand(); cosNz=(1-2*rand()); sinNz=sqrt(1-cosNz^2);
            cosNx=sinNz*cos(phi); cosNy=sinNz*sin(phi);
        else
            break;
        end
    end
end

% Plot location of graphite points 
figure(2) %second fig
for i = 1:(size(x,2)-1) %run through all events (s and eventually a)
    plot3(x(:,i),y(:,i),z(:,i))
    hold on; grid on
end
%Create x at end of neutron travel
for i =(size(x,2))
    plot3(x(:,i), y(:,i), z(:,i), '-bX')
end
title('100 Neutron Paths in Graphite 3D')
xlabel('X (cm)');ylabel('Y (cm)');zlabel('Z (cm)'); %cross section units are (cm^-1)
saveas(gcf,'Graphite Neutron Path 3D Projection.pdf')

%% Plot location of graphite points in 2D
figure(5);
for n=1:plottedN;
    plot(x(:,n),y(:,n)); hold on; grid on
end
xlabel('X (cm)')
ylabel('Y (cm)')
title('100 Graphite Neutron Paths on 2D Projection')
saveas(gcf,'Graphite Neutron Path on 2D Projection.pdf')

%% %% Table Calculations for graphite
%Max number of collisions
graphiteMaxCol = size(x,1);
graphiteX = x; graphiteY = y; graphiteZ = z;

%Mean and Max # of collisions to absorption
%Assume neutron doesn't ever return to 0,0,0
collisionCounter = 0;
collisionsPerN = [];
%Relate each point in trace to previous point
distanceBetweenPoints = 0; 
distanceToOrigin = 0; originPoint=[0,0,0]
distanceBetweenN = [];
distanceFromOriginTracker = [];
% Reminds which point is under watch
currentPoint=[]; previousPoint=[]; twoPoints=[];

%Run through graphiteX, ...Y, ...Z
for j = 1:size(waterX,2) %1:100 number of saved neutrons
    collisionCounter=0; %reset iterator
    
    for i = 2:size(graphiteX,1) %happens to be max(col.) -1 (origin))
        if graphiteX(i,j) ~= 0
            collisionCounter=collisionCounter+1;
        end
        %Define point at each spot in arrays
        currentPoint=[graphiteX(i,j),graphiteY(i,j),graphiteZ(i,j)]; previousPoint=[graphiteX(i-1,j),graphiteY(i-1,j),graphiteZ(i-1,j)];
        twoPoints=[currentPoint;previousPoint];
        %Use pdist function to figure distance along tracks
        distanceBetweenPoints = pdist(twoPoints);
        distanceToOrigin = pdist([currentPoint;originPoint]);
        %Fill arrays    
        distanceBetweenN(i-1,j) = distanceBetweenPoints;
        distanceFromOriginTracker(i-1,j) =distanceToOrigin;
        
    end
    %Save collisionCounter
    collisionsPerN(j) = (collisionCounter);
end
%Using written lists
graphiteMeanCollisions = mean(collisionsPerN);
graphiteMaxCollisions = max(collisionsPerN);

%% Calculating mean distance for the individual points, and for each
%neutron's path
%Runs a find function for non-zeros, only accounts them for each column
[~,ii,v] = find(distanceBetweenN);
graphiteMeanPntPntDistancePerN = accumarray(ii,v,[],@mean);
graphiteMeanPntPntDistance = mean(graphiteMeanPntPntDistancePerN);

%%  Derive 1/6 "mean square crow flight distance"
%Take list of distances from origin, square them, and find the mean (*1/6)
graphiteMeanSquareCrowFlight = max(distanceFromOriginTracker);
graphiteMeanSquareCrowFlight = waterMeanSquareCrowFlight*1/6;

%% Mean distance travelled between birth and absorption
%Runs a find function for non-zeros, only accounts them for each column
[~,ii,v] = find(distanceFromOriginTracker);
graphiteMeanPntOriginDistancePerN = accumarray(ii,v,[],@mean);
graphiteMeanPntOriginDistance = mean(graphiteMeanPntPntDistancePerN);


%% Maximum distance of any neutron during flight (cm)
graphiteMaxFromOriginDistancePerN = max(distanceFromOriginTracker);
graphiteMaxFromOriginDistance = max(graphiteMaxFromOriginDistancePerN);



%% Run for heavy water
for i=1:neutrons
    %Generate Angles of initial direction
    phi=2*pi*rand(); cosNz=(1-2*rand()); sinNz=sqrt(1-cosNz^2);
    cosNx=sinNz*cos(phi); cosNy=sinNz*sin(phi);
    distanceX=0.0; distanceY=0.0; distanceZ=0.0; %Start at zero
    
    for j=2:1:scatLimit
        %For each neutron up until scatter limit, calculate some random
        %distance until next interaction with medium. Define distance in x,
        %y, and z by random incidence angles.
        distanceTravel=-log(rand())/waterSt;
        distanceX=distanceX+cosNx*distanceTravel;
        distanceY=distanceY+cosNy*distanceTravel;
        distanceZ=distanceZ+cosNz*distanceTravel;
        if i<=plottedN
            x(j,i)=distanceX;y(j,i)=distanceY;z(j,i)=distanceZ;
        end
        if rand() < real(D2OSigmaS/D2OSt)
            phi=2*pi*rand(); cosNz=(1-2*rand()); sinNz=sqrt(1-cosNz^2);
            cosNx=sinNz*cos(phi); cosNy=sinNz*sin(phi);
        else
            break;
        end
    end
end


% % Plot location of heavy water points
figure(3) %third fig
for i = 1:(size(x,2)-1) %run through all events (s and eventually a)
    plot3(x(:,i),y(:,i),z(:,i))
    hold on; grid on
end
%Create x at end of neutron travel
for i =(size(x,2))
    plot3(x(:,i), y(:,i), z(:,i), '-bX')
end
title('100 Neutron Paths in Heavy Water 3D')
xlabel('X (cm)');ylabel('Y (cm)');zlabel('Z (cm)'); %cross section units are (cm^-1)
saveas(gcf,'Heavy Water Neutron Path 3D Projection.pdf')

%% Plot location of graphite points in 2D
figure(6);
for n=1:plottedN;
    plot(x(:,n),y(:,n)); hold on; grid on
end
xlabel('X (cm)')
ylabel('Y (cm)')
title('100 Heavy Water Neutron Paths on 2D Projection')
saveas(gcf,'Heavy Water Neutron Path 2D Projection.pdf')


%% %% Table Calculations for heavy water
%Max number of collisions
D2OMaxCol = size(x,1);
D2OX = x; D2OY = y; D2OZ = z;

%Mean and Max # of collisions to absorption
%Assume neutron doesn't ever return to 0,0,0
collisionCounter = 0;
collisionsPerN = [];
%Relate each point in trace to previous point
distanceBetweenPoints = 0; 
distanceToOrigin = 0; originPoint=[0,0,0]
distanceBetweenN = [];
distanceFromOriginTracker = [];
% Reminds which point is under watch
currentPoint=[]; previousPoint=[]; twoPoints=[];

%Run through D2OX, ...Y, ...Z
for j = 1:size(D2OX,2) %1:100 number of saved neutrons
    collisionCounter=0; %reset iterator
    
    for i = 2:size(D2OX,1) %happens to be max(col.) -1 (origin))
        if D2OX(i,j) ~= 0
            collisionCounter=collisionCounter+1;
        end
        %Define point at each spot in arrays
        currentPoint=[D2OX(i,j),D2OY(i,j),D2OZ(i,j)]; previousPoint=[D2OX(i-1,j),D2OY(i-1,j),D2OZ(i-1,j)];
        twoPoints=[currentPoint;previousPoint];
        %Use pdist function to figure distance along tracks
        distanceBetweenPoints = pdist(twoPoints);
        distanceToOrigin = pdist([currentPoint;originPoint]);
        %Fill arrays    
        distanceBetweenN(i-1,j) = distanceBetweenPoints;
        distanceFromOriginTracker(i-1,j) =distanceToOrigin;
        
    end
    %Save collisionCounter
    collisionsPerN(j) = (collisionCounter);
end
%Using written lists
D2OMeanCollisions = mean(collisionsPerN);
D2OMaxCollisions = max(collisionsPerN);

%% Calculating mean distance for the individual points, and for each
%neutron's path
%Runs a find function for non-zeros, only accounts them for each column
[~,ii,v] = find(distanceBetweenN);
D2OMeanPntPntDistancePerN = accumarray(ii,v,[],@mean);
D2OMeanPntPntDistance = mean(D2OMeanPntPntDistancePerN);

%%  Derive 1/6 "mean square crow flight distance"
%Take list of distances from origin, square them, and find the mean (*1/6)
D2OMeanSquareCrowFlight = max(distanceFromOriginTracker);
D2OMeanSquareCrowFlight = D2OMeanSquareCrowFlight*1/6;

%% Mean distance travelled between birth and absorption
%Runs a find function for non-zeros, only accounts them for each column
[~,ii,v] = find(distanceFromOriginTracker);
D2OMeanPntOriginDistancePerN = accumarray(ii,v,[],@mean);
D2OMeanPntOriginDistance = mean(D2OMeanPntPntDistancePerN);


%% Maximum distance of any neutron during flight (cm)
D2OMaxFromOriginDistancePerN = max(distanceFromOriginTracker);
D2OMaxFromOriginDistance = max(D2OMaxFromOriginDistancePerN);


toc