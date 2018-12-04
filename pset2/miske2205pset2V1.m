% Jacob Miske
% 22.05 Pset 2
% Fall 2018

%All neutrons are born at 0,0,0. Medium is homogeneous. Only S and A.
%Cross sections independent of energy (monoenergetics)

%Clear variables
clc; clear all
%Set problem constants
rng(1) %Select a random seed
neutrons = 100; %1e6; %number of neutrons to simulate
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
        if i<=plottedN
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
for i = 1:(size(x,2)-1) %run through all events (s and eventually a)
    plot3(x(:,i),y(:,i),z(:,i))
    hold on
end
%Create x at end of neutron travel
for i =(size(x,1))
    plot3(x(:,i), y(:,i), z(:,i), 'x')
end
title('100 Neutron Paths in Light Water')
xlabel('cm');ylabel('cm');zlabel('cm'); %cross section units are (cm^-1)

%% Plot location of water points in 2D
figure(4);
for n=1:plottedN;
    plot(x(:,n),y(:,n),'r'); hold on
end
xlabel('X')
ylabel('Y')
title('Light Water Neutron Path on 2D Projection')
saveas(gcf,'Light Water Neutron Path on 2D Projection.pdf')


%% Calculations for water
%Max number of collisions
waterMaxCol = size(x,1);
waterX = x; waterY = y; waterZ = z;

%Mean # of collisions to absorption
%Assume neutron doesn't ever return to 0,0,0
collisionCounter = 0;
collisionsPerN = [];

for i = 2:size(waterX,1)
    
    for j = 1:size(waterX,2)
        if waterX(i,j) ~= 0
            collisionCounter=collisionCounter+1;
        end
    end
    collisionsPerN(j) = (collisionCounter);
    collisionCounter = 0;
end

%Calculating mean distance for the individual points
for i = 1:size(waterX,1)
    currentPoint = [0,0,0];
    for j = 1:size(waterX,2)
        currentPoint(1) = waterX(i,j);currentPoint(2) = waterY(i,j);currentPoint(3) = waterZ(i,j);
        
    end
    
end
% Derive 1/6 "mean square crow flight distance"


%Mean distance travelled between birth and absorption


% Maximum distance of any neutron during flight (cm)




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
    hold on
end
%Create x at end of neutron travel
for i =(size(x,2))
    plot3(x(:,i), y(:,i), z(:,i), 'x')
end

title('100 Neutron Paths in Graphite')
xlabel('cm');ylabel('cm');zlabel('cm'); %cross section units are (cm^-1)





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
    hold on
end
%Create x at end of neutron travel
for i =(size(x,2))
    plot3(x(:,i), y(:,i), z(:,i), 'x')
end

title('100 Neutron Paths in Heavy Water')
xlabel('cm');ylabel('cm');zlabel('cm'); %cross section units are (cm^-1)

