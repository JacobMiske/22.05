%% Jacob Miske
% 22.05 2015 PSET 1

%Remove written variables
clc; clear;
tic
%% Create initial variables
gamblers = ones(1,10000)*10000; %Each person starts with 10,000
gamblingTimes = 1:1:10000; %Assume each person goes for 10,000 bets/day
gamblerTimeSum = []; %for summing winnings over time

%Let's say if "round(rand)" produces 0, the gambler is wrong at that
%instance of the game

%This for loop is for one day
for i = 1:1:size(gamblingTimes, 2) %for each time possible spent gambling
    for j = 1:1:size(gamblers, 2) %for each gambler in each time slot, create a new 50/50 probability
        Outcome(j) = round(rand); %create the 50/50 chance
        if Outcome(j) == 0 && gamblers(j) >99 %gambler loses and can pay out
            gamblers(j) = gamblers(j) - 100; %If lose, detract
        end
        if Outcome(j) == 1 && gamblers(j) >99 %gambler wins and hasn't gone broke
            gamblers(j) = gamblers(j) + 100; %If win, increase
        end
    end
    %Calculate over-time total gambler 
    gamblerTimeSum.append(i) = sum(gamblers); %run sum of all gamblers money per time unit
end




%% Create an array from the appended structure
gamblerTimeSumArray = struct2array(gamblerTimeSum);
%Create plot of total gambler earnings at each time interval
figure(1)
plot(gamblingTimes, gamblerTimeSumArray, 'b');
savefig('Gambler Earnings over time.fig')
title('Total Gambler Winnings over Time')
xlabel('Times to Gamble (1 = 8.64seconds)')
ylabel('Dollars')


%Probability that some gambler ends day with >$20000
%Check through gamblers array for those who ended with >$20000
greater50kGamblers = 0;
for i = 1:1:10000;
    if gamblers(i) >20000
        greater50kGamblers = greater50kGamblers+1;
    end
end
Pgt20k = greater50kGamblers/10000

%% Question 1: How much cash should casino start with each day to have a 99%
%chance of making necessary payouts, given a fair game? Higher gameOdds hurt casino!
AdvantageFactor=0.99
Advantage = (1+AdvantageFactor)/2

%If game winning probability = 0.5, casino should break even. Assuming 10,000 gamblers each
%making $250 wagers, 10000 times... Assuming nobody goes broke beforehand.
expectedDailyPayout = 0.5*100*10000*10000 - 8000000000 %Large number is for gameOdds = 0.5
potentialDailyPayoutArray = []
AdvantageArray = 0.000:0.1111:0.999 %Different house advantages to examine
for i =1:1:9;
    potentialDailyPayoutArray(i) = ((AdvantageArray(i)+1)/2)*100*10000*10000 - 8000000000;
end
%Plotting potentialDailyPayoutArray to different house advantage factors
figure(2)
plot(AdvantageArray, potentialDailyPayoutArray, 'ro')
title('House Return versus House Advantage Factor')
xlabel('Advantage Factors')
ylabel('Dollars of House Money At Risk (Starting Cash Minimum')

%% Question 2: The casino needs +$20mil for loan payments, what is the
%minimum house advantage the casino needs?
expectedDailyPayoutArray = [];
gameOddsArray = 0.1:0.1:0.9;
%plot expected daily revenue versus house advantage factor
for i = 1:1:9
    expectedDailyPayoutArray(i) = 8000000000 - gameOddsArray(i)*250*10000*10000;
end
figure(3)
plot(AdvantageArray,expectedDailyPayoutArray)
savefig('Expected Revenue.fig')
title('Expected Revenue versus House Advantage')
twentymillion = 1:1:9;
twentymillion = twentymillion*20000000
plot(AdvantageArray,twentymillion)


%% Question 3: Two ounces of vodka per gambler, not broke, per hour

%Calculate once again the main iteration 
AdvantageFactor = 0.001
Advantage = (1+AdvantageFactor)/2;

%Calculate complimentary vodka
%Casino is 24hr, two ounces of vodka per hour per gambler (who is broke)
%Every ~417 times of (i), one hour passes, this is a conservative estimate,
%players could go broke in the middle of an hour and this effects vodka
%served and the estimate of vodka consumption

% Create initial variables
gamblersP3 = ones(1,10000)*10000; %Each person starts with 10,000
gamblingTimesP3 = 1:1:10000; %Assume each person goes for 10,000 bets/day
gamblerTimeSumP3 = [];
vodkaCounter = 0; %For counting up vodka consumption times
vodkaVolume = 0; %ounces; For dispensing complimentary drinks

%This for loop is for one day
for i = 1:1:10000 %for each time possible spent gambling
    vodkaCounter = vodkaCounter +1; %Related time spent to timer that detates vodka consumption
    for j = 1:1:10000 %for each gambler in each time slot, create a new 50/50 probability
        Outcome(j) = round(rand - AdvantageFactor); %create the 50/50 chance
        if Outcome(j) == 0 && gamblersP3(j) >99 %gambler loses and can pay out
            gamblersP3(j) = gamblersP3(j) - 100; %If lose, detract
        end
        if Outcome(j) == 1 && gamblersP3(j) > 99 %gambler wins and hasn't gone broke
            gamblersP3(j) = gamblersP3(j) + 100; %If win, increase
        end
    end
    %Calculate over-time total gambler
    gamblerTimeSumP3.append(i) = sum(gamblersP3); %run sum of all gamblers money per time unit
    
    
    if vodkaCounter == 417
        for k = 1:1:10000;
            if gamblersP3(k) < 100 %Go through list of gamblers
                vodkaVolume = vodkaVolume + 2; %ounces
            end
        end
        vodkaCounter = 0; %reset the counter that corresponds to time
    end
end

%Total ounces consumed
print('vodkaVolume')


toc