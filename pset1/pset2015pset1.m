%Jacob Miske
%22.05 2015 PSET 1

%Create variables
gamblers = ones(1,10000)*25000; %Each person starter
gamblingTimes = 1:1:10000; %Assume each person goes for 10,000 bets/day
gamblerTimeSum = [];
%Let's say if "round(rand)" produces 0, the gambler is wrong at that
%instance of the game

vodkaCounter = 0; %For registering a new hour
vodkaVolume = 0; %ounces; For dispensing complimentary drinks
for i = 1:1:10000 %for each time possible spent gambling
    vodkaCounter = vodkaCounter +1;
    for j = 1:1:10000 %for each gambler in each time slot
        Outcome(j) = round(rand); %create a 50/50
        if Outcome(j) == 0 && gamblers(j) >249
            gamblers(j) = gamblers(j) - 250; %If lose, detract
        end
        if Outcome(j) == 1 && gamblers(j) > 249
            gamblers(j) = gamblers(j) + 250; %If win, increase
        end
    end
    %Calculate over-time total gambler 
    gamblerTimeSum.append(i) = sum(gamblers); %run sum of all gamblers money per time unit
    
    %Calculate complimentary vodka
    %Casino is ~24hr, two ounces of vodka per hour per gambler (not broke)
    %Every ~417 times of (i), one hour passes, this is a conservative
    %estimate of vodka consumption
    if vodkaCounter == 417
        for k = 1:1:10000;
            if gamblers(k) > 0 %Go through list of gamblers
                vodkaVolume = vodkaVolume + 2; %ounces
            end
        end
        vodkaCounter = 0; %reset
    end
        
end
%Create an array from the appended structure
gamblerTimeSumArray = struct2array(gamblerTimeSum);
%Create plot of total gambler earnings at each time interval
plot(gamblingTimes, gamblerTimeSumArray)
title('Total Gambler Winnings over Time')
xlabel('Dollars')
ylabel('Times to Gamble')
%Question 1: Probability that some gambler ends day with >$50000
%Check through gamblers array for those who ended with >$50000
greater50kGamblers = 0;
for i = 1:1:10000;
    if gamblers(i) >50000
        greater50kGamblers = greater50kGamblers+1;
    end
end

%Question 2: Two ounces of vodka per gambler, not broke, per hour
print('vodkaVolume')

%Question 3: How much cash should casino start with each day to have a 99%
%chance of making necessary payouts? Higher gameOdds hurt casino!!!
gameOdds = 0.5;
Advantage = 2*gameOdds - 1;
%If gameOdds = 0.5, casino should break even. Assuming 10,000 gamblers each
%making $250 wagers, 10000 times... Assuming nobody goes broke beforehand.
expectedDailyPayout = gameOdds*250*10000*10000 - 12500000000 %Large number is for gameOdds = 0.5
potentialDailyPayoutArray = []
AdvantageArray = 0.1:0.1:0.9 %Different house advantages to examine
for i =1:1:9;
    potentialDailyPayoutArray(i) = ((AdvantageArray(i)+1)/2)*250*10000*10000 - 12500000000;
end
%Plotting potentialDailyPayoutArray to different house advantage factors
plot(AdvantageArray, potentialDailyPayoutArray, 'ro')
title('House Return versus House Advantage Factor')
xlabel('Advantage Factors')
ylabel('Dollars of House Money At Risk (Starting Cash Minimum')

%Question 4: The casino needs +$50mil for loan payments, what is the
%minimum house advantage the casino needs?
expectedDailyPayoutArray = [];
gameOddsArray = 0.1:0.1:0.9;
%plot expected daily revenue versus house advantage factor
for i = 1:1:9
    expectedDailyPayoutArray(i) = 12500000000 - gameOddsArray(i)*250*10000*10000;
end
plot(AdvantageArray,expectedDailyPayoutArray)
title('Expected Revenue versus House Advantage')
hold on
fiftymillion = 1:1:9;
fiftymillion = fiftymillion*50000000
plot(AdvantageArray,fiftymillion)