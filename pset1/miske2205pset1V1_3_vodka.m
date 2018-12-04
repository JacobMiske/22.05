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
vodkaTimes =[]

%This for loop is for one day
for i = 1:1:10000 %for each time possible spent gambling
    vodkaCounter = vodkaCounter +1; %Related time spent to timer that detates vodka consumption
    for j = 1:1:10000 %for each gambler in each time slot, create a new 50/50 probability
        Outcome(j) = round(rand - AdvantageFactor/2); %create the 50/50 chance
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
    vodkaTimes.append(i) = vodkaVolume
end
%resolves type issue
gamblerTimeSumArrayP3 = struct2array(gamblerTimeSumP3);
%Plot similar to Figure 1, for house advantage of 0.001
figure(4)
plot(gamblingTimesP3, gamblerTimeSumArrayP3, 'b');
savefig('Gambler Earnings over time with Vodka Serving.fig')
title('Total Gambler Winnings over Time with Vodka Serving')
xlabel('Times to Gamble (1 = 8.64seconds)')
ylabel('Dollars')
saveas(gcf,'Total Gambler Winnings over Time with Vodka Serving.pdf')

%Total ounces consumed
vodkaVolume
vodkaTimes = struct2array(vodkaTimes)
figure(8)
plot(gamblingTimesP3, vodkaTimes, 'g');
savefig('Vodka Ounces Served over time.fig')
title('Vodka Ounces Served over time')
xlabel('Times to Gamble (1 = 8.64seconds)')
ylabel('Ounces')
saveas(gcf,'Vodka Ounces Served over time.pdf')
