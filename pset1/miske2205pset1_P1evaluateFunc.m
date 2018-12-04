%% Plot relation for P1
%Function described in writeup
HotelDayHoldings = @(A) -7000000000*A + 40000
HotelDayHoldingsUpper = @(A) -7000000000*A + 40000 + 22612
HotelDayHoldingsLower = @(A)-7000000000*A + 40000 - 22612

figure(9)
fplot(HotelDayHoldings, 'k'); hold on
fplot(HotelDayHoldingsUpper, 'g-'); hold on
fplot(HotelDayHoldingsLower, 'b-');

xlim([0 0.02])
ylim([-0.2e9 0.1e9])
title('Start of Day Hotel Holdings verse Advantage factor')
xlabel('Advantage Factor')
ylabel('House Holdings at beginning of day, without t-value')
saveas(gcf,'P1 relation.pdf')