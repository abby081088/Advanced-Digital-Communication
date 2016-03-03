clear all
close all

epsilon = 0.1:0.05:0.8;

q = linspace(0,1,1000);

figure(1)
plot(q,1-(1-q).^(1/(6-1)), 'k', 'linewidth',4)
hold on 
for k = 1:length(epsilon)
    plot(q,epsilon(k).*q.^2,'linewidth',2)
end
plot(0.5*ones(length(q),1),q,'k', 'linewidth',2)
hold off
legend('g-1(q))',num2str(epsilon(1)),num2str(epsilon(2)),num2str(epsilon(3)),num2str(epsilon(4)),num2str(epsilon(5)),num2str(epsilon(6)),num2str(epsilon(7)),num2str(epsilon(8)),num2str(epsilon(9)),num2str(epsilon(10)),num2str(epsilon(11)),num2str(epsilon(12)),num2str(epsilon(13))),num2str(epsilon(14)),num2str(epsilon(15));        
xlabel('q')