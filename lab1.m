clc;
close all;
clear all;
%1.a
load("SysIdenData_StudentVersion.mat")
t = LogData.time;
y_act = LogData.signals(1).values(:,2);
y_actm = LogData.signals(1).values(:,1);
u_act = LogData.signals(2).values;
Ts = t(2,1)-t(1,1);

%1.b
figure;
subplot(2,1,1);
plot(t,y_act,'b');
hold on;
plot(t,y_actm,'r');
hold off;
title('Actual Output Signal')
xlabel('Time (sec)');
ylabel('Water Level (V)');
ylim([1 4])
legend("Noise-Reduced Output","Measured Output");
subplot(2,1,2);
plot(t,u_act,'b');
title('Actual Input Signal')
xlabel('Time (sec)');
ylabel('Pump Voltage (V)');
ylim([1 3])
legend("Actual Input");

%1.c
i = 1;
y_offset = 0;
while(u_act(i+1) == u_act(i))
    i = i+1;
end
y_offset = mean(y_act(1:i))
%y_offset = y_act(1,1);
u_offset = u_act(1,1);
y = y_act-y_offset;
u = u_act-u_offset;

figure;
subplot(2,1,1);
plot(t,y,'r');
title('Actual Offset-Free Output Signal')
xlabel('Time (sec)');
ylabel('Water Level (V)');
ylim([-2,1])
legend("Actual Output");
subplot(2,1,2);
plot(t,u,'b');
title('Actual Offset-Free Input Signal')
xlabel('Time (sec)');
ylabel('Pump Voltage (V)');
ylim([-0.5 0.5])
legend("Actual Input");

%2.a
Kstart = 3;
Ncount = ceil(size(y, 1)/2);
fan = zeros(Ncount,4);
row_index = 0;
bigY =zeros(Ncount,1);
for N = Kstart:Kstart + Ncount
    row_index = row_index + 1;
    yN = [y(N-1,1),y(N-2,1),u(N-1,1),u(N-2,1)];
    fan(row_index,:) = yN;
    bigY(row_index,:) = y(N);
end
N = Ncount;
%2.b
Theta_j = inv(transpose(fan) * fan) * transpose(fan) * bigY;
a1 = Theta_j(1,1);
a2 = Theta_j(2,1);
b1 = Theta_j(3,1);
b2 = Theta_j(4,1);
G = [0,1;a2,a1];
H = [0;1];
C = [b2,b1];
D = 0;

%2.c
Gz = tf([b1 b2],[1 -a1 -a2] , Ts)
sys = ss(G,H,C,D,Ts)

%3.a
SrcDataToSimulink = [t, u];
A = [1, -a1, -a2];
B = [b1, b2];
%y_test = filter(B, A, u);
y_test = lsim(sys, u, t);
disp('MSE for the whole scequence');
MSE = mean((y_test-y).^2)
disp('MSE for the second half scequence');
MSE3 = mean((y_test(N+1:end)-y(N+1:end)).^2)

%3.b 
figure;
y_sim = filter(B, A, u(N+1:end,:));
y_act = y(N+1:end,:);
t_tr = t(1:N-1,:);
subplot(2,1,1)
plot(t_tr, y_sim,'--', t_tr, y_act);
grid on;
legend('Simulated Output','Actual Output');
xlabel('Time(sec)');
ylabel('Water Level(V)');
title('Offset5-Free Model Verification(2^{nd} Half)');
text(10, 0.7, strcat('MSE = ', num2str(MSE3)));

subplot(2,1,2)
plot(t, y_test, '--', t, y);
grid on;
legend('Simulated Output','Actual Output');
xlabel('Time(sec)');
ylabel('Water Level(V)');
title('Offset-Free Model Verification(Entire)');
text(10, 0.7, strcat('MSE = ', num2str(MSE)));

%3.c
% Post-lab Exercise 1
% y(k) = -a1y(k-1)+b1u(k-1)
%starting from t = 2, N = 452

start = 10;
Matrix2(:, 1) = y(start-1:N-1, :);
Matrix2(:, 2) = u(start-1:N-1, :);
est_val2 = inv(Matrix2'*Matrix2)*Matrix2'*y(start:N,:);
sys3 = tf([est_val2(2)], [1, -est_val2(1)])
G = [0, 1; 0, est_val2(1)];
H = [0;1];
C = [0, est_val2(2)];
D = 0;
sys4 = ss(G, H, C, D, Ts)
% sys4 = ss(-est_val(1), est_val(2))
A = [1, -est_val2(1)];
B = [est_val2(2)];
% y_test2 = filter(B, A, u);
y_test2 = lsim(sys4, u, t);
MSE2 = mean((y_test2-y).^2)
figure(4)
plot( t, y, 'r', t, y_test2, 'g.', t, y_test,'b--');
grid on;
legend('Actual Output', '1^{st}Order Model Response', '2^{nd}Order Model Response');
xlabel('Time(sec)');
ylabel('Water Level(V)');
title('Comparison of Different Offset-Free Models');
text(10, 0.65, strcat('MSE1 = ', num2str(MSE2)));
text(10, 0.35, strcat('MSE2 = ', num2str(MSE)));
% Post-lab Exercise 2
