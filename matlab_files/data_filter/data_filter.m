N=5; %滤波器阶数
fs=100; %1000hz 系统采样周期
fc=3; %3hz 截止频率

[m,n]=size(q);

%%
if(0)
q_filter=zeros(m,n);

for i=1:n
    q_filter(:,i)=digital_filter(N, fs, 3, q(:,i));
end

close all;
figure;
plot(t,q(:,1),t,q_filter(:,1));
legend('$q_{1}$', '$q_{1}$ filter','Interpreter','Latex');

figure;
plot(t,q(:,2),t,q_filter(:,2));
legend('$q_{2}$', '$q_{2}$ filter','Interpreter','Latex');

figure;
plot(t,q(:,3),t,q_filter(:,3));
legend('$q_{3}$', '$q_{3}$ filter','Interpreter','Latex');

figure;
plot(t,q(:,4),t,q_filter(:,4));
legend('$q_{4}$', '$q_{4}$ filter','Interpreter','Latex');

figure;
plot(t,q(:,5),t,q_filter(:,5));
legend('$q_{5}$', '$q_{5}$ filter','Interpreter','Latex');

figure;
plot(t,q(:,6),t,q_filter(:,6));
legend('$q_{6}$', '$q_{6}$ filter','Interpreter','Latex');

figure;
plot(t,q(:,7),t,q_filter(:,7));
legend('$q_{7}$', '$q_{7}$ filter','Interpreter','Latex');
end

%%
if(0)
qDot_filter=zeros(m,n);

for i=1:n
    qDot_filter(:,i)=digital_filter(N, fs, 3, qDot(:,i));
end

close all;
figure;
plot(t,qDot(:,1),t,qDot_filter(:,1));
legend('$\dot{q}_{1}$', '$\dot{q}_{1}$ filter','Interpreter','Latex');

figure;
plot(t,qDot(:,2),t,qDot_filter(:,2));
legend('$\dot{q}_{2}$', '$\dot{q}_{2}$ filter','Interpreter','Latex');

figure;
plot(t,qDot(:,3),t,qDot_filter(:,3));
legend('$\dot{q}_{3}$', '$\dot{q}_{3}$ filter','Interpreter','Latex');

figure;
plot(t,qDot(:,4),t,qDot_filter(:,4));
legend('$\dot{q}_{4}$', '$\dot{q}_{4}$ filter','Interpreter','Latex');

figure;
plot(t,qDot(:,5),t,qDot_filter(:,5));
legend('$\dot{q}_{5}$', '$\dot{q}_{5}$ filter','Interpreter','Latex');

figure;
plot(t,qDot(:,6),t,qDot_filter(:,6));
legend('$\dot{q}_{6}$', '$\dot{q}_{6}$ filter','Interpreter','Latex');

figure;
plot(t,qDot(:,7),t,qDot_filter(:,7));
legend('$\dot{q}_{7}$', '$\dot{q}_{7}$ filter','Interpreter','Latex');
end

%%
if(0)
qDDot_diff=zeros(m,n);
qDDot_filter=zeros(m,n);

for i=1:n
    qDDot_diff(:,i)=digital_differentiator(fs, qDot_filter(:,i));
    qDDot_filter(:,i)=digital_filter(N, fs, 3, qDDot_diff(:,i));
end

close all;
figure;
plot(t,qDDot_diff(:,1),t,qDDot_filter(:,1));
legend('$\ddot{q}_{1}$', '$\ddot{q}_{1}$ filter','Interpreter','Latex');

figure;
plot(t,qDDot_diff(:,2),t,qDDot_filter(:,2));
legend('$\ddot{q}_{2}$', '$\ddot{q}_{2}$ filter','Interpreter','Latex');

figure;
plot(t,qDDot_diff(:,3),t,qDDot_filter(:,3));
legend('$\ddot{q}_{3}$', '$\ddot{q}_{3}$ filter','Interpreter','Latex');

figure;
plot(t,qDDot_diff(:,4),t,qDDot_filter(:,4));
legend('$\ddot{q}_{4}$', '$\ddot{q}_{4}$ filter','Interpreter','Latex');

figure;
plot(t,qDDot_diff(:,5),t,qDDot_filter(:,5));
legend('$\ddot{q}_{5}$', '$\ddot{q}_{5}$ filter','Interpreter','Latex');

figure;
plot(t,qDDot_diff(:,6),t,qDDot_filter(:,6));
legend('$\ddot{q}_{6}$', '$\ddot{q}_{6}$ filter','Interpreter','Latex');

figure;
plot(t,qDDot_diff(:,7),t,qDDot_filter(:,7));
legend('$\ddot{q}_{7}$', '$\ddot{q}_{7}$ filter','Interpreter','Latex');
end

%%
if (1)
tau_filter=zeros(m,n);

for i=1:n
    tau_filter(:,i)=digital_filter(N, fs, 2, tau(:,i));
end

close all;
figure;
plot(t,tau(:,1),t,tau_filter(:,1));
legend('$\tau_{1}$', '$\tau_{1}$ filter','Interpreter','Latex');

figure;
plot(t,tau(:,2),t,tau_filter(:,2));
legend('$\tau_{2}$', '$\tau_{2}$ filter','Interpreter','Latex');

figure;
plot(t,tau(:,3),t,tau_filter(:,3));
legend('$\tau_{3}$', '$\tau_{3}$ filter','Interpreter','Latex');

figure;
plot(t,tau(:,4),t,tau_filter(:,4));
legend('$\tau_{4}$', '$\tau_{4}$ filter','Interpreter','Latex');

figure;
plot(t,tau(:,5),t,tau_filter(:,5));
legend('$\tau_{5}$', '$\tau_{5}$ filter','Interpreter','Latex');

figure;
plot(t,tau(:,6),t,tau_filter(:,6));
legend('$\tau_{6}$', '$\tau_{6}$ filter','Interpreter','Latex');

figure;
plot(t,tau(:,7),t,tau_filter(:,7));
legend('$\tau_{7}$', '$\tau_{7}$ filter','Interpreter','Latex');
end
