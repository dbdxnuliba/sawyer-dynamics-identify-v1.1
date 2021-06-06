%%
% filepath = fullfile('no1_dyn_identify_v1.bag');
% bag = rosbag(filepath);
% message = select(bag,'MessageType','std_msgs/Float64MultiArray');
% data = readMessages(message);

count = 1001;

q = zeros(count,7);
qDot = zeros(count,7);
tau = zeros(count,7);
t = zeros(count,1);

for i=1:count
    for j=1:7
        q(i,j) = data{i}.Data(j);
        qDot(i,j) = data{i}.Data(7+j);
        tau(i,j) = data{i}.Data(7*2+j);
    end
    t(i) = data{i}.Data(7*3+1);
end

% %
close all;

figure;
plot(t,q(:,1),t, q(:,2),t, q(:,3),t, q(:,4),t, q(:,5),t, q(:,6),t, q(:,7),'linewidth',1);
xlabel('time (sec)');
ylabel('position (rad)');
legend('$q_{1}$','$q_{2}$','$q_{3}$','$q_{4}$','$q_{5}$','$q_{6}$','$q_{7}$','Interpreter','Latex');

figure;
plot(t,qDot(:,1),t, qDot(:,2),t, qDot(:,3),t, qDot(:,4),t, qDot(:,5),t, qDot(:,6),t, qDot(:,7),'linewidth',1);
xlabel('time (sec)');
ylabel('velocity (rad)');
legend('$\dot{q}_{1}$','$\dot{q}_{2}$','$\dot{q}_{3}$','$\dot{q}_{4}$','$\dot{q}_{5}$','$\dot{q}_{6}$','$\dot{q}_{7}$','Interpreter','Latex');

figure;
plot(t,tau(:,1),t, tau(:,2),t, tau(:,3),t, tau(:,4),t, tau(:,5),t, tau(:,6),t, tau(:,7),'linewidth',1);
xlabel('time (sec)');
ylabel('velocity (rad)');
legend('$\tau_{1}$','$\tau_{2}$','$\tau_{3}$','$\tau_{4}$','$\tau_{5}$','$\tau_{6}$','$\tau_{7}$','Interpreter','Latex');











