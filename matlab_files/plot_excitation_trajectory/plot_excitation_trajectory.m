%%
dof=7;
wf=0.1*2*pi;
Tf=10;
h=0.1;
opt_x=[-0.176324, -0.0653999, 0.0700704, 0.388769, -0.217115, -0.159513, 0.104284, 0.407091, 0.249853, -0.453948, -0.186278, 0.182229, -0.148129, 0.0635953, 0.0885826, -0.132611, 0.171417, 0.167529, 0.217469, -0.316537, -0.0596368, -0.154235, 0.202112, 0.118757, -0.106997, -0.160041, 0.0940674, 0.57065, -0.0845296, -0.280385, 0.416422, -0.278377, -0.226542, -0.00275119, 0.0912474, -0.160754, 0.378351, -0.00661289, -0.0352553, -0.0870176, -0.169239, -0.0028689, 0.275517, 0.228386, -0.331796, -0.219385, 0.247307, 0.494129, 0.0147409, -0.363316, -0.463405, 0.290573, 0.0951751, 0.280039, -0.202382, 0.00476106, -0.119392, 0.138698, 0.177614, -0.178505, 0.164718, -0.373482, 0.249517, 0.0619327, -0.102686, -0.076923, -0.0482629, 0.402511, 0.091639, -0.280128];


count=Tf/h;
opt_q=zeros(count+1,dof);
opt_dq=zeros(count+1,dof);
opt_ddq=zeros(count+1,dof);

t=-h;
for k=1:(count+1)
    t=t+h;
    [opt_q(k,:),opt_dq(k,:),opt_ddq(k,:)]=generate_fourier_trajectory(opt_x,wf,t); 
end

time=0:h:Tf;
close all;
figure;
plot(time,opt_q(:,1),time,opt_q(:,2),time,opt_q(:,3),time,opt_q(:,4),time,opt_q(:,5),time,opt_q(:,6),time,opt_q(:,7),'Linewidth',2);
grid on;
legend('$q_{1}$','$q_{2}$','$q_{3}$','$q_{4}$','$q_{5}$','$q_{6}$','$q_{7}$','Interpreter','Latex');

figure;
plot(time,opt_dq(:,1),time,opt_dq(:,2),time,opt_dq(:,3),time,opt_dq(:,4),time,opt_dq(:,5),time,opt_dq(:,6),time,opt_dq(:,7),'Linewidth',2);
grid on;
legend('$\dot{q}_{1}$','$\dot{q}_{2}$','$\dot{q}_{3}$','$\dot{q}_{4}$','$\dot{q}_{5}$','$\dot{q}_{6}$','$\dot{q}_{7}$','Interpreter','Latex');

figure;
plot(time,opt_ddq(:,1),time,opt_ddq(:,2),time,opt_ddq(:,3),time,opt_ddq(:,4),time,opt_ddq(:,5),time,opt_ddq(:,6),time,opt_ddq(:,7),'Linewidth',2);
grid on;
legend('$\ddot{q}_{1}$','$\ddot{q}_{2}$','$\ddot{q}_{3}$','$\ddot{q}_{4}$','$\ddot{q}_{5}$','$\ddot{q}_{6}$','$\ddot{q}_{7}$','Interpreter','Latex');

