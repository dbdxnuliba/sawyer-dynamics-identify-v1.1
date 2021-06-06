%%
function [q,dq,ddq] = generate_fourier_trajectory(x,wf,t)

dof=evalin('base','dof');
N=5;

q=zeros(dof,1);
dq=zeros(dof,1);
ddq=zeros(dof,1);

for i=1:dof
    for j=1:N
        a=x((i-1)*(2*N)+j);
        b=x((i-1)*(2*N)+N+j);
        q(i)=q(i)+a/(wf*j)*sin(wf*j*t)-b/(wf*j)*cos(wf*j*t);
        dq(i)=dq(i)+a*cos(wf*j*t)+b*sin(wf*j*t);
        ddq(i)=ddq(i)-wf*a*j*sin(wf*j*t)+wf*b*j*cos(wf*j*t);
    end
end

end
