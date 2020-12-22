function [Q,R]=qrsim(A)

[m,n]=size(A);
Q=zeros(m);
R=zeros(m);

alpha=zeros(m,n);
beta=zeros(m,n);
K=eye(m,n);
for i=1:m
    alpha(:,i)=A(:,i);
end
for i=1:m
    beta(:,i)=alpha(:,i);
    for j=1:i-1
        K(j,i)=dot(alpha(:,i),beta(:,j))/dot(beta(:,j),beta(:,j)) ;
        beta(:,i)=beta(:,i)-K(j,i)*beta(:,j);
    end
end
for i=1:n
    Q(:,i)=beta(:,i)/norm(beta(:,i));
    Beta(i)=norm(beta(:,i));
end
R=diag(Beta)*K; 
