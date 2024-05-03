function [ index, XT ] = cfar_os( xc, N, k, pro_N, PAD)
%   It is assumed that the echo follows a Gaussian distribution
%   alpha There are some problems with assignment, a more complex higher-order function

%% calculate alpha
% syms alpha PFA;
% PFA(alpha)=gamma(N-1).*gamma(N-k+alpha-1)./gamma(N-k-1)./gamma(N+alpha-1);
% [alpha,~,~]=solve(PFA(alpha)==PAD,'ReturnConditions', true) ;

alpha=N.*(PAD.^(-1./N)-1);

index=1+N/2+pro_N/2:length(xc)-N/2-pro_N/2;
XT=zeros(1,length(index));

for i=index
    cell_left=xc(1,i-N/2-pro_N/2:i-pro_N/2-1);
    cell_right=xc(1,i+pro_N/2+1:i+N/2+pro_N/2);
    cell_all=cat(2,cell_left,cell_right);
    cell_sort=sort(cell_all);

    Z=cell_sort(1,k);

    XT(1,i-N/2-pro_N/2)=Z.*alpha;
end

end