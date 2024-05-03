function [inhale_to_exhale, exhale_to_inhale, point30Percent, point70Percent]=exhale_inhale_area_ratio(respiration,Tf)

%
extrMaxIndex = find(diff(sign(diff(respiration)))==-2)+1;
extrMinIndex = find(diff(sign(diff(respiration)))==+2)+1;

%Find the boundary points of the 30%-70% area
extrIndex = sort([extrMaxIndex,extrMinIndex]);
[~,n]=size(extrIndex);
for i=1:n-1
    Point30Percent(i) = 0.7*respiration(extrIndex(i+1))+0.3*respiration(extrIndex(i));
    Point70Percent(i) = 0.3*respiration(extrIndex(i+1))+0.7*respiration(extrIndex(i));
end

Point30PercentOdd = Point30Percent(1:2:end);
Point30PercentEven = Point30Percent(2:2:end);
Point70PercentOdd = Point70Percent(1:2:end);
Point70PercentEven = Point70Percent(2:2:end);

point30Percent=[];
point70Percent=[];

point30Percent(1:2:2*length(Point30PercentOdd))=Point30PercentOdd;
point30Percent(2:2:2*length(Point70PercentEven))=Point70PercentEven;
point70Percent(1:2:2*length(Point70PercentOdd))=Point70PercentOdd;
point70Percent(2:2:2*length(Point30PercentEven))=Point30PercentEven;

point30Percent(2,:)=point30Percent;
point70Percent(2,:)=point70Percent;

for i=1:n-1
    sum=abs(respiration(extrIndex(i))-point30Percent(2,i));
    for j=extrIndex(i):extrIndex(i+1)
        if abs(respiration(j)-point30Percent(2,i))<sum
           sum=abs(respiration(j)-point30Percent(2,i));
           point30Percent(1,i)=j;
        end
    end
end

for i=1:n-1
    sum=abs(respiration(extrIndex(i))-point70Percent(2,i));
    for j=extrIndex(i):extrIndex(i+1)
        if abs(respiration(j)-point70Percent(2,i))<sum
           sum=abs(respiration(j)-point70Percent(2,i));
           point70Percent(1,i)=j;
        end
    end
end

point30Percent(2,:)=respiration(point30Percent(1,:));
point70Percent(2,:)=respiration(point70Percent(1,:));
point30Percent(1,:)=point30Percent(1,:)*Tf;
point70Percent(1,:)=point70Percent(1,:)*Tf;

%find quadrilateral area
for i=1:n-2
    AB(i)=sqrt((point30Percent(1,i+1)-point30Percent(1,i))^2+(point30Percent(2,i+1)-point30Percent(2,i))^2);
    BC(i)=sqrt((point70Percent(1,i+1)-point30Percent(1,i+1))^2+(point70Percent(2,i+1)-point30Percent(2,i+1))^2);
    CD(i)=sqrt((point70Percent(1,i+1)-point70Percent(1,i))^2+(point70Percent(2,i+1)-point70Percent(2,i))^2);
    AD(i)=sqrt((point70Percent(1,i)-point30Percent(1,i))^2+(point70Percent(2,i)-point30Percent(2,i))^2);
    BD(i)=sqrt((point30Percent(1,i+1)-point70Percent(1,i))^2+(point30Percent(2,i+1)-point70Percent(2,i))^2);
    s(i)=1/2*AB(i)*AD(i)*sqrt(1-((AB(i)^2+AD(i)^2-BD(i)^2)/(2*AB(i)*AD(i)))^2)+1/2*BC(i)*CD(i)*sqrt(1-((BC(i)^2+CD(i)^2-BD(i)^2)/(2*BC(i)*CD(i)))^2);
end

%Classification of expiratory and inspiratory regions
if extrMaxIndex(1)>extrMinIndex(1)
    s_inhale_area=s(1:2:end);
    s_exhale_area=s(2:2:end);
else
    s_exhale_area=s(1:2:end);
    s_inhale_area=s(2:2:end);
end

%breath area ratio
[~,m]=size(s_exhale_area);
[~,n]=size(s_inhale_area);
sum_exhale_area=0;
sum_inhale_area=0;
for i=1:m
    sum_exhale_area=sum_exhale_area+s_exhale_area(i);
end
for j=1:n
    sum_inhale_area=sum_inhale_area+s_inhale_area(j);
end
inhale_to_exhale=sum_inhale_area/sum_exhale_area*m/n;
exhale_to_inhale=sum_exhale_area/sum_inhale_area*n/m;

end