P=Data001(:,2);
X=Data001(:,1);
CP=zeros(size(X,1),1);
for i=1:size(X,1)
    for j=1:(i-1)
        CP(i,:)=CP(i,:)+P(j,:);
    end;
end;
CP_trace=[X CP];
plot(CP_trace(:,1),CP_trace(:,2));
