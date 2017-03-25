%this function is going to help us with finding the points
%we have to connect our cluster1 corners to.
function [cpoints1, cpoints2] = find_where_to_go(pcluster1, pcluster2)
%first we have the pcluster1, now:
%initialization of the points:
cpoints1 = zeros(5,2);
cpoints2 = zeros(5,2);

for i = 1:size(pcluster1,1)
    %for each point in pcluster1 (cluster 1), we are going to find 
    %the two closest points in cluster2, and manage to connect them.
    differences = pcluster2-repmat(pcluster1(i,:),[5,1]);
    distances = zeros(size(differences,1),1);
    for j = 1:size(distances,1)
        distances(j,1) = norm(differences(j,:));
    end
    [a,b] = min(distances);
    cpoints1(i,:) = pcluster2(b,:);
    distances(b,1)=Inf;
    [a,b] = min(distances);
    cpoints2(i,:) = pcluster2(b,:);
    
end

end