function [flyClusterCountDivZs] = barcluster_panSpec_maxZout(T, datalist, flyfoldersUnique, flyClusters, sortFliesGenoType)

snK = length(unique(T));
NaF = length(flyfoldersUnique);
flyClusterCountDivZs = zeros(NaF, snK);
for f = 1 : NaF 
    flyfolder = flyfoldersUnique{f};
    disp(flyfolder)
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);
    flyClusters_fly = flyClusters(:,:,zs_fly);
    flyClusters_fly = permute(flyClusters_fly, [1,3,2]);
%     disp( max (permute(sum(flyClusters_fly, 1), [2,3,1]), [], 1) )
    flyClusterCountDivZs(f,:) = max (permute(sum(flyClusters_fly, 1), [2,3,1]), [], 1);
end

if sortFliesGenoType
    sortFMatrix = [1     2     3     4     5     6    11     8    17     9    12    13    15     7    18    19    16    22    23    24    14   10    20    21];
    [~,idx_sortingFlies] = sort(sortFMatrix);
    flyClusterCountDivZs = flyClusterCountDivZs(idx_sortingFlies,:);
end

