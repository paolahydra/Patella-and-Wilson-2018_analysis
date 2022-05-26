% goal here is to refine the tdtomato masks calculated from tdt+ flies.
% I have run the first block(s) in WEDAMMC_AllMapstdTomatoAlignment.m up to
% 
% tdFlies = [119 122 126 121 138];
% for f = 1:length(tdFlies)
%     maskst = allStacks(allStacks_flyNums==tdFlies(f)).masked_alignedstackR;
%     zof = unique(dataTable.zOffsetAlStack(dataTable.fly == tdFlies(f)));
%     maxiTDTMaskStack(:,:, 1+zof : size(maskst,3)+zof, f) = maskst;
% end



% This code instead comes originally from
% final_TdTomatoBoundary_stack2fly150_Registration.m where the masks were made

%% added for running
load('/Users/galileo/Dropbox (HMS)/Data/TDTstacks/newReferenceStack_fly119_cropped.mat'); %reference reference_MP


%%
strelSizeDilate = 8;
    se = strel('disk',strelSizeDilate,0);
strelSizeErode = 2;
    seErode = strel('disk',strelSizeErode,0);
strelSizeClose = 11;
arbitraryThreshold = 0.6; %between 0 and 1
    sizeMAP = size(reference_MP);
    
    
for f = 1:length(tdFlies) % length(allStacks)
    disp(f)
    flyNum_i = allStacks_flyNums==tdFlies(f);
% ... deleted code

    %% tdtomato mask
    masked_alignedstackR = zeros([sizeMAP, size(allStacks(flyNum_i).stackG,3)]);
    boundaries = cell([1, size(allStacks(flyNum_i).stackG,3)]);
    for z = allStacks(flyNum_i).MPrange      
figure;
subplot(1,2,1)
        tdTomPositivePxls = false(sizeMAP);
        tdTomPositivePxls(allStacks(flyNum_i).alignedstackR(:,:,z)>=arbitraryThreshold) = 1;
        tdTomPositivePxls = logical(tdTomPositivePxls);
imshow(tdTomPositivePxls);
        
        tdTomPositivePxls = imdilate(tdTomPositivePxls, se);
%         figure;  title(sprintf('dilate %d',strelSizeDilate))

        tdTomPositivePxls = imclose(tdTomPositivePxls, strel('disk',strelSizeClose,0));
%         figure; imshow(tdTomPositivePxls); title(sprintf('close %d',strelSizeClose))
        
        %new:
        tdTomPositivePxls = imfill(tdTomPositivePxls, 'holes');
%         figure; imshow(tdTomPositivePxls); title('imfill holes')
        
        tdTomPositivePxls = imerode(tdTomPositivePxls, seErode);
subplot(1,2,2)
imshow(tdTomPositivePxls); title(sprintf('erode %d',strelSizeErode))
               
%         figure; imagesc(tdTomPositivePxls); axis image; colormap(gray);
        masked_alignedstackR(:,:,z) = tdTomPositivePxls;
       %         figure; imagesc(masked_alignedstackR(:,:,z)); axis image; colormap(gray);
        [B,L] = bwboundaries(tdTomPositivePxls,'noholes');
        boundaries{z} = B;
%         boundary = B{1};
%         plot(boundary(:,2), boundary(:,1), 'LineWidth', 1)
    end
    size(allStacks(flyNum_i).masked_alignedstackR)
    size(masked_alignedstackR)
    allStacks(flyNum_i).masked_alignedstackR = masked_alignedstackR;
    allStacks(flyNum_i).aligned_Rboundaries = boundaries;
% allStacks(flyNum_i).masked_alignedstackR = [];
% allStacks(flyNum_i).aligned_Rboundaries = [];
    
end        
%%
save(saveName, 'allStacks', '-v7.3')