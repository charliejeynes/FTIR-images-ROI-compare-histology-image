%%% Charlie Jeynes - 14 / 11/ 2018
%%% this scripts takes a FTIR image from a cary agilent 630, where it has
%%% already been normalised. It then kemans the image and displays it next
%%% to an H&E section of an adjacent section of the FTIR imaged slice. The
%%% user can then use a freehand tool to draw around parts of the kmeaned
%%% image while the H&E section is displayed next to it. 
%%% the PURPOSE is so the user can look at the H&E section and pick out ROI
%%% on the FTIR image based on the H&E image

clc
close
clear all
%%
filenames = dir('*.mat');  
%% get all the subdirectories
[subfolders] = subdir('/Volumes/Charlie_mac_backup/TAU Project');%subdir('/Users/jcgj201/Documents/MATLAB/FTIR data/Testing_moveResultsTau_data'); % get all the subdirectories
pattern = ["TG" , "WT"]; 
logicalSub = contains(subfolders, pattern); 
subfolders1 = subfolders(logicalSub); 
subfoldersT = subfolders1';

%% load .mat files 
function D = loadData()
D = load('normal_one_EMSC.mat');

end
%%
function [statsTable] = extractData()
 
%% extract data
normData = D.normal_one_EMSC.Data; 
%% squeeze data
normDataSq = squeeze(normData); 
%% plot check
figure, plot(1:801, normDataSq(:, 200)); 
%% transpose for kmeans: rows instances, columns variables
Y = normDataSq'; 
%% Perform k-means
tic
[dataKmean_idx_6,Clustered_spectra] = kmeans(Y,6); 
toc
%% transform kmeans data to image format and show 
x = D.normal_one_EMSC.xy(2:end, 1);
y = D.normal_one_EMSC.xy(2:end, 2);
originalImageSize = zeros(max(y), max(x));
for k = 1 : length(x)
	originalImageSize(y(k), x(k)) = dataKmean_idx_6(k);
end
figure, imagesc(originalImageSize);
cmap = jet(6); 
cmap = flipud(cmap(1:6,:));
% cmap(1,:) = [1,1,1]; set to white, mycolors = [1 0 0; 1 1 0; 0 0 1]; red, green, blue, scale 0,1 
colormap(cmap);
colorbar
axis('on', 'image');
%%
%%%%%% Draw on the kmeans image and get indices%%%%%%%%
%% Put D.data into original image format, so that the ROI getter works
%%%%%% gets inputs for function "ROIcapture" %%%%%%%%
normDataSqT = normDataSq';  %transpose the normalised data so 59000x801 format
x = D.normal_one_EMSC.xy(2:end, 1);
y = D.normal_one_EMSC.xy(2:end, 2);
OriginalData_Image_format = zeros(max(y), max(x), size(normData,1)); % initialise zero array the same size as the original data
for k = 1 : length(x)
	OriginalData_Image_format(y(k), x(k), :) = normDataSqT(k, :);
end
%%
ImgName = D.normal_one_EMSC.FileName; 

%% Call the makeStatsTable on the image
statsTable = makeStatsTable(originalImageSize,OriginalData_Image_format, ImgName); 

end

%% Plot to check that sensible data is coming  - they are all noramlised to 43.9
%%%%%%%%%%%%%this doesn't work for some
%%%%%%%%%%%%%reason%%%%%%%%%%%%%%%%%%%%%%%%%%
% OriginalData_Image_format_Reshape = reshape(OriginalData_Image_format, 65536, 801); 
% figure, plot(1:801, OriginalData_Image_format_Reshape(200, :)); 
% figure, imagesc(OriginalData_Image_format);
% colorbar
%%
function [statsTable] = makeStatsTable(originalImageSize,OriginalData_Image_format, ImgName)
region = []; 
imagename = {}; 
ROIdata = {}; 
statsTable = table(ROIdata, imagename, region); 


% while true
%   drawnow()
%   stop_state = get(YourPushbuttonHandle, 'Value');
%   if stop_state
%     break;
%   end
j = 1;
    while true
    % for k=1:2

        statsTable.ROIdata{j} = ROIcapture(originalImageSize,OriginalData_Image_format); %This functions passes in the image data % name
        statsTable.region(j) = j; 
        statsTable.imagename{j} = ImgName; 
        j = j+1; 
        fig = figure; 
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        descr = {'click me for another ROI, keyboard press to quit'};
        axes(ax1) % sets ax1 to current axes
        text(.025,0.6,descr)
        w = waitforbuttonpress;
        if w == 0
            disp('Button click')
            disp('continue')
        else
            disp('Key press & break')
            break
        end
    end
end 
%%

function [ROIspectra] = ROIcapture(originalImageSize,OriginalData_Image_format)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is David's impoly bit which allows user to draw a polygon around
% image

%     figure('Name', name); 
%     imagesc(sum(data,3));
    figure, 
    subplot(1,2,2);
    imshow('normal_one_H&E.jpg');
    subplot(1,2,1);
    imagesc(originalImageSize);
    cmap = jet(6); 
    cmap = flipud(cmap(1:6,:));
    % cmap(1,:) = [1,1,1]; set to white, mycolors = [1 0 0; 1 1 0; 0 0 1]; red, green, blue, scale 0,1 
    colormap(cmap);
    colorbar
    axis('on', 'image');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % make it full screen
    
%     ax = axes('Position',[0 0 0.8 1]);
            %"FINISHED-NEXT IMAGE" BUTTON
%         finishButton = uicontrol('String', 'NEXT IMAGE', 'Style', ...
%             'pushbutton', 'Units', 'normalized', 'FontSize', 14);
%         finishButton.Position = [0.8 0.6 0.13 0.05];
%         finishButton.Callback = 'set(gca,''Tag'',''finished'');';
%         finishButton.Enable = 'on';

    ROI  = impoly(gca); %imfreehand imrect(gca);
    wait(ROI); 
    ROIBW = ROI.createMask; 
    [yCoordinates, xCoordinates] = find(ROIBW);

%     saveas(gcf,sprintf('fullImage%d.png',i));
%     close

    x_min = min(xCoordinates); x_max = max(xCoordinates); x_range = x_max-x_min+1;
    y_min = min(yCoordinates); y_max = max(yCoordinates); y_range = y_max-y_min+1;
%     dataROI = -1000*ones(y_range, x_range, size(data,3));
    dataROI = -1000*ones(y_range, x_range, size(OriginalData_Image_format,3)); %intialises a onesROI pixelsxpixelsx801
    for k = 1:numel(yCoordinates)
%         dataROI(yCoordinates(k)-y_min+1, xCoordinates(k)-x_min+1, :) = data(yCoordinates(k), xCoordinates(k), :);
          dataROI(yCoordinates(k)-y_min+1, xCoordinates(k)-x_min+1, :) = OriginalData_Image_format(yCoordinates(k), xCoordinates(k), :);
    end
    dataROI2_im = sum(dataROI,3);
    min_val = min(dataROI2_im(dataROI2_im>=-1000));
    for k1 = 1:y_range
        for k2 = 1:x_range
            if dataROI2_im(k1,k2)<=-1000
                dataROI2_im(k1,k2) = min_val;
            end
        end
    end
    
    ROIspectra = zeros(numel(yCoordinates),size(OriginalData_Image_format,3));
    for k = 1:numel(yCoordinates)
        ROIspectra(k, :) = OriginalData_Image_format(yCoordinates(k), xCoordinates(k), :);
    end

%     ROIspectraMasterNUM1 = size(ROIspectra);
%     ROIspectraMasterNUM2 = ROIspectraMasterNUM1(1); 
%     ROIspectraMasterNUM{i} = ROIspectraMasterNUM2; 
%     ROIspectraMaster{i} = ROIspectra; 
%     ROIspectraMaster1 = cell2mat(ROIspectraMaster'); % create matrix of spectra
%     
%     figure('Name', name);

    figure, imagesc(sum(dataROI2_im,3));
    colorbar
%     saveas(gcf,sprintf('magnified%d.png',i));
end


