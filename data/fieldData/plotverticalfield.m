%% Import data from text file.
% Script for importing data from the following text file:
%
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2021/05/11 15:06:31

NBR = 1019;
START_ET = 1;
END_ET = NBR;
loops = END_ET - START_ET;
Mov(loops) = struct('cdata',[],'colormap',[]);
index = 1;

for j = 1:NBR
    fprintf('Field layer:%d\n', j);
    %% Initialize variables
    
%############# USER DEFINED ###############

    % Absolute path to samplelayer files
    %fileBase = 'C:\github\bremmer_master\vsproject\blood3D\data\fieldData\migratedfield';
    %fileBase2 = 'C:\github\bremmer_master\vsproject\blood3D\data\fieldData\imaginaryfield';
    fileBase = 'C:\github\bremmer_master\data\fieldData\migratedfield';
    fileBase2 = 'C:\github\bremmer_master\data\fieldData\imaginaryfield';
    
%##########################################
    %% Initialize variables.
    
    fileNumber = string(j-1);
    fileType = '.txt';
    filename = strcat(fileBase, fileNumber, fileType);
    filename2 = strcat(fileBase2, fileNumber, fileType);
    delimiter = '\t';
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

    %% Open the text file.
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

    fileID2 = fopen(filename2,'r');
    dataArray2 = textscan(fileID2, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    
    %% Close the text file.
    fclose(fileID);
    fclose(fileID2);

    %% Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.

    %% Create output variable
    migratedfield = [dataArray{1:end-1}];
    imaginaryfield = [dataArray2{1:end-1}];
   
    % Cut out coloumn
    M(j,:)=migratedfield(:,512);
    K(j,:)=migratedfield(512,:);
    O(j,:)=imaginaryfield(:,512);
    P(j,:)=imaginaryfield(512,:);
    
    % Save field in focal line
    E_focalline(j) = migratedfield(512,j);
    
    if j >= START_ET
        if j <= END_ET
            s=pcolor(migratedfield);
            s.FaceColor='interp';
            set(s, 'EdgeColor', 'none')
            colorbar
            caxis([-1.5 1.5])
            %drawnow
            Mov(index) = getframe;
            index=index+1;
        end
    end

    clearvars migratedfield;
    
end

figure
movie(Mov)

figure
s=pcolor(M);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
colorbar
caxis([0 1.5])
s.FaceColor='interp';

hold on
contour(verticalcut, 'black')

figure
s=pcolor(K);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
colorbar
caxis([0 1.5])
s.FaceColor='interp';


hold on
contour(verticalcut2, 'black')

figure
s=pcolor(O);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
colorbar
caxis([0 1.5])
s.FaceColor='interp';


hold on
contour(verticalcut, 'black')

figure
s=pcolor(P);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
colorbar
caxis([0 1.5])
s.FaceColor='interp';


hold on
contour(verticalcut2, 'black')

% Save movie
myVideo = VideoWriter('myfile.avi');
uncompressedVideo = VideoWriter('myfile.avi', 'Uncompressed AVI');
open(myVideo);
writeVideo(myVideo, Mov);
close(myVideo);