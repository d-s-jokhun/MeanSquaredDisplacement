%%%%% MSD in each of the 2 dimensions %%%%% 

clc;
clear all;

sample='4FPS_DN-KASH';

duration=2 %in min

MSD_tab=[];
MSD_by_Tau=[];

filenames = dir (['nuc_centroid_',sample,'*.mat.xls']);

for c = 1:size(filenames,1)   %%% c will be being processed  %%%  
 
    MSD_n = [];
    MSD_c = [];
    position = [];
    
    [~,~,telo_pos_xls]=xlsread(filenames(c).name);
    
    num_of_tracks = size(telo_pos_xls,2)
    FPS=(size(telo_pos_xls,1)-1)/(duration*60)
    time_lapse=1000/FPS
    
    position(1:size(telo_pos_xls,1)-1,1)=0:size(telo_pos_xls,1)-2; %writes frame number in 1st column
    
    for count_track=1:num_of_tracks/2
    for count_frame=1:size(telo_pos_xls,1)-1
        position(count_frame,count_track*2)=(telo_pos_xls{count_frame+1,(count_track*2)-1}); %fills the x coordinates 
        position(count_frame,(count_track*2)+1)=(telo_pos_xls{count_frame+1,count_track*2}); %fills the y coordinates 
    end
    end
    
%     for count_track=1:num_of_tracks
%     for count_frame=1:size(telo_pos_xls,1)-2
%         position(count_frame,count_track*2)=str2num(telo_pos_xls{count_frame+2,(count_track*3)-1}); %fills the x coordinates 
%         position(count_frame,(count_track*2)+1)=str2num(telo_pos_xls{count_frame+2,count_track*3}); %fills the y coordinates 
%     end
%     end


 MSD_n (1,1:num_of_tracks) = 1:num_of_tracks ;    % column 1 to column N of 1st row will indicate the particle number (Header)
 MSD_n (2,1:num_of_tracks) = 0 ;    % column 1 to column N of 2nd row will be 0 (MSD at time lag zero is zero)
 
 
for n = 1:num_of_tracks

particle_n = position(:,n+1); 

dt_max = floor (size(particle_n,1)/3);   % dt_max is the largest time lag unit (in terms of data point) (number of frames) which will give accurate msd for particle n
MSD_tab (3:dt_max+2,1) = ((1:dt_max)*time_lapse)/1000 ;    % Row 2 to row dt_max+1 of 1st column will indicate the time lag (in mseconds) being considered (Header)

MSD_c (3:dt_max+2,1) = ((1:dt_max)*time_lapse)/1000 ; % putting time lag in 1st column of MSD_c

    for dt = 1:dt_max

        d = particle_n(1+dt:end,:) - particle_n(1:end-dt,:); % dxyz will be a matrix containing all the z,y translations of particle n within for time lag=dt
    
        MSD_n(dt+2,n) = mean(d.^2) ;  % MSD of particle n will be stored in the nth column
                                  % MSD for time lag dt will be stored in the dt+2th row

    end




end



w=size(MSD_n);
t=w(1,1);
y=w(1,2);
   lines_to_be_filled = size(MSD_tab) - size(MSD_n);  % the number of extra lines in MSD_tab
                                                      % these rows have to be filled with 0s before we can use horzcat
   lines_to_be_filled = lines_to_be_filled (1,1);
   
   
   MSD_c = horzcat (MSD_c , MSD_n) ;
   filename = strcat ('MSD_',sample,'_',num2str(c),'.xls') ;
   

   
   
   if lines_to_be_filled > 1
   MSD_n(t+1:t+lines_to_be_filled,:)=0 ;
   end
   
   MSD_tab = horzcat (MSD_tab,MSD_n);

   
   
   
    

end

e = size(MSD_tab) ;
MSD_tab ( 2:e(1,1) , (e(1,2))+1 ) = sum (MSD_tab(2:end,2:end) ,2) ;    % Does a horizontal addition such that the last column of MSD_tab is the sum of all the previous columns except the 1st (which is time lag)

for r = 2:e(1,1) ; 
    MSD_tab( r , (e(1,2))+1 ) = ( MSD_tab(r,(e(1,2))+1) ) / nnz(MSD_tab(r,2:e(1,2))) ; % Finds the average by dividing the horizontal sums by the number of non-zero MSD values in that row
end
    
    
% xlswrite(['MSD_',sample,'_all.xls'] , MSD_tab)

MSD_by_Tau=MSD_tab;
for count=2:size(MSD_by_Tau,2)
    MSD_by_Tau(3:end,count)=MSD_by_Tau(3:end,count)./MSD_by_Tau(3:end,1);
end

MSD_by_Tau_X(:,1)=MSD_by_Tau(:,1);
MSD_by_Tau_Y(:,1)=MSD_by_Tau(:,1);
for count = 1:((size(MSD_by_Tau,2))-1)/2
    MSD_by_Tau_X(:,count+1)=MSD_by_Tau(:,count*2);
    MSD_by_Tau_Y(:,count+1)=MSD_by_Tau(:,count*2+1);
end

MSD_by_Tau_X(1,:)=0:size(MSD_by_Tau_X,2)-1;
MSD_by_Tau_Y(1,:)=0:size(MSD_by_Tau_Y,2)-1;

MSD_by_Tau_X(2:end,end+1)=mean(MSD_by_Tau_X(2:end,2:end),2);
MSD_by_Tau_Y(2:end,end+1)=mean(MSD_by_Tau_Y(2:end,2:end),2);



xlswrite(['MSD by Tau_X_',sample,'.xls'] , MSD_by_Tau_X) ;
xlswrite(['MSD by Tau_Y_',sample,'.xls'] , MSD_by_Tau_Y) ;