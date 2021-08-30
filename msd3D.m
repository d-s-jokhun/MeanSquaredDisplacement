%%%%% MSD in 3-Dimensions %%%%% 
%%% position_1.csv contains all the dots from cell 1 %%%
%%% Results for all particles from all cells are in MSD_tab %%%
%%% 1st column gives the time lag and last column gives the average for all particles %%%



clc;
clear all;

for c = 1:2   %%% c will be 1 to the number of cells (csv files) which have to be processed  %%%  
 
    MSD_n = [];
    MSD_c = [];
    
    cell= strcat ('position_',num2str(c),'.csv')    %%% cell will be the name of the file to be processed (corresponds to cell c)
 position = importdata(cell, ',' , 3); % Imports the data from the csv file 'position.csv'
                                            % The number 3 is telling the program to load only as from the 4th line of the file (1st 3 lines are headers)

 N = size (position.data) ;
 N = ( N(1,2) - 1 ) /3    % n = number of particles being tracked
 image_rate = position.data(2,1) - position.data(1,1)   % 2nd time value - 1st time value will give the time duration bewteen successive time points
 MSD_n (1,1:N) = 1:N ;    % column 1 to column N of 1st row will indicate the particle number (Header)
 
 
 
for n = 1:N

X = position.data(:, 3*(n-1)+2 ); % 3*(n-1)+2 is the column corresponding to x cordinates of particle n
Y = position.data(:, 3*(n-1)+3 ); % 3*(n-1)+3 is the column corresponding to y cordinate 
Z = position.data(:, 3*(n-1)+4 ); % 3*(n-1)+4 is the column corresponding to z cordinate 


particle_n = vertcat([X Y Z]); % Puts the x,y,z coordinates of particle n in table part_n
                            % row 1 is 1st time point, row 2 is 2nd time point, ...

                            
%%% Finding the rows (empty) which need to be deleted from particle_n %%%
q=isfinite(particle_n);
q=sum(sum(q));
q = q/3 ; % total number of data points for particle n

particle_n(q+1:end,:) = [];   % deletes the rows which don't have any data (NaN)

dt_max = floor (q/3);   % dt_max is the largest time lag unit (in terms of data point) (number of frames) which will give accurate msd for particle n
MSD_tab (2:dt_max+1,1) = (1:dt_max)*image_rate ;    % Row 2 to row dt_max+1 of 1st column will indicate the time lag (in seconds) being considered (Header)

MSD_c (2:dt_max+1,1) = (1:dt_max)*image_rate ; % putting time lag in 1st column of MSD_c

    for dt = 1:dt_max

        dxyz = particle_n(1+dt:end,1:3) - particle_n(1:end-dt,1:3); % dxyz will be a matrix containing all the z,y,z translations of particle n within for time lag=dt
        d = dxyz.^2 ;
        d = d(1:end,1) + d(1:end,2) + d(1:end,3);
        d = sqrt(d);  % d will be a matrix containing all the displacements of particle n for time lag = dt
    
        MSD_n(dt+1,n) = mean(d.^2) ;  % MSD of particle n will be stored in the n+1th column
                                  % MSD for time lag dt will be stored in the dt+1th row

    end




end



w=size(MSD_n);
t=w(1,1);
y=w(1,2);
   lines_to_be_filled = size(MSD_tab) - size(MSD_n);  % the number of extra lines in MSD_tab
                                                      % these rows have to be filled with 0s before we can use horzcat
   lines_to_be_filled = lines_to_be_filled (1,1);
   
   
   MSD_c = horzcat (MSD_c , MSD_n) ;
   filename = strcat ('MSD_',num2str(c),'.csv') ;
   
   
   i = size(MSD_c) ;
MSD_c ( 2:i(1,1) , (i(1,2))+1 ) = sum (MSD_c(2:end,2:end) ,2) ;    % Does a horizontal addition such that the last column of MSD_tab is the sum of all the previous columns except the 1st (which is time lag)

for u = 2:i(1,1) ; 
    MSD_c( u , (i(1,2))+1 ) = ( MSD_c(u,(i(1,2))+1) ) / nnz(MSD_c(u,2:i(1,2))) ; % Finds the average by dividing the horizontal sums by the number of non-zero MSD values in that row
end
   
   
   csvwrite(filename , MSD_c) ;
   
   
   
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
    
    
csvwrite('MSD_all.csv' , MSD_tab)

