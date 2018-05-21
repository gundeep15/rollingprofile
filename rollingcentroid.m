% Assumption(s): small frame rate
% Note: 'A' matrix can be uncommented if there is a need to actually have
% the data about the movement of each of the centroids is required (for rolling, 
% or velocity, or orientation). If only number of centroids with time is
% needed then no need to utilize the storage of the system

function rollingcentroid
video='500_rot.avi';
vid=VideoReader(video);
T=vid.Duration;
%T=14;
c=100; % number of frames
dt=T/c;
numbercells = zeros (1, c);
time = zeros (1, c);
t=0;
vid.CurrentTime = t;
n=1;
     

     frames1=readFrame(vid);
     frames1=rgb2gray(frames1);
     imgThresh1 = frames1 >90;
     imgfilter1 = bwareaopen (imgThresh1, 95);
     se = strel('disk',4,0);
     imgfilter1= imdilate(imgfilter1,se);
     %figure,imshow (imgfilter1)
     s = regionprops(imgfilter1,'centroid');
     centroids1 = cat(1, s.Centroid);
     centroids1= sortrows(centroids1,2); % will sort the y's in ascending order
     numbercells(n) = size(centroids1,1); % number of centroids/rows in centroid xy matrix 
     
     sz=size(centroids1);

     A = zeros (sz(1,1),sz(1,2),c);
        
        for i=1:sz(1,1)
            A (i,1,n) = centroids1(i,1);
            A (i,2,n) = centroids1(i,2);
        end
    
 
vid.CurrentTime = t+dt;
n=n+1;
while hasFrame(vid)
    frames2=readFrame(vid);
    frames2=rgb2gray(frames2);
    imgThresh2 = frames2 >90;
    imgfilter2 = bwareaopen (imgThresh2, 95);
    se = strel('disk',4,0);
    imgfilter2= imdilate(imgfilter2,se);
    s = regionprops(imgfilter2,'centroid');
    %figure, imshow (imgfilter2)
    centroids2 = cat(1, s.Centroid);
    centroids2= sortrows(centroids2,2); % will sort the y's in ascending order
    numbercells(n) = size(centroids2,1);% number of rows
    
    sz2=size(centroids2);    
    
    if(numbercells(n)<=numbercells(n-1)) % if number of cells remain same no noise -> no filter needed 
        numbercells(n) = numbercells(n);
        
        for i=1:sz2(1,1)
            A (i,1,n) = centroids2(i,1);
            A (i,2,n) = centroids2(i,2);
        end
        
    else % that is noise is there, find the deviation delete it i.e. filter out the noise
         centroids2 = centroidFilter (centroids1,centroids2);
         numbercells(n) = size(centroids2,1); 
         
         sz2=size(centroids2); % size of centroid matrix changed after filtering
              
         for i=1:sz2(1,1)
             A (i,1,n) = centroids2(i,1);
             A (i,2,n) = centroids2(i,2);
         end
    end
     
     imshowpair(imgfilter2,frames2,'montage') %remove 'montage' option for overlaid comparison     
     hold on
     plot(centroids2(:,1),centroids2(:,2), 'b*')
     text(10,10,strcat('\color{red}Objects Found:',num2str(numbercells(n))))
     hold off
     time(n)=t;
     t=t+dt;
     if t<T
     vid.CurrentTime= min(t,T);
     else
         break;
     end
     n=n+1;
end
    %% plot against time
    %subplot(2,1,1); 
    %plot (time,numbercells,'r-')
    %title('Patient ID:xxxx FlowRate:yyyy')
    %xlabel('t (s)')
    %ylabel('N (Cell Count)')
    %hold off
    
    %% plot tracked cells
    %subplot(2,1,2);
    imshow (frames1)
    hold on
    xf = centroids2(:,1);
    yf = centroids2(:,2);
    plot(xf,yf,'r*',...
        'MarkerSize',15)
    hold on
    
    for f =1:c 
        x = A(:,1,f);
        y = A(:,2,f);
        hold on
        plot(x,y,'bo',...
            'MarkerSize',4)
    end
    hold off
    


%% Velocity profile : Translational velocity of rolling cells

%% METHOD 1 : velocity profile final - initial

% I have just calculated the velocity using final - initial for now, 
% we can get better value by taking averages, or other techniques
s=1;

for h = 1:c
    if nnz(A(:,:,h+1)) ~= nnz(A(:,:,h))
        length1 = nnz(A(:,:,h))/2;
        length2 = nnz(A(:,:,h+1))/2;
        for i = 1 : length1
            flag=0;
            for j = 1 : length2
                dist = ((A(i,1,h)-A(j,1,h+1))^2+(A(i,2,h)-A(j,2,h+1))^2)^0.5;
                
                if dist<10
                   flag=1; 
                   break;
                end
                
            end
                if (flag == 0)
                        pos(s,1) = A(i,1,h);
                        pos(s,2) = A(i,2,h);
                        pos(s,3) = h;
                        s=s+1;
                end
         end
     end
   
end

finaltime = ones (size(centroids2,1),1);
finaltime = c * finaltime;
cent2 = cat(2, centroids2, finaltime);
%finalpos = cat(1, cent2, pos);
%finalpos= sortrows(finalpos,2);

%for i=1:size(centroids1,1)
%    vx(i) = (centroids1 (i,1) - finalpos (i,1))/ (finalpos (i,3)*dt);
%    vy(i) = (centroids1 (i,2) - finalpos (i,2))/ (finalpos (i,3)*dt);    
%end   

%% METHOD 2 : determine rolling using continuous velocity trend
% avg of rolling velocities of all rolling cells

r=1;
q=1;
for h=1:4:c-3
    if nnz(A(:,:,h)) == nnz(A(:,:,h+4)) %if no cell has disappeared, determine if any cell is rolling
       % assumption things have only moved significantly in x direction
       for i=1: nnz(A(:,:,h))/2 
           dist = ((A(i,1,h)-A(i,1,h+4))^2+(A(i,2,h)-A(i,2,h+4))^2)^0.5;
           if (dist > 5 && dist< 25) %rolling cell
               vxx(r)= (A(i,1,h)-A(i,1,h+4))/(4*dt);
               vyy(r)= (A(i,2,h)-A(i,2,h+4))/(4*dt);
               r=r+1;
           end
       end
    else % a particular cell is gone from the frame
        % temp = pos(:,3) - h
        % [aa,bb] = min (temp)
        % use bb as the index value, but for now it is arranged
        B = A(:,:,h);% making a copy to keep the original centroid matrix intact
        for w = 1:nnz(A(:,:,h))/2 
            dist = ((A(w,1,h) - pos(q,1))^2+ (A(w,2,h) - pos(q,2))^2)^0.5;
            if dist < 5 % element found, delete it
                rowToDelete = w; 
                B(rowToDelete, :) = [];
                q=q+1;
                break;
            end
        end
           for i=1: nnz(B(:,:))/2 
           dist = ((B(i,1)-A(i,1,h+4))^2+(B(i,2)-A(i,2,h+4))^2)^0.5;
               if (dist > 5 && dist< 25) %rolling cell
                   vxx(r)= (B(i,1)-A(i,1,h+4))/(4*dt);
                   vyy(r)= (B(i,2)-A(i,2,h+4))/(4*dt);
                    r=r+1;
               end
           end

    end
    
end
mean(vxx)

%% Method 3 : Hopefully more accurate and concrete, edit : definitely better
flag = 0;
ss=0;
rr = 1;
for h=1:c
    length1 = size(centroids1,1); 
    for i =1: length1
        for j=1:length1
            tempA(j) = ((A(i,1,h) - A(j,1,h+1))^2 +(A(i,2,h) - A(j,2,h+1))^2)^0.5;
        end
        %[tempdist,index] = min(tempA) if want vx and vy -> extra work
        tempdist = min(tempA);
        if tempdist < 30
           if flag == 0 
                v(i,h) = tempdist/dt;
           else 
               if h == ha
                  v(i+ss-1,h) = tempdist/dt;
               else
                  len = length(tt);
                  
                  if(len==1)
                        if i<tt(1)
                            v(i,h) = tempdist/dt;
                        elseif i== tt(1)
                            v(i+1,h) =tempdist/dt;
                        elseif tt(1)<i
                            v(i+1,h) = tempdist/dt;                        
                        end
                        
                  elseif len > 1
                      if i <tt(1)
                          v(i,h) = tempdist/dt;
                      end
                      
                      if i ==tt(1)
                         count = 0;
                            ii=i;
                            while (nnz(ii==tt)==1)
                                count = count+1;
                                ii=ii+1;   
                            end
                            v(i+count,h) = tempdist/2;       
                      end
                      
                      for a = 1:len-1
                          if (i>tt(a) && i<tt(a+1))
                              v(i+a,h) = tempdist/dt;
                          end
                          if (i==tt(a+1))
                              count = a;
                              ii=i;
                            while (nnz(ii==tt)==1)
                                count = count+1;
                                ii=ii+1;   
                            end
                              v(i+count,h) = tempdist/2;
                          end 
                      end
                      
                      if i>tt(a+1)
                          v(i+a+1,h) =tempdist/dt;
                      end
                  end
                  
               end
           end
        else
           v(i,h) = 0; 
           tt (rr) = i+ss;
           tt= sort (tt);
           flag = 1;
           ha = h;
           rr=rr+1;
           ss=ss+1;
        end 
    end    
end

rowMean = sum(v,2) ./ sum(v~=0,2);
rowMean = rowMean(~isnan(rowMean));
velchar = num2str(rowMean);
for i = 1:size(v,1) 
    for j = 1: size (v,2)
        if(v(i,j)~=0)
            stdvel(i,j) = v(i,j);
        else 
            stdvel(i,j) =nan;
        end
    end
end



text(centroids1(:,1),centroids1(:,2),velchar,'Color','white','FontSize',8)



