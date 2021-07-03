clc;
clear;
close all;
%% Define variables
step=0.100;% DEFINE TIME STEP/resolution
save = 1;  % To save the profile use 1
n=0.4;     %in degrees/s noise MAX amplitude added to a amplitude calculated as noise=(n/4)+rand(1,1)*(3*n/4)
%noise is pseudorandomized
%f_noise=0; %in HZ noise MAX frequency calculated as frequency=(f_noise/4)+rand(1,1)*(3*f_noise/4)
%% *YOU ARE ADVISED NOT TO CHANGE ANY OTHER VALUES (Except line 12 -save directory*
for p=20:24   % HOW MANY PARTICIPANTS
if save==1
% DEFINE YOUR PREFERED SAVE DIRECTORY Change both mkdir and cd (below are my computers :D)
%mkdir (['C:\Users\Administrator\OneDrive - Universiteit Utrecht\TNO\MatLab\ROOT\participant ' int2str(p)])
%cd (['C:\Users\Administrator\OneDrive - Universiteit Utrecht\TNO\MatLab\ROOT\participant ' int2str(p)])
mkdir (['X:\Dont Touch\Docs\OneDrive\OneDrive - Universiteit Utrecht\TNO\MatLab\ROOT\participant ' int2str(p)])
cd (['X:\Dont Touch\Docs\OneDrive\OneDrive - Universiteit Utrecht\TNO\MatLab\ROOT\participant ' int2str(p)])
end
set_400_500=[400:1:500];
number_to_count=set_400_500(randperm(100,24));
num_of_run=1;
%Randomize profiles
for block = 1:2

% Define the runs

if block == 1
    
    Baseline_runs = {'1.  Baseline'; '2.  Baseline'};
    Surprise_runs = {'12.  Surprise-High-left'; '11.  Surprise-Low-right'};
    Nosurprise_runs = {'8.  No surprise-High-left'; '7.  No surprise-High-right'; '3.  No surprise-Low-left'; '4.  No surprise-Low-right';...
    '9.  No surprise-High-left'; '10.  No surprise-High-right'; '5.  No surprise-Low-left'; '6.  No surprise-Low-right'};
else
    
    Baseline_runs = {'1.  Baseline'; '2.  Baseline'};
    Surprise_runs = {'12  Surp-High-right'; '11.  Surp-Low-left'};
    Nosurprise_runs = {'8.  No surprise-High-left'; '7.  No surprise-High-right'; '3.  No surprise-Low-left'; '4.  No surprise-Low-right';...
    '9.  No surprise-High-left'; '10.  No surprise-High-right'; '5.  No surprise-Low-left'; '6.  No surprise-Low-right'};
end

% Shuffle the runs

Baseline_runs = Baseline_runs(randperm(2));
Surprise_runs = Surprise_runs(randperm(2));
Nosurprise_runs = Nosurprise_runs(randperm(8));

% Create 4 sections within the block, each starting with 3 No surprise.
% Alternate the order of these sections between participants

if rem(p,2)  
    Runorder(1:12,block) = [Nosurprise_runs(1:2); Baseline_runs(1);...
    Nosurprise_runs(3:4); Surprise_runs(1);...
    Nosurprise_runs(5:6); Baseline_runs(2);...
    Nosurprise_runs(7:8); Surprise_runs(2)];
else
    Runorder(1:12,block) = [Nosurprise_runs(1:2); Surprise_runs(1);...
    Nosurprise_runs(3:4); Baseline_runs(1);...
    Nosurprise_runs(5:6); Surprise_runs(2);...
    Nosurprise_runs(7:8); Baseline_runs(2)];
end

end
%Produce profiles
for b=1:2
directR=0;
directL=0;
for i=1:12 %1-2 Baseline| 3-6 No Low| 7-10 No High| 11 Yes Low| 12 Yes High
%time0=2;
time0=round(rand(1,1)*10,0); %Dead time
%time1=45or75  % acceleration (predifined due to predifined magnitudes)
time2=0;       % constanst speed ! if vel is constant, perception decays at after 15.8s-26s(St George, 2011)
time3=2;       % diceleration to lower speed
time4=12;      % constant lower speed
time5=2;       % ILLUSION DECELERATION TIME
time6=0;       % Dead time 2
run=i;
random_run=str2double(regexprep(Runorder{i,b},'\D',''));
%random_run=6;
%i=12;
if random_run<=2
    ID="BASELINE";
    disp (">>"+ID+" dead time " + time0)
    a=25;
    time1=75;
    a_illusion= 0.5;
elseif random_run==12
    a=25;
    time1=75;
    a_illusion= 0.5;
    ID="INCO_HIGH";
    disp (">>"+ID+" dead time " + time0)
elseif random_run<=10 && random_run>6
     a=2;
    time1=75;
    a_illusion= -1.8;
    ID="CONG_HIGH";
    disp (">>"+ID+" dead time " + time0)
elseif random_run==11
    a=15;
    time1=75;
    a_illusion= 1;
    ID="INCO_LOW";
    disp (">>"+ID+" dead time " + time0)
elseif random_run<=6 && random_run>2
    a=2;
    time1=75;
    a_illusion= -1.6;
    ID="CONG_LOW";
    disp (">>"+ID+" dead time " + time0)
end
acc_t_low=0.28; 
acc_t_med=0.72; %Seemungal,2004 (M=1.18, SD=0.46)
acc_t_high=1.64;

filename=p+"_block_"+ b +"_run_"+run+"_"+ID+"_"+number_to_count(num_of_run);
num_of_run=num_of_run+1;
if save==1
diary(filename+".txt");
end
total_time=(time0+time1+time2+time3+time4+time5+time6);
%t=0:step:(total_time);
%% Zero step of the motion (dead time phase)
t0=0:step:time0;
v0(1:round(length(t0)))=0;
%% First of the motion (subthreshold ACCELERATION phase)
t1=0:step:time1;
f1=(1/time1);
v1=(-a)*(cos((pi/2)*f1*t1)-1);
noise1=((3*n/4)+(rand(1,1)*(n/4)))*(sin((2*pi)*(17.23*t1(1:end-50))));
noise2=((3*n/4)+(rand(1,1)*(n/4)))*(sin((2*pi)*(17.76*t1(1:end-50))));
noise3=((3*n/4)+(rand(1,1)*(n/4)))*(sin((2*pi)*(17.48*t1(1:end-50))));
noise4=((3*n/4)+(rand(1,1)*(n/4)))*(sin((2*pi)*(17.53*t1(1:end-50))));
noise5=((3*n/4)+(rand(1,1)*(n/4)))*(sin((2*pi)*(17.87*t1(1:end-50))));
noise=noise1+noise2+noise3+noise4+noise5;
pacifier_start=(-0.5)*(cos(pi*(1/10)*t1(1:101))-1);
pacifier(1:length(t1)-151)=1;
pacifier_end=(-0.5)*(cos(pi*(1/10)*t1(end-150:end-50)+pi)-1);
pacifier_total=[pacifier,pacifier_end];
pacifier_total(1:101)=pacifier_total(1:101).*pacifier_start;
%plot(pacifier_total)
noise=noise.*pacifier_total;

if random_run<=2
    v1(1:end)=0;
end 
v1(1:end-50)=v1(1:end-50)+noise;
%% Calculate acc (first step)
dt1=diff(t1);
dv1=diff(v1);
acc=dv1./dt1; %devide each point
Max_acc=round(max(acc),2);
disp("The Max Acceleration is " +Max_acc+"°/s^2")
disp("With Frequency " +round(f1,2)+"Hz")
%% Distance travelled (first step)
distance1=round((sum((v1)*step)),2);
%% Second step of the motion (CONSTANT SPEED)
t2=0:step:time2;
v2(1:round(length(t2)))=v1(end);
%% Distance travelled (second step)
distance2=round((sum((v2)*step)),2);
%% Third step of the motion (DECELERATION to a lower speed phase)
t3=0:step:time3;
f3=(1/time3);
v3=(-v2(end)/(2+a_illusion))*(cos((pi*f3*t3+pi))-1-a_illusion);
%% Calculate dicceleration
dt3=diff(t3);
dv3=diff(v3);
dec=dv3./dt3; %devide each point
Max_dec=round(min(dec),2);
%% Distance travelled (third step)
distance3=round((sum((v3)*step)),2);
%% Fourth of the motion (Lower constant speed - Post yaw illusion phase)
t4=(0:step:time4);
f4=(1/time4);
v4(1:round(length(t4)))=v3(end);
%% Calculate acc (fourth step/illusion)
dt4=diff(t4);
dv4=diff(v4);
acc_postyaw=dv4./dt4; %devide each point
Max_acc_postyaw=round(max(acc_postyaw),2);
%% Distance travelled (fourht step)
distance4=round((sum((v4)*step)),2);
%% Fifth step of the motion (deceleration to stop after illusion)
t5=(0:step:time5);
f5=(1/time5);
v5=((-v4(end))/2)*(cos(pi*f5*t5+pi)-1);
% v5=(-(a-(2*(a/(2+a_illusion))))/2)*(cos((pi*f5*t5)+pi)-1);
%% Zero step of the motion (2nd dead time phase)
t6=0:step:time6;
v6(1:round(length(t6)))=0;
%% Total distances
distance=distance1+distance2+distance3+distance4;
%% Total times
t1_temp=t0(end)+t1;
t2_temp=t1_temp(end)+t2;
t3_temp=t2_temp(end)+t3;
t4_temp=t3_temp(end)+t4;
t5_temp=t4_temp(end)+t5;
t6_temp=t5_temp(end)+t6;
t=[t0,t1_temp,t2_temp,t3_temp,t4_temp,t5_temp,t6_temp];
%% Total velocities
v=[v0,v1,v2,v3,v4,v5,v6];
%% displays
disp("The Max Deceleration is " +Max_dec+"°/s^2")
disp("With Frequency " +round(f3,2)+"Hz")
disp("The total displacement is " +distance+"°")
disp("With " +round(distance/360,2)+" turns")
disp("In " + total_time+"s Total time")
%% Calculate total acc
dt=diff(t);
dv=diff(v);
total_acc=dv./dt; %devide each point
total_acc(end+1) = 0;
%% Acceleration threshold in t
acc_threshold_low=zeros(size(t));
acc_threshold_low(:)=acc_t_low;
acc_threshold_med=zeros(size(t));
acc_threshold_med(:)=acc_t_med; 
acc_threshold_high=zeros(size(t));
acc_threshold_high(:)=acc_t_high; 
%% Plot total (t,v) & (acc)
figure(i);
subplot(2,1,1);
hold on 
plot(t,v,'Color', 'black','LineWidth',1); %plot time and velocity
yline(0,'Color', 'blue','LineWidth',.5, 'Linestyle','--');
ylim([-25 36]);
xlim([0 total_time+1]);
grid on;
xlabel('Time in s');
ylabel('Velocity in °/s');
subplot(2,1,2);
hold on 
%plot(t,acc_threshold_med,'Color', 'red','LineWidth',.5, 'Linestyle','--');
    %plot(t,acc_threshold_high,'Color', 'magenta','LineWidth',.5, 'Linestyle','--');
%plot(t,-acc_threshold_med,'Color', 'red','LineWidth',.5, 'Linestyle','--');
    %plot(t,-acc_threshold_high,'Color', 'magenta','LineWidth',.5, 'Linestyle','--');
plot (t, total_acc,'Color', 'black','LineWidth',1);
%legend('0.72°/s^2');
ylim([-26 26]);
grid on;
%ylim([Max_dec-2 Max_acc+2]);
xlim([0 total_time+1]);
%grid on;
xlabel('Time in s');
ylabel('Δv/Δt in °/s^2');
hold off
%subplot(3,1,3);
%plot(t1,noise)
%xlabel('Time in s');
%ylabel('Velocity in °/s');
%% Writing to Excel .csv
t=0:step:(length(v)-1)*step';
t=round(t,3);
v=round(v,3);
TVCombined=transpose(cat(1,t,v));
if save==1
saveas(gcf,filename+'.png');
if random_run<=2
writematrix(TVCombined,filename+".csv",'Delimiter',';');
end
direction=round(rand(1,1),1);
if direction>=0.5 && random_run>=3 && directR<5
writematrix(TVCombined,filename +"_R.csv",'Delimiter',';');
directR=directR+1;
elseif random_run>=3
TVCombined=transpose(cat(1,t,-v));
writematrix(TVCombined,filename+ "_L.csv" ,'Delimiter',';');
directL=directL+1;
end
end
diary off
diary("Print for Participant "+p+".xls");
disp (filename)
diary off
clearvars -except i p save step b Runorder directR directL n f_noise number_to_count num_of_run
end
close all
diary("Print for Participant "+p+".xls");
disp Block2
diary off
end
end