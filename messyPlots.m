%% 6,2000
close all
ts = 0:0.001:2;
freq = 6;
mag = 2000;

sig= mag*cos(ts*2*pi*freq);
beta = 0.8;
alpha = 0.2;
r = mag*0.4;
noise = 500;
for k = 2:length(sig)
    noise(k) = alpha*noise(k-1)+beta*r*rand(1,1);
end
k=plot(ts,sig+noise);
set(gcf,'Position',[1 22 1920 978])
xkcdify(gca)
set(k,'LineWidth',2)
set(k,'Color','k')
set(gcf,'Position',[-1454         489        1357         282])
ylim([-6000 6000])
%% 6,5000
figure
ts = 0:0.001:2;
freq = 6;
mag = 5000;

sig= mag*cos(ts*2*pi*freq);
beta = 0.8;
alpha = 0.2;
r = mag*0.4;
noise = 500;
for k = 2:length(sig)
    noise(k) = alpha*noise(k-1)+beta*r*rand(1,1);
end
k=plot(ts,sig+noise);
set(gcf,'Position',[1 22 1920 978])

xkcdify(gca)
set(k,'LineWidth',2)
set(k,'Color','k')
set(gcf,'Position',[-1454         489        1357         282])
ylim([-6000 6000])

%% 10,2000
figure
ts = 0:0.001:2;
freq = 10;
mag = 2000;

sig= mag*cos(ts*2*pi*freq);
beta = 0.8;
alpha = 0.2;
r = mag*0.4;
noise = 500;
for k = 2:length(sig)
    noise(k) = alpha*noise(k-1)+beta*r*rand(1,1);
end
k=plot(ts,sig+noise);
set(gcf,'Position',[1 22 1920 978])
xkcdify(gca)
set(k,'LineWidth',2)
set(k,'Color','k')
set(gcf,'Position',[-1454         489        1357         282])
ylim([-6000 6000])

%% 10, 5000
figure
ts = 0:0.001:2;
freq = 10;
mag = 5000;

sig= mag*cos(ts*2*pi*freq);
beta = 0.8;
alpha = 0.2;
r = mag*0.4;
noise = 500;
for k = 2:length(sig)
    noise(k) = alpha*noise(k-1)+beta*r*rand(1,1);
end
k=plot(ts,sig+noise);
set(gcf,'Position',[1 22 1920 978])

xkcdify(gca)
set(k,'LineWidth',2)
set(k,'Color','k')
set(gcf,'Position',[-1454         489        1357         282])
ylim([-6000 6000])

savefigs('messy',1,1)