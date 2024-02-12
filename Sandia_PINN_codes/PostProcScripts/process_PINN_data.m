function [] = process_PINN_data(Pe)
figure(); 
a = csvread('ParamResults_FOM_None_No_Text.csv'); 
[i,v] = find(a(:,2) == Pe);
data = a(i,:);
[i_2sd_noSnap_weakBC,v] = find((data(:,3) == 0) & (data(:,4) == 0) & (data(:,5) == 2) ...
    & (data(:,8) < 20) );
cpu = data(i_2sd_noSnap_weakBC,7); 
for i=1:length(i_2sd_noSnap_weakBC)
    mse(i) = mean(data(i_2sd_noSnap_weakBC(i),9:end));
end
if (length(cpu) > 0)
  loglog(cpu, mse, 'ro'); 
end
[i_3sd_noSnap_weakBC,v] = find((data(:,3) == 0) & (data(:,4) == 0) & (data(:,5) == 3) ...
    & (data(:,8) < 20) );
cpu = data(i_3sd_noSnap_weakBC,7);
clearvars mse; 
for i=1:length(i_3sd_noSnap_weakBC)
    mse(i) = mean(data(i_3sd_noSnap_weakBC(i),9:end));
end
hold on;
if (length(cpu) > 0) 
  loglog(cpu, mse, 'bo'); 
end
[i_4sd_noSnap_weakBC,v] = find((data(:,3) == 0) & (data(:,4) == 0) & (data(:,5) == 4) ...
    & (data(:,8) < 20) );
cpu = data(i_4sd_noSnap_weakBC,7); 
clearvars mse; 
for i=1:length(i_4sd_noSnap_weakBC)
    mse(i) = mean(data(i_4sd_noSnap_weakBC(i),9:end));
end
hold on;
if (length(cpu) > 0) 
  loglog(cpu, mse, 'go'); 
end
[i_5sd_noSnap_weakBC,v] = find((data(:,3) == 0) & (data(:,4) == 0) & (data(:,5) == 5) ...
    & (data(:,8) < 20) );
cpu = data(i_5sd_noSnap_weakBC,7); 
for i=1:length(i_5sd_noSnap_weakBC)
    mse(i) = mean(data(i_5sd_noSnap_weakBC(i),9:end));
end
hold on;
if (length(cpu) > 0) 
  loglog(cpu, mse, 'ko'); 
end

[i_2sd_noSnap_strongBC,v] = find((data(:,3) == 0) & (data(:,4) == 1) & (data(:,5) == 2)   ...
    & (data(:,8) < 20) );
cpu = data(i_2sd_noSnap_strongBC,7); 
for i=1:length(i_2sd_noSnap_strongBC)
    mse(i) = mean(data(i_2sd_noSnap_strongBC(i),9:end));
end
if (length(cpu) > 0) 
  loglog(cpu, mse, 'rv'); 
end
[i_3sd_noSnap_strongBC,v] = find((data(:,3) == 0) & (data(:,4) == 1) & (data(:,5) == 3)  ...
    & (data(:,8) < 20) );
cpu = data(i_3sd_noSnap_strongBC,7); 
clearvars mse; 
for i=1:length(i_3sd_noSnap_strongBC)
    mse(i) = mean(data(i_3sd_noSnap_strongBC(i),9:end));
end
hold on; 
if (length(cpu) > 0) 
  loglog(cpu, mse, 'bv'); 
end
[i_4sd_noSnap_strongBC,v] = find((data(:,3) == 0) & (data(:,4) == 1) & (data(:,5) == 4)  ...
    & (data(:,8) < 20) );
cpu = data(i_4sd_noSnap_strongBC,7); 
clearvars mse; 
for i=1:length(i_4sd_noSnap_strongBC)
    mse(i) = mean(data(i_4sd_noSnap_strongBC(i),9:end));
end
hold on;
if (~isempty(cpu)) 
  loglog(cpu, mse, 'gv'); 
end
[i_5sd_noSnap_strongBC,v] = find((data(:,3) == 0) & (data(:,4) == 1) & (data(:,5) == 5)  ...
    & (data(:,8) < 20) );
cpu = data(i_5sd_noSnap_strongBC,7); 
clearvars mse; 
for i=1:length(i_5sd_noSnap_strongBC)
    mse(i) = mean(data(i_5sd_noSnap_strongBC(i),9:end));
end
hold on;
if (length(cpu) > 0) 
  loglog(cpu, mse, 'kv'); 
end

[i_2sd_snap_weakBC,v] = find((data(:,3) == 1) & (data(:,4) == 0) & (data(:,5) == 2)  ...
    & (data(:,8) < 20) );
cpu = data(i_2sd_snap_weakBC,7); 
for i=1:length(i_2sd_snap_weakBC)
    mse(i) = mean(data(i_2sd_snap_weakBC(i),9:end));
end
if (length(cpu) > 0) 
  loglog(cpu, mse, 'r^'); 
end
[i_3sd_snap_weakBC,v] = find((data(:,3) == 1) & (data(:,4) == 0) & (data(:,5) == 3)  ...
    & (data(:,8) < 20) );
cpu = data(i_3sd_snap_weakBC,7); 
clearvars mse; 
for i=1:length(i_3sd_snap_weakBC)
    mse(i) = mean(data(i_3sd_snap_weakBC(i),9:end));
end
hold on;
if (length(cpu) > 0) 
  loglog(cpu, mse, 'b^'); 
end
[i_4sd_snap_weakBC,v] = find((data(:,3) == 1) & (data(:,4) == 0) & (data(:,5) == 4)  ...
    & (data(:,8) < 20) );
cpu = data(i_4sd_snap_weakBC,7); 
clearvars mse; 
for i=1:length(i_4sd_snap_weakBC)
    mse(i) = mean(data(i_4sd_snap_weakBC(i),9:end));
end
hold on;
if (length(cpu) > 0)
  loglog(cpu, mse, 'g^'); 
end
[i_5sd_snap_weakBC,v] = find((data(:,3) == 1) & (data(:,4) == 0) & (data(:,5) == 5)  ...
    & (data(:,8) < 20) );
cpu = data(i_5sd_snap_weakBC,7); 
clearvars mse; 
for i=1:length(i_5sd_snap_weakBC)
    mse(i) = mean(data(i_5sd_snap_weakBC(i),9:end));
end
hold on;
if (length(cpu) > 0) 
  loglog(cpu, mse, 'k^'); 
end

[i_2sd_snap_strongBC,v] = find((data(:,3) == 1) & (data(:,4) == 1) & (data(:,5) == 2)  ...
    & (data(:,8) < 20) );
cpu = data(i_2sd_snap_strongBC,7); 
for i=1:length(i_2sd_snap_strongBC)
    mse(i) = mean(data(i_2sd_snap_strongBC(i),9:end));
end
if (length(cpu) > 0)
  loglog(cpu, mse, 'rx');
end
[i_3sd_snap_strongBC,v] = find((data(:,3) == 1) & (data(:,4) == 1) & (data(:,5) == 3)  ...
    & (data(:,8) < 20) );
cpu = data(i_3sd_snap_strongBC,7); 
clearvars mse; 
for i=1:length(i_3sd_snap_strongBC)
    mse(i) = mean(data(i_3sd_snap_strongBC(i),9:end));
end
hold on;
if (length(cpu) > 0) 
  loglog(cpu, mse, 'bx'); 
end
[i_4sd_snap_strongBC,v] = find((data(:,3) == 1) & (data(:,4) == 1) & (data(:,5) == 4)  ...
    & (data(:,8) < 20) );
if (length(cpu) > 0) 
  cpu = data(i_4sd_snap_strongBC,7); 
end
clearvars mse; 
for i=1:length(i_4sd_snap_strongBC)
    mse(i) = mean(data(i_4sd_snap_strongBC(i),9:end));
end
hold on;
if (length(cpu) > 0) 
  loglog(cpu, mse, 'gx'); 
end
[i_5sd_snap_strongBC,v] = find((data(:,3) == 1) & (data(:,4) == 1) & (data(:,5) == 5)  ...
    & (data(:,8) < 20) );
cpu = data(i_5sd_snap_strongBC,7); 
for i=1:length(i_5sd_snap_strongBC)
    mse(i) = mean(data(i_5sd_snap_strongBC(i),9:end));
end
hold on;
if (~isempty(cpu)) 
  loglog(cpu, mse, 'kx'); 
end
xlabel('CPU time (s)'); 
ylabel('Average MSE over all Domains');
title(['Pe = ', num2str(Pe)]); 