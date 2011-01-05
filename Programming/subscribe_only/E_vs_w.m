radius = [20,18,16,14,12,10];
Num  = 5;

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h1_K2_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h1_k2_e(i) = energy;
  std_h1_k2_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h2_K2_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h2_k2_e(i) = energy;
  std_h2_k2_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h3_K2_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h3_k2_e(i) = energy;
  std_h3_k2_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h1_K3_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h1_k3_e(i) = energy;
  std_h1_k3_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h2_K3_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h2_k3_e(i) = energy;
  std_h2_k3_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h3_K3_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h3_k3_e(i) = energy;
  std_h3_k3_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h1_K4_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h1_k4_e(i) = energy;
  std_h1_k4_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h2_K4_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h2_k4_e(i) = energy;
  std_h2_k4_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h3_K4_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h3_k4_e(i) = energy;
  std_h3_k4_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h1_KG_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h1_kg_e(i) = energy;
  std_h1_kg_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h2_KG_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h2_kg_e(i) = energy;
  std_h2_kg_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\STD_',num2str(radius(i)),'_840_h3_KG_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  std_h3_kg_e(i) = energy;
  std_h3_kg_w(i) = 2*pi/t_star;
end




for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h1_K2_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h1_k2_e(i) = energy;
  ISO_h1_k2_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h2_K2_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h2_k2_e(i) = energy;
  ISO_h2_k2_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h3_K2_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h3_k2_e(i) = energy;
  ISO_h3_k2_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h1_K3_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h1_k3_e(i) = energy;
  ISO_h1_k3_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h2_K3_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h2_k3_e(i) = energy;
  ISO_h2_k3_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h3_K3_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h3_k3_e(i) = energy;
  ISO_h3_k3_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h1_K4_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h1_k4_e(i) = energy;
  ISO_h1_k4_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h2_K4_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h2_k4_e(i) = energy;
  ISO_h2_k4_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h3_K4_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h3_k4_e(i) = energy;
  ISO_h3_k4_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h1_KG_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h1_kg_e(i) = energy;
  ISO_h1_kg_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h2_KG_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h2_kg_e(i) = energy;
  ISO_h2_kg_w(i) = 2*pi/t_star;
end

for i = 1:Num
  file = ['C:\School\Programming\subscribe_only\ISO_',num2str(radius(i)),'_840_h3_KG_r3g_targ10.eph'];
  A = load(file);
  [t_star,energy] = find_period(A);
  ISO_h3_kg_e(i) = energy;
  ISO_h3_kg_w(i) = 2*pi/t_star;
end

clear A;
clear i;
clear t_star;
clear energy;
clear file;
clear radius;

