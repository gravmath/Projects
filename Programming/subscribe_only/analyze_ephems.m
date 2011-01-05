clear all;
diary on;

mc = ['STD','ISO'];
hs = ['h1','h2','h3'];
ks = ['K2','K3','K4','KG'];
rs = ['20','18','16','14','12','10'];

parms.M             = 1;
parms.dt            = 0.1;
parms.freq          = 50;
parms.choice        = 1;
num_orbs            = 3;
num_orbs_run        = 10;
fid = fopen('results.dat','w');


for i = 1:3:max(size(mc))
   for j = 1:2:max(size(rs))
      for k = 1:2:max(size(ks))
         for m = 1:2:max(size(hs))
            clear ephem;
            filename  = [mc(i:i+2),'_',rs(j:j+1),'_840_',hs(k:k+1),'_',ks(m:m+1),'_r3g_targ10.eph'];
            filename  = ['e:\kd_runs\',mc(i:i+2),'_',rs(j:j+1),'_1000_',ks(k:k+1),'_',hs(m:m+1),'_r3g_kd.eph'];
            cur_eph   = load(filename);
            parms.num = ( max(size(cur_eph)) - 2)/10*3;
            R         = sqrt( cur_eph(:,2).*cur_eph(:,2) + cur_eph(:,3).*cur_eph(:,3) + cur_eph(:,4).*cur_eph(:,4) ); 
            mR        = mean(R);
            sR        = std(R);
            E         = cur_eph(:,8);
            mE        = mean(E);
            if ( i == 1 );
               parms.metric_choice = 'STD_SCHW';
               v                   = 1/sqrt(mR - 3);
               parms.state         = [mR;0;0;0;v;0];
            end
            if ( i == 4 );
               parms.metric_choice = 'ISO_SCHW';
               v                   = sqrt( (2*mR + 1)^4 / mR^3 / (-8*mR + 1 +4*mR^2 ) ) /2;
               parms.state         = [mR;0;0;0;v;0];
            end
            ephem  = gr_geos_v2(parms);
            anal_E = ephem.point(2).H;
            save analytic.mat ephem;
            deg_orbit = comp_anal_and_fp('analytic.mat',filename,parms.metric_choice,1);
            fprintf(fid,'%s\t%s\t%s\t%s\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',mc(i:i+2),rs(j:j+1),hs(k:k+1),...
                    ks(m:m+1),mR,sR,deg_orbit,anal_E,mE);
         end
      end
   end
end
close(fid);
diary off;