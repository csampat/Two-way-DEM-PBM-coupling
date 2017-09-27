clear
close all
clc

%cd /home/chai/Documents/One_way_PBM/fixed_256/comp_16/mpi_1/;
cores = [128]; diameter = [200]; % omp = [1;2;4;8];
lc = length(cores); ld = length(diameter);
%mpicount = [2;4;8;16];
%for j = 1:length(omp)
    %a = num2str(mpicount(j));
    %b = strcat('omp_',num2str(omp(j)),'/csvDump/');
    b = 'csvDump/';
    cd(b);
    for i = 1:lc
        for k = 1:ld 
            filen = strcat('d50_',num2str(cores(i)),'_',num2str(diameter(k)),'.csv');
            q = importdata(filen);
            dat = q.data();
            dat(isnan(dat))=0.0;
            len = length(dat(1,:));
            figure
            for m = 3:len
                plot(dat(:,2),dat(:,m));
                hold on;
            end
            filena = strtok(filen,'.')
            title(filena)
            xlabel('Time (s)')
            ylabel('d-50 (mm)')
            print(filena,'-dpdf','-bestfit')
            close
        end
    end
    cd ../..;
 %end
