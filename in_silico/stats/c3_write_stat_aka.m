function c3_write_stat_aka(pop1,pop2)
%% This code finds the best model for each gene


for n1 = 1:100
    for n2 = 1:8

        filename = ['./reg10/results_aka0001/res_aka_' num2str(n1) '_' num2str(n2) '_' num2str(pop1) '_' num2str(pop2) '.txt'];


        data = dlmread(filename);



        chi = data(1:8,end);

    	[~, min_chi] = min(chi);


    	dlmwrite(['reg10/stat_aka0001_' num2str(pop1) '_' num2str(pop2) '.txt'],[n1,n2,min_chi],'-append');

    end
end


end