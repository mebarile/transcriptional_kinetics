function c2_write_aka(pop1,pop2)
%% This code computes the akaike index per each gene

dof_base = [24,32,32,32,40,40,40,48]+8;


n = 80;

for n1 = 1:100
    for n2 = 1:8

        filename = ['./reg10/results0001/gene_' num2str(n1) '_' num2str(n2) '_' num2str(pop1) '_' num2str(pop2) '.txt'];

        data = dlmread(filename);


        for line = 1:8

            chi = data(line,4);

            k = dof_base(line);


            aka = 2*k + chi + (2*k^2+2*k)/(n-k-1);

            dlmwrite(['./reg10/results_aka0001/res_aka_' num2str(n1) '_' num2str(n2) '_' num2str(pop1) '_' num2str(pop2) '.txt'],[line,chi,k,aka],'-append');


        end
    end
end


