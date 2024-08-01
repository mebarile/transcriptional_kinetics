function c4_write_para(pop1,pop2)
%% This code writes the parameters estimated for each gene


for n1 = 1:100
    for n2 = 1:8


        for model = 1:8

            filename = ['./reg10/data_output/parameters_'  num2str(n1) '_' num2str(n2) '_model' num2str(model)  '_reg10_' num2str(pop1) '_' num2str(pop2) '.mat'];
            data = load(filename);
            para = data.parameters.MS.par(:,1)';


            dlmwrite(['./reg10/parameters_.txt'],[n1,n2,model,para],'-append')

        end

    end
end
