function c1_analyse_genes
%% this code puts together the results for each gene

addpath(genpath('./matlab_simulations'))

reg = 10;

pop1 =1;
pop2 =2;

for n1 = 1:100
    for n2 = 1:8
        for model = 1:8

            filename = ['./reg10/data_output/parameters_' num2str(n1) '_' num2str(n2) '_model' num2str(model) '_reg' num2str(reg) '_' num2str(pop1) '_' num2str(pop2) '.mat'];

            data  = load(['./data_input/data_gene_' num2str(n1) '_model' num2str(n2) '.mat']);


            data_input = data.D.pop.mean;
            data_input = data_input(:,[(pop1-1)*2+1,pop1*2,(pop2-1)*2+1,pop2*2]);


            t = data.D.pop.t;

            err_regulation = 0.001*max(data_input);
            err_regulation  = repmat(err_regulation,20,1);

            para = load(filename);


            name_m =  ['simulate_kinetics_model' num2str(model)];

            modelfun = str2func(name_m);


            theta = para.parameters.MS.par(:,1);


            npara = length(theta)-4;

            [~,~,~,N,~,~] = modelfun(t,theta,[],[],[]);

            ind_w(1,1) = 1;
            ind_k(1,1) = 1;

            time_all = [data.D.pop.tw,data.D.pop.tm];

            for iHelp = 2:20

                ind_w(iHelp,1) = find(data.D.pop.t == time_all(iHelp,pop1));
                ind_k(iHelp,1) = find(data.D.pop.t == time_all(iHelp,pop2));

            end




            N = [N(ind_w,1:2),N(ind_k,3:4)];

            negLogL = 0;

            diff =0;

            for it = 1:20
                for col = 1:4

                    erri = exp(theta(npara+col));

                    negLogL = negLogL + 0.5*log(2*pi*erri^2) + ...
                        0.5*(data_input(it,col)-N(it,col)).^2/...
                        (erri^2);



                    diff = diff + (data_input(it,col)-N(it,col)).^2/(err_regulation(it,col)^2);

                end
            end

            dlmwrite(['./reg10/results0001/gene_' num2str(n1) '_' num2str(n2) '_' num2str(pop1) '_' num2str(pop2) '.txt'],[model,para.parameters.MS.logPost(1),negLogL,diff],'-append')

        end
    end
end

