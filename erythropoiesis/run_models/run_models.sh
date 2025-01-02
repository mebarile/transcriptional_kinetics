# run this line in the command line, changing the input: 
# gene_number: taken from a list
# reg: regularization parameter
# pop1: usually 1 unless you have more than 2 conditions
# pop2: usually 2 unless you have more than 2 conditions

matlab -nodisplay -nojvm -r "main_model1($gene_number,$reg,$pop1,$pop2)"

