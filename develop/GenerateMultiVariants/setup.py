# python /Users/pgreisen/pythonscripts/develop/GenerateMultiVariants/generate_mvariant_fasta.py -f csa.fasta -n L,V -p 356,379 -m F,L
# vars.txt
import sys
vars = sys.argv[1]

gen = []

with open(vars,'r') as f:
    for line in f:
        tmp_ = line.strip().split('_')
        n_ = ""
        p_ = ""
        m_ = ""
        for j in tmp_:        
            n_ += j[0]+","
            m_ += j[-1]+","
            p_ += j[1:-1]+","
        template_ = "python /Users/pgreisen/pythonscripts/develop/GenerateMultiVariants/generate_mvariant_fasta.py -f csa.fasta -n "+n_[0:-1]+" -p "+p_[0:-1]+" -m "+m_[0:-1]
        gen.append(template_)
with open('gen.sh','w') as f:
    for line in gen:
        f.write(line+"\n")
