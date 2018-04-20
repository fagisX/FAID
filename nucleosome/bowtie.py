import os

def bowtie(samples, name):
    fastq_dump_cmd = "/home/tmhbxx3/archive/tools/sratoolkit/bin/fastq-dump "

    bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/sacCer3/bowtie/sacCer3 '

    pbs = open(name + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N bowtie_" + name + '\n')
    pbs.write("#PBS -q highmem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=32:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    for sample_id in samples:
        pbs.write("wget https://sra-download.ncbi.nlm.nih.gov/srapub/" + sample_id + '\n')
        pbs.write(fastq_dump_cmd + os.getcwd()+'/'+sample_id + '\n')
    pbs.write('cat ' + ' '.join([s+'.fastq' for s in samples])+' >' + name + '.fastq'+'\n')
    pbs.write(bowtie_cmd + name + ".fastq " + name + ".bowtie\n")
    pbs.close()
    # os.system('qsub '+sample_id+".pbs")
    # break
    return

# samples = ['SRR023806', 'SRR023807', 'SRR023808']
# name = 'Galactose'
# bowtie(samples, name)

samples = ['SRR023800', 'SRR023801', 'SRR023802', 'SRR023803', 'SRR023804', 'SRR023805']
name = 'Dextrose'
bowtie(samples, name)

python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpos Galactose.bowtie -o Galactose -u 1 --smooth_width 0 -c 25000000 --frsz 200 --extend 200
