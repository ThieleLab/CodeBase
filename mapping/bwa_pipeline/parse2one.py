import re
import os
import glob
from Bio import SeqIO
os.chdir('P:/MAPPING/Genomes')

#fo = open("joined_contigs.txt", "wb")
#for i in glob.iglob('*.fasta'):
#    fasta_sequences = SeqIO.parse(open(i),'fasta')
#    orgnam = i.replace('.fasta','')
#    print orgnam
#    for fasta in fasta_sequences:
#        name = fasta.description
#        fo.write('>' + orgnam + ' _:_ ' + name + '\n');
#        fo.write(str(fasta.seq) + '\n\n');
        
        #print name
    #print orgnam
#fo.close()

print 'Starting...'

fo = open("P:/MAPPING/Genomes/joint/joined_genomes.txt", "wb")
for i in glob.iglob('*.fasta'):
    fasta_sequences = SeqIO.parse(open(i),'fasta')
    orgnam = i.replace('.fasta','')
    orgseq = ''
    print orgnam
    for fasta in fasta_sequences:
        orgseq = orgseq + str(fasta.seq) + 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
       #print name
        #print orgnam
    fo.write('>' + orgnam + '\n');
    fo.write(orgseq[1:(len(orgseq)-50)] + '\n');
fo.close()

print 'Finished'
