import peak_clusters as pk
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord as SR
import pdb
import time


class differential_motif:
    
    def __init__(self):
        seq_list = []
        p = pk.peak_clusters()
        
        for i in range(10):
            seq_list.append(p.form_cluster_sequences(i))
        
        self.seqs = seq_list
            
        
    def differential(self, index):
        seqs = self.seqs
        rseqs = seqs[index]
        
        cseqs = []
        for i in range(len(seqs)):
            if not i == index:
                cseqs = cseqs + seqs[i]
                
        rSR = [SR(Seq(rseqs[i]), id="response_str" + str(i), description="") 
               for i in range(len(rseqs))]
        cSR = [SR(Seq(cseqs[i]), id="control_str" + str(i), description="") 
               for i in range(len(cseqs))]
        
        SeqIO.write(rSR, "response.fasta", "fasta")
        SeqIO.write(cSR, "control.fasta", "fasta")
        
    def meme_ame(self, motif_file,
                 ame_dir="/Users/sr286/meme/"):
        
        start = time.time()

        bpath = ame_dir + "bin"
        lpath = ame_dir + "libexec/meme-5.1.0"
        path = os.environ["PATH"]
       
        os.environ["PATH"] = bpath + ":" + lpath + ":" + path
        os.environ["mf"] = motif_file
        os.system("ame --control control.fasta response.fasta $mf")
        
        end = time.time()
        print(end - start)
        # MEME-AME is too slow!  Let's use homer!

 
        
       