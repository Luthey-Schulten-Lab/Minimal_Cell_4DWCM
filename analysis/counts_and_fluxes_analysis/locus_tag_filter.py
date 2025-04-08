import numpy as np
import os

class locus_tag_info:

    def __init__(self):

        self.valid_info = False
        self.info = dict()

        return

    def load_info_file(self,
                       locus_tags,
                       info_file):

        self.info = dict()

        if os.path.isfile(info_file) == False:

            raise ValueError('File does not exist')

        with open(info_file, 'r') as f:

            genes = np.loadtxt(f,
                              delimiter=',',
                              usecols=(0),
                              dtype=np.str_)

            f.seek(0)
            
            gene_lengths = np.loadtxt(f,
                                      delimiter=',',
                                      usecols=(1),
                                      dtype=np.int32)

            f.seek(0)
            
            gene_tags = np.loadtxt(f,
                                   delimiter=',',
                                   usecols=(2),
                                   dtype=np.str_)

        for i in range(genes.shape[0]):
            genes[i] = genes[i][-4:]

        lt_array = np.array(locus_tags,dtype=np.str_)
        overlap = np.intersect1d(lt_array,genes)

        if overlap.shape[0] != len(locus_tags):

            missing = np.setdiff1d(lt_array,genes)
            print(missing)
            raise ValueError('Incompatible info - missing loci')

        for i in range(genes.shape[0]):

            self.info[genes[i]] = dict()

            self.info[genes[i]]['length'] = gene_lengths[i]
            self.info[genes[i]]['tags'] = (gene_tags[i].tolist()).split(':')

        self.valid_info = True
            
        return

    def info_test(self):

        if self.valid_info == False:

            raise ValueError('No info')

        return

    def get_tagged_loci(self,
                        tag):

        self.info_test()

        tagged_loci = []

        for key in self.info.keys():

            if any(t == tag for t in self.info[key]['tags']):

                tagged_loci.append(key)

        return tagged_loci

    def get_lengths(self,
                    loci):

        self.info_test()

        lengths = np.zeros((len(loci)),dtype=np.int32)

        for i in range(len(loci)):

            lengths[i] = self.info[loci[i]]['length']

        return lengths
        
