
class Precursor():

    def __init__(self, name, mirna1, mirna2, precursor):

        # The prefix of the filename for the pre-miRNA
        self.name = name

        # The precursor sequence
        self.premirna = precursor
        
        # The tuple containing the positions of the miRNA within the pre-miRNA
        self.pos1 = self.pos(precursor, mirna1)

        # The tuple containing the positions of the other miRNA within the pre-miRNA
        self.pos2 = self.pos(precursor, mirna2)

        # The length of the pre-miRNA
        self.prelen = len(precursor)

        # The Minimum Free Energy for the precursor as predicted by RNAfold
        self.premfe = 0


    def pos(self, premirna, mirna):

        if mirna:

            if mirna in premirna:

                if premirna.count(mirna) > 1:
                    print('WARNING! miRNA {} was found more than once in the precursor sequence {}..., but its last occurrence will be used!'.format(mirna, premirna[25:]))

                return (premirna.find(mirna) + 1, premirna.find(mirna) + len(mirna))
            
            else:
                
                raise Exception('ERROR! Could not find sequence {} inside {}, please correct this'.format(mirna, premirna))

        else:

            return None

