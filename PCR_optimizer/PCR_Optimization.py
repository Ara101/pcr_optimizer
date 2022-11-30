# Authored Lily Torp, K. Lionel Tukei 
# Fall 2022

class pcr(object): 
    def __init__ (self, gene, forward_primer, reverse_primer, template_type = "plasmid"):
        """
        gene: target gene being amplified
        forward_primer: primer used for forward reaction ***should only be portion of primer that overlaps to gene***
        reverse_primer: primer used for reverse reaction ***should only be portion of primer that overlaps to gene***
        template_type: type of DNA being amplified, either "plasmid" or "genomic

        """
        self.sequence = gene
        self.length = len(gene)
        self.fp = forward_primer
        self.rp = reverse_primer
        self.template_type = template_type
      
    def countDNAbase(self): 
        """
        Counts the number of each DNA base present in a sequence

        Parameters
        ----------
        dna: input gene 

        Returns
        ----------
        The number of each base 
        GC content : str

        """
        gene = self.sequence 
        
        a_count = 0
        c_count = 0
        g_count = 0
        t_count = 0
        for char in gene: 
          if char == "a": 
            a_count += 1 
          if char == "c": 
            c_count += 1 
          if char == "g": 
            g_count += 1 
          if char == "t": 
            t_count += 1
        gc_content = 100*(g_count + c_count)/(g_count + c_count + a_count + t_count)
        return "The GC Content is {}".format(gc_content), gc_content

    def checkPrimerGeneCompatability(self, startr = 0, stopr = 0 , startf = 0, stopf = 0): 
        """
        Check the forward and reverse primer's compatibility to the gene/sequence.
        If the start positions are not at the begining on the gene you can enter the 
        start and stop positions in acordance with Primer-blast
        
        Parameters
        ----------
        gene: str
        forward_primer: str 
        reverse_primer: str
        startr: int
        stopr: int
        startf: int
        stopf: int
        
        Returns
        -------
        A statment on primer compatibility : str
        """
        
        forward_primer = self.fp
        reverse_primer = self.rp
        gene = self.sequence

        # Preparing inputs 
        gene = gene.lower()
        forward_primer = forward_primer.lower()
        reverse_primer = reverse_primer.lower()

        # Create complement squence for comparison with forward_primer
        complement = gene.upper().replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")

        # Locating the start and stop pointsoption
        if startr == 0 and stopr == 0:
          gene_reverse = gene[-len(reverse_primer):]
        else:
          gene_reverse = gene[stopr-1:startr]
        if startf == 0 and stopf == 0:
          gene_forward = complement[:len(forward_primer)]
        else:
          gene_forward = complement[startf-1:stopf]

        # Create gene comparison points
        gene_reverse = gene_reverse[::-1] # Reverse Primers are entered in Reverse so this is ment to flip the gene to compare

        # Ensure gene and primer compatibility
        position_gene = 0
        Primer1 = 0
        Primer2 = 0
        for base_gene in gene_forward:
          if base_gene == "a" and forward_primer[position_gene] == "t":
            position_gene += 1
            continue
          elif base_gene == "t" and forward_primer[position_gene] == "a":
            position_gene += 1
            continue
          elif base_gene == "c" and forward_primer[position_gene] == "g":
            position_gene += 1
            continue
          elif base_gene == "g" and forward_primer[position_gene] == "c":
            position_gene += 1
            continue
          else: 
            Primer1 = 1
            break
            
        position_gene = 0
        for base_gene in gene_reverse:
          if base_gene == "a" and reverse_primer[position_gene] == "t":
            position_gene += 1
            continue
          elif base_gene == "t" and reverse_primer[position_gene] == "a":
            position_gene += 1
            continue
          elif base_gene == "c" and reverse_primer[position_gene] == "g":
            position_gene += 1
            continue
          elif base_gene == "g" and reverse_primer[position_gene] == "c":
            position_gene += 1
            continue
          else: 
            Primer2 = 1
            break
            
        if Primer1 + Primer2 == 2:
          return ("Both Primers and gene are incompatible")
        elif Primer1 == 1:
          return ("Forward Primer and gene are incompatible")   
        elif Primer2 == 1:
          return ("Reverse Primer and gene are incompatible")
        else:
          return ("Primers and Gene are compatible")

    def iProofAnalyzer(self): 
        """
        Acquires iProof DNA Polymerase Results 
        Based on protocol from iProof: "https://www.bio-rad.com/webroot/web/pdf/lsr/literature/10002298B.pdf"

        Parameters: 
        ----------
        a pcr object (contains gene seq, f primer, and r primer)

        Returns: 
        ----------
        enzyme_df : DataFrame
          enzyme volume 
          enzyme cost
        stats : DataFrame 
          Annealing temperature
          Annealing time 
          Extension time 
          Total PCR reaction time

        """
        #Gene length and elongation time
        #For iProof, extension time depends on the type of template DNA 
        #Protocol recommends 15 sec per kb if plasmid, and 30 sec per kb if genomic
        from Bio.SeqUtils import MeltingTemp as mt
        from Bio.Seq import Seq
        import pandas as pd

        gene_length = len(self.sequence)
        if self.template_type == "plasmid" or self.template_type == "lambda" or self.template_type == "BAC DNA": 
          extension_time_seconds = (gene_length/1000)*15
          extension_time_minutes = extension_time_seconds/60 
        elif self.template_type == "genomic": 
          extension_time_seconds = round((gene_length/1000)*30, 2)
          extension_time_minutes = round(extension_time_seconds/60, 2)

        #iProof Enzyme 
        #iProof can range from 0.5-2 units per 50mL 
        #Chose the 4 most common PCR volumes. With a higher volume you get more PCR product 
        
        enzyme_dict = {"20 uL": [], "20 uL Cost" : [], "50 uL": [], "50 uL Cost" : [], "100 uL" : [], 
                    "100 uL Cost": [], "200 uL" : [], "200 uL Cost" : []} # Store cost for each volume
        enzyme_cost = 1.53

        for volume in [20,50,100,200]: 
          enzyme_dict[str(volume) + " uL"].append(0.5*(volume/50)) 
          enzyme_dict[str(volume) + " uL Cost"].append(0.5*(volume/50)*enzyme_cost)
        enzyme_df = pd.DataFrame(enzyme_dict)
        enzyme_df.index = ['Enzyme amount']
        
        #iProof Primer Annealing Temperature
        annealing_temp_1 = mt.Tm_NN(Seq(self.fp))
        annealing_temp_2 = mt.Tm_NN(Seq(self.rp))
        if len(self.fp) or len(self.fp) > 20: 
          if annealing_temp_1 > annealing_temp_2: 
            annealing_temp_final = round(annealing_temp_2 + 3, 2)
          elif annealing_temp_1 < annealing_temp_2: 
            annealing_temp_final = round(annealing_temp_1 +3, 2)

        elif len(self.fp) or len(self.rp) <=20: 
          if annealing_temp_1 > annealing_temp_2: 
            annealing_temp_final = round(annealing_temp_2,2)
          elif annealing_temp_1 < annealing_temp_2: 
            annealing_temp_final = round(annealing_temp_1,2)
        annealing_time_seconds = 30 #this is the same regardless of the input
        annealing_time_minutes = 0.5 

        total_pcr_time_minutes = round((30 + (10 + annealing_time_seconds + extension_time_seconds)*35 + 10*60)/60, 2)
        total_pcr_time_hours = round (total_pcr_time_minutes/60, 2)
        
        factor = ['Annealing Temperature:','Annealing Time:','Extention Time:','Total PCR reaction time is:' ]
        stats_data = {'Result': ["{} degrees Celcius".format(annealing_temp_final), "{} seconds or {} minutes".format(annealing_time_seconds, annealing_time_minutes), "{} seconds or {} minutes".format(extension_time_seconds, extension_time_minutes), "{} minutes or {} hours".format(total_pcr_time_minutes, total_pcr_time_hours)]}
        stats = pd.DataFrame(data = stats_data, index = factor)

        return enzyme_df, stats 
      
    def taqAnalyzer (self): 
        """
        Acquires Taq DNA Polymerase Results 
        Based on protocol from Taq: "https://www.neb.com/protocols/0001/01/01/taq-dna-polymerase-with-standard-taq-buffer-m0273"

        Parameters: 
        ----------
        a pcr object (contains gene seq, f primer, and r primer)

        Returns: 
        ----------
        enzyme_df : DataFrame
          enzyme volume 
          enzyme cost
        stats : DataFrame 
          Annealing temperature
          Annealing time 
          Extension time 
          Total PCR reaction time

        """
        #Gene length and elongation time
        #For iProof, extension time depends on the type of template DNA 
        #Protocol recommends 15 sec per kb if plasmid, and 30 sec per kb if genomic
        from Bio.SeqUtils import MeltingTemp as mt
        from Bio.Seq import Seq
        import pandas as pd

        gene_length = len(self.sequence)
        if self.template_type == "plasmid" or self.template_type == "lambda" or self.template_type == "BAC DNA": 
          extension_time_seconds = (gene_length/1000)*15
          extension_time_minutes = extension_time_seconds/60 
        elif self.template_type == "genomic": 
          extension_time_seconds = round((gene_length/1000)*30, 2)
          extension_time_minutes = round(extension_time_seconds/60, 2)

        #Taq DNA Polymerase Enzyme 
        #Taq DNA Polymerase is 0.25 units per 50 
        #Chose the 4 most common PCR volumes. With a higher volume you get more PCR product 
        
        enzyme_dict = {"20 uL": [], "20 uL Cost" : [], "50 uL": [], "50 uL Cost" : [], "100 uL" : [], 
                    "100 uL Cost": [], "200 uL" : [], "200 uL Cost" : []} # Store cost for each volume
        enzyme_cost = 1.78 # dollar cost per reaction
        

        for volume in [20,50,100,200]: 
          enzyme_dict[str(volume) + " uL"].append(0.25*(volume/50)) 
          enzyme_dict[str(volume) + " uL Cost"].append(0.25*(volume/50)*enzyme_cost)
        enzyme_df = pd.DataFrame(enzyme_dict)
        enzyme_df.index = ['Enzyme amount']

        #iProof Primer Annealing Temperature
        annealing_temp_1 = mt.Tm_NN(Seq(self.fp))
        annealing_temp_2 = mt.Tm_NN(Seq(self.rp))
        if len(self.fp) or len(self.fp) > 20: 
          if annealing_temp_1 > annealing_temp_2: 
            annealing_temp_final = round(annealing_temp_2 + 3, 2)
          elif annealing_temp_1 < annealing_temp_2: 
            annealing_temp_final = round(annealing_temp_1 + 3, 2)

        elif len(self.fp) or len(self.rp) <=20: 
          if annealing_temp_1 > annealing_temp_2: 
            annealing_temp_final = round(annealing_temp_2,2)
          elif annealing_temp_1 < annealing_temp_2: 
            annealing_temp_final = round(annealing_temp_1,2)
        annealing_time_seconds = 60 
        annealing_time_minutes = annealing_time_seconds/60 

        # Time for each Processes in seconds
        if self.countDNAbase()[1] < 65:
          initial_Denaturation = 30 
        else:
          initial_Denaturation = 4*60

        if annealing_temp_1 > 65 and annealing_temp_2 > 65:
          cycles = 30
        else:
          cycles = 35

        Denaturation = 30  
        final_extension = 5 * 60 

        # Total Time calculation
        total_pcr_time_minutes = round((initial_Denaturation + (Denaturation + annealing_time_seconds + extension_time_seconds)*cycles + final_extension)/60, 2)
        total_pcr_time_hours = round(total_pcr_time_minutes/60, 2)
        
        factor = ['Annealing Temperature:','Annealing Time:','Extention Time:','Total PCR reaction time is:' ]
        stats_data = {'Result': ["{} degrees Celcius".format(annealing_temp_final), "{} seconds or {} minutes".format(annealing_time_seconds, annealing_time_minutes), "{} seconds or {} minutes".format(extension_time_seconds, extension_time_minutes), "{} minutes or {} hours".format(total_pcr_time_minutes, total_pcr_time_hours)]}
        stats = pd.DataFrame(data = stats_data, index = factor)

        return enzyme_df, stats 
    
