# Authored Lily Torp, K. Lionel Tukei 
# Fall 2022

class pcr(object): 
  
  def __init__ (self, gene, forward_primer, reverse_primer, template_type = "plasmid", enzymes = [], factor = "cost"): 
    self.sequence = gene
    self.length = len(gene)
    self.factor = factor
    self.template_type = template_type
    self.fp = forward_primer
    self.rp = reverse_primer

  def countDNAbase(gene): 
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
    return "The GC Content is",gc_content, "%"

  def checkPrimerGeneCompatability(gene, forward_primer, reverse_primer, startr = 0, stopr = 0 , startf = 0, stopf = 0): 
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

  def meltingTemperaturePhusion(seq): 
    """
    checks primer melting temperature (Tm)
    Adapted from Sigma Aldrich protocol: https://www.sigmaaldrich.com/US/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp
    Santa Lucia calculations
    Tm = temp where 50% of DNA is melted, and 50% of DNA is double stranded

    Parameters: 
    ------------
    seq : Primer sequnce as a string

    Returns: 
    ------------
    tm : int
      DNA melting temperature in degrees celsius

    """
    import math as math 
    H = 0 
    S = 0 

    for x in range(0,len(seq)+1): 
      check = str(seq[x:x+2]).lower()
  
      if check == "aa" or check  == "tt": 
        H += -9.1
        S += -0.0240

      elif check == "at": 
        H += -8.6
        S += -0.0239

      elif check == "ta": 
        H += -6.0
        S += -0.0169

      elif check == "ca" or check == "tg": 
        H += -5.8
        S += -0.0129

      elif check == "gt" or check == "ac": 
        H += -6.5 
        S += -0.0173

      elif check == "ct" or check == "ag": 
        H += -7.8 
        S += -0.0208

      elif check == "ga" or check == "tc": 
        H += -5.6 
        S += -0.0135

      elif check == "cg": 
        H += -11.9
        S += -0.0278

      elif check == "gc": 
        H += -11.1
        S += -0.0267

      elif check == "gg" or check == "cc": 
        H += -11.0
        S += -0.0266

    tm = (H/(-0.0108+S+0.00199*math.log(0.0000005/4)))-273.15+16.6*math.log10(0.05)
    return tm
