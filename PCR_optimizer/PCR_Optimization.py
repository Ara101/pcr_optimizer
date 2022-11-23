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
