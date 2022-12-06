# PCR_protocol_optimizer
A Python package to aid researchers in optimizing pcr protocols

Polymerase chain reaction (PCR) is a common molecular biology technique used to amplify DNA. Each enzyme has different optimal reaction conditions such as extension times and primer annealing temperatures depending on the properties of the target gene, the primers used to amplify it, and buffers present in the reaction. PCR enzymes are expensive reagents, so it is also important to consider cost before choosing a protocol. Additionally, PCRs can be ineffective due to user design errors. Primers must anneal at the correct sites on a gene and have a carefully calculated annealing temperature to ensure proper binding conditions. Checking for any errors before conducting PCR is also crucial to avoid wasting time and money due to a failed reaction. Considering all these variables can make designing an effective PCR experiment difficult. 

This tool optimizes the PCR protocol decision making process. First the gene and primer sequences are checked for errors and evaluated to ensure the primers are binding in the correct location on the target sequence. Next, enzyme amount, cost per reaction, annealing temperature, annealing time, extension time, and total PCR reaction time for two common PCR enzymes (iProof High-Fidelity Polymerase and Taq Polymerase) are reported to the user. The user can also designate a factor (time or cost) to focus their optimization. 

## Setup 

```pip install pcr_protocol_optimizer```

## Defining PCR object 

```test = PCR(gene, forward_primer, reverse_primer)```

