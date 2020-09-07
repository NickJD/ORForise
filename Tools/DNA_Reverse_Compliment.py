import collections

def revCompIterative(watson): #Gets Reverse Complement
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                   'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev:
        try:
            crick += complements[nt]
        except KeyError:
            crick += nt # Do not modify error
    return crick

def sortORFs(tool_ORFs):
    tool_ORFs_tmp = collections.OrderedDict()
    tool_ORFs_Sorted = collections.OrderedDict()
    for cds,data in tool_ORFs.items(): # TransDecoder CDS must be ordered.
        cds = cds.split(',')
        data.insert(0,cds[1])
        tool_ORFs_tmp.update({int(cds[0]):(data)})
    tool_ORFs_tmp = sorted(tool_ORFs_tmp.items())
    for orf in tool_ORFs_tmp:
        key = str(orf[0])
        data = orf[1]
        pos = key+','+data[0]
        value = [data[1],data[2],data[3]]
        tool_ORFs_Sorted.update({pos:value})
    return tool_ORFs_Sorted