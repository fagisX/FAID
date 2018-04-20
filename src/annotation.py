import pandas as pd

def int_to_roman(input):
   """
   Convert an integer to Roman numerals.

   Examples:
   >>> int_to_roman(0)
   Traceback (most recent call last):
   ValueError: Argument must be between 1 and 3999

   >>> int_to_roman(-1)
   Traceback (most recent call last):
   ValueError: Argument must be between 1 and 3999

   >>> int_to_roman(1.5)
   Traceback (most recent call last):
   TypeError: expected integer, got <type 'float'>

   >>> for i in range(1, 21): print int_to_roman(i)
   ...
   I
   II
   III
   IV
   V
   VI
   VII
   VIII
   IX
   X
   XI
   XII
   XIII
   XIV
   XV
   XVI
   XVII
   XVIII
   XIX
   XX
   >>> print int_to_roman(2000)
   MM
   >>> print int_to_roman(1999)
   MCMXCIX
   """
   if type(input) != type(1):
      raise TypeError, "expected integer, got %s" % type(input)
   if not 0 < input < 4000:
      raise ValueError, "Argument must be between 1 and 3999"
   ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
   nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
   result = ""
   for i in range(len(ints)):
      count = int(input / ints[i])
      result += nums[i] * count
      input -= ints[i] * count
   return result

def roman_to_int(input):
    """
    Convert a roman numeral to an integer.

    >>> r = range(1, 4000)
    >>> nums = [int_to_roman(i) for i in r]
    >>> ints = [roman_to_int(n) for n in nums]
    >>> print r == ints
    1

    >>> roman_to_int('VVVIV')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: VVVIV
    >>> roman_to_int(1)
    Traceback (most recent call last):
     ...
    TypeError: expected string, got <type 'int'>
    >>> roman_to_int('a')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: A
    >>> roman_to_int('IL')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: IL
    """
if __name__ == "__main__":
    if type(input) != type(""):
        raise TypeError, "expected string, got %s" % type(input)
    input = input.upper()
    nums = ['M', 'D', 'C', 'L', 'X', 'V', 'I']
    ints = [1000, 500, 100, 50, 10, 5, 1]
    places = []
    for c in input:
        if not c in nums:
            raise ValueError, "input is not a valid roman numeral: %s" % input
    for i in range(len(input)):
        c = input[i]
        value = ints[nums.index(c)]
        # If the next place holds a larger number, this value is negative.
        try:
            nextvalue = ints[nums.index(input[i + 1])]
            if nextvalue > value:
                value *= -1
        except IndexError:
            # there is no next place.
            pass
        places.append(value)
    sum = 0
    for n in places: sum += n
    # Easiest test for validity...
    if int_to_roman(sum) == input:
        return sum
    else:
        raise ValueError, 'input is not a valid roman numeral: %s' % input

   gtf = pd.read_csv('./ref_data/sacCer3.txt', sep='\t')
   genes = set(gtf.name.str.strip().str.upper().unique())
   print len(genes)
   GPL = pd.read_csv('./expression_data/GPL884.txt', sep='\t', header=10)

   not_in = 0
   for name in list(GPL.NAME.values):
     if name == 'BrightCorner' or name == 'NegativeControl':
         continue
     if name not in genes:
         # print name
         not_in +=1

   print not_in
   #
   GPL = pd.read_csv('./expression_data/GPL2529-9333.txt', sep='\t', header=None)
   GPL = GPL[GPL.iloc[: , 2] == 'Saccharomyces cerevisiae']
   not_in = 0
   for name in list(GPL.iloc[:, 1].values):
     if name == 'BrightCorner' or name == 'NegativeControl':
         continue
     if name not in genes:
         print name
         not_in +=1

   print not_in


   not_in = 0
   for name in genes:
     if name not in GPL.iloc[:, 1].unique():
         print name
         not_in += 1
   print not_in

   essential_genes = './ref_data/yeast_essential_genes.xlsx'
   eg = pd.read_excel(essential_genes)
   eg['ORF_name'] = eg['ORF_name'].str.strip().str.upper()


   insertion_gene_df = pd.read_csv('./expression_result/exp_insertion_total_for_exp_trascribed_region.xls', sep='\t')

   gtf_1kb = gtf[gtf.name.isin(insertion_gene_df['1kb'].unique())]
   gtf_hotspot = gtf[gtf.name.isin(insertion_gene_df['hotspot'].unique())]
   gtf_hotspot_1kb = gtf[gtf.name.isin(insertion_gene_df['hotspot_1kb'].unique())]
   gtf_total = gtf[gtf.name.isin(insertion_gene_df['total'].unique())]

   not_in = 0
   for protein in insertion_gene_df['total'].unique():
     if protein not in eg['ORF_name'].unique():
         not_in +=1
   print not_in, gtf_total.shape, eg['ORF_name'].unique().shape[0]-46

   not_in = 0
   for protein in gtf_1kb.name.unique():
     if protein not in eg['ORF_name'].unique():
         not_in +=1
   print not_in, gtf_1kb.shape
   #
   not_in = 0
   for protein in gtf_hotspot.name.unique():
     if protein not in eg['ORF_name'].unique():
         not_in +=1
   print not_in, gtf_hotspot.shape

   not_in = 0
   for protein in gtf_hotspot_1kb.name.unique():
     if protein not in eg['ORF_name'].unique():
         not_in +=1
   print not_in, gtf_hotspot_1kb.shape

   gtf = pd.read_csv('./ref_data/sacCer3.txt', sep='\t')
   ARSes = pd.read_csv('./ref_data/Confirmed_OriDB.xls', sep='\t')
   # print ARSes
   ARSes['start'] = ARSes.start - 1000
   ARSes['end'] = ARSes.end + 1000

   genes_in_1kb = []

   for i in range(gtf.shape[0]):
     name, chr, start, end = gtf.ix[i, 'name'], gtf.ix[i, 'chrom'], gtf.ix[i, 'txStart'], gtf.ix[i, 'txEnd']
     chr = roman_to_int(chr.replace('chr', ''))

     cur_ARSes = ARSes[ARSes.chr == chr]
     cur_ARSes = cur_ARSes[((cur_ARSes.start <= start)&(cur_ARSes.end>=start)) | \
                           ((cur_ARSes.start <= end)&(cur_ARSes.end>=end)) | \
                           ((cur_ARSes.start >= start)&(cur_ARSes.end<=end))]
     if cur_ARSes.shape[0] > 0:
         genes_in_1kb.append(name)
   f = open('genes_in_1kb.txt', 'w')
   for g in genes_in_1kb:
     f.write(g+'\n')
   f.close()


