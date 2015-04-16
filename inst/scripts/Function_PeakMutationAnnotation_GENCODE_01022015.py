###Add parameter: output nearest one or two from two directions
###Add hg38
###Add internal promoter promoter_radius2
###V2: USE GENCODEhg19 or mm10, mm9 as reference
###Add stop at intron(adjacent)
###only nearest one if peak is not located in bidirectional region, instead of nearest left and right
from bisect import *
import ast, math, sys, string, math, shutil, math, os, gzip, time, glob, multiprocessing
from datetime import datetime
from shutil import rmtree


def FindPeakMutation(inputfile,outputfile,search_radius,promoter_radius,promoter_radius2,genome,adjacent,pwd,SNP,PromoterStop,NearestTwoDirection):
   print 'python process start:', datetime.now()
   print 'Load Reference'
   # Formatted output.
   def out_string(peak, M, record, comment):
      if record[6] == '+':
         distance = a - M
      else:
         distance = b - M
      return out_string_v2(peak, distance, record, comment)

   def out_string_exon(peak, M, record, comment):
      '''
      if m<exona and n<exonb:
         distance = n-exona
      if m>=exona and n<exonb:
         distance = n-m
      if m>=exona and n>=exonb:
         distance = exonb-m
      if m<exona and n>=exonb:
         distance = exonb-exona
      '''
      distance = [exona,exonb]
      return out_string_v2(peak, distance, record, comment)
   
   def out_string_intron(peak, M, record, comment):
      '''
      if m<introna and n<intronb:
         distance = n-introna
      if m>=introna and n<intronb:
         distance = n-m
      if m>=introna and n>=intronb:
         distance = intronb-m
      if m<introna and n>=intronb:
         distance = intronb-introna
      '''
      distance = [introna,intronb]
      return out_string_v2(peak, distance, record, comment)
   def out_string_cds(peak, M, record, comment):
      distance = [cdsa,cdsb]
      return out_string_v2(peak, distance, record, comment)
   def out_string_utr(peak, M, record, comment):
      distance = [utra,utrb]
      return out_string_v2(peak, distance, record, comment)   
   
   def out_string_v2(peak, distance, record, comment):
      chromo, source, TSS, TTS, strand = record[0], record[1], record[3], record[4], record[6]
      gene_name = record[8].strip().split(';')[1].strip().split(' ')[1][1:-1]
      transcript_id = record[8].strip().split(';')[0].strip().split(' ')[1][1:-1]
      gnote = reduce(lambda x, y: x + '\t' + y, [chromo, TSS, TTS, strand, gene_name, source, transcript_id])
      return (peak + '\t' + str(n - m) + '\t' + str(distance) + '\t' + comment + '\t' + gnote)
   
   ######load exon, intron, geneid index information
   geneid_set =[]
   exonband =[]
   intronband =[]
   cdsband =[]
   utrband =[]
   for line in open(pwd + r'/GENCODE_' + genome + '_EXONinfo.txt','r'):
      line = ast.literal_eval(line.strip().split('\t')[1])
      exonband.append(line)
   
   for line in open(pwd + r'/GENCODE_' + genome + '_INTRONinfo.txt','r'):
      line = ast.literal_eval(line.strip().split('\t')[1])
      intronband.append(line)
      
   for line in open(pwd + r'/GENCODE_' + genome + '_CDSinfo.txt','r'):
      line = ast.literal_eval(line.strip().split('\t')[1])
      cdsband.append(line)

   for line in open(pwd + r'/GENCODE_' + genome + '_UTRinfo.txt','r'):
      line = ast.literal_eval(line.strip().split('\t')[1])
      utrband.append(line)
      
   for line in open(pwd + r'/GENCODE_' + genome + '_GENEID_index.txt','r'):
      line = ast.literal_eval(line.strip().split('\t')[0])
      geneid_set.append(line)

   #index gemeid_set
   genemap=dict()
   geneset=set()
   for i in xrange(len(geneid_set)):
      key=geneid_set[i][1]
      genemap[key] = geneid_set[i][0]
      geneset.add(key)
   


   print 'Check Reference files'

   '''
   PURPOSE
      Read and parse transcript file.
   PARAMETERS
      TRANSCRIPT       List of chromosomes; each chromosome itself is a list of transcript strand.
      TRANSCRIPTLeft   List of chromosomes; each chromosome itself is a list of lower index of transcript strand.
   '''
   TRANSCRIPT = []
   TRANSCRIPTLeft = []
   transcript_chrm = {}
   ichrm = newChrmId = -1
   for line in open(pwd + 'gencode.' + genome + '.annotation_GENE_GTF.txt', 'r'):
      line = line.strip().split('\t')
      if line[0][3:] in transcript_chrm:
         ichrm = transcript_chrm[ line[0][3:] ]
      else:
         newChrmId += 1
         ichrm = newChrmId
         transcript_chrm[ line[0][3:] ] = ichrm
         TRANSCRIPT.append([])
         TRANSCRIPTLeft.append([])
      TRANSCRIPT[ichrm].append( line )
      TRANSCRIPTLeft[ichrm].append( int(line[3]) )

   '''
   PURPOSE
      Sort transcript by strand end id.
   PARAMETERS
      TRANSCRIPTRightSort   List of chromosomes; each chromosome sorted by strand right location.
      TRANSCRIPTRight       List of chromosomes; each chromosome itself is a list of right location.
      TRANSCRIPTRightId     List of chromosomes; each chromosome itself is a list of exon id in the file.
   '''
   TRANSCRIPTRightSort = []
   TRANSCRIPTRight     = []
   TRANSCRIPTRightID   = []

   for transcript in TRANSCRIPT:

      TRANSCRIPTRightSort.append([])
      TRANSCRIPTRight.append([])
      TRANSCRIPTRightID.append([])

      for i in xrange(len(transcript)):
         TRANSCRIPTRightSort[-1].append( (int(transcript[i][4]), i) )

      TRANSCRIPTRightSort[-1].sort(key = lambda x: x[0])
      for key in TRANSCRIPTRightSort[-1]:
         TRANSCRIPTRight[-1].append(key[0])
         TRANSCRIPTRightID[-1].append(key[1])

   '''
   PURPOSE
      Prepare lists of begin id & key for all genes and for protein genes.
   PARAMETERS
      TSBEGIN      List of chromosome; each chromosome itself a list of (arrow begin index, line number in the transcript).
      TSBEGINKEY   List of chromosome; each chromosome itself a list of (arrow begin index).
      TSLINEID     List of chromosome; each chromosome itself a list of (line number in the transcript).

      PCBEGIN, PCBEGINKEY, PCLINEID       Same as above, but for protein codings.
   '''
   TSBEGIN    = []
   TSBEGINKEY = []
   TSLINEID   = []

   PCBEGIN = []
   PCBEGINKEY = []
   PCLINEID = []

   for transcript in TRANSCRIPT:

      TSBEGIN.append([])
      TSBEGINKEY.append([])
      TSLINEID.append([])

      PCBEGIN.append([])
      PCBEGINKEY.append([])
      PCLINEID.append([])

      for i in xrange(len(transcript)):
         record = transcript[i]

         # Append to gene list.
         if record[6] == '+':
            begin = int(record[3])
         else:
            begin = int(record[4])
         TSBEGIN[-1].append( (begin, i) )

         # Append to protein coding gene list.
         if record[1] == 'protein_coding':
            PCBEGIN[-1].append( (begin, i) )

      TSBEGIN[-1].sort(key = lambda x: x[0])
      for key in TSBEGIN[-1]:
         TSBEGINKEY[-1].append(key[0])
         TSLINEID[-1].append(key[1])

      PCBEGIN[-1].sort(key = lambda x: x[0])
      for key in PCBEGIN[-1]:
         PCBEGINKEY[-1].append(key[0])
         PCLINEID[-1].append(key[1])

   '''
   PURPOSE
      Parse the peak file.
   ALGORITHM
      (1) Identify exon;
      (2) Identify intron; CDS;UTR
      (3) Find nearest protein coding;
      (4) Find neighbors within a range.
   '''
   print 'fixed reference done:', datetime.now()

   fout = open(outputfile,'w')
   print 'Start Annotation'

   count = -1
   for line in open(inputfile, 'r'):

      # Neglect the comment line.
      count += 1
      if count == 0:
         sline =line.strip().split('\t')
         fout.write(sline[0] + '\t' +sline[1] + '\t' +sline[2] + '\t' +sline[3] + '\t' )
         fout.write('PeakLength' + '\t' + 'peakMtoStart_Overlap' + '\t' + 'type' + '\t' + 'BidirenctionalRegion' + '\t')
         fout.write('Chr' + '\t' + 'TSS'  + '\t' + 'TTS' + '\t' + 'strand' + '\t' + 'gene_name'  + '\t' + 'source'+ '\t' + 'transID' + '\n')
         continue

      # ---------------------------------------
      # Parse the peak information (m < n).
      #
      #           middle
      #             |
      #     ---m---------n---
      #
      # ---------------------------------------
      line = line.strip().split()
      pkchrm = line[1]
      if pkchrm[0:3].upper() == 'CHR':
         pkchrm = pkchrm[3:]
      m, n = int(line[2]), int(line[3])
      middle = (m + n) / 2
      peakLeft  = m
      peakRight = n

      # The information about the peak to be printed.
      pkhd = reduce(lambda x, y: x + '\t' + y, line[0:4])

      # ---------------------------------------
      # Annotate exon;intron
      # ---------------------------------------

      # Check if the chromosome has been registered.
      if pkchrm not in transcript_chrm:
         print pkchrm, 'Chromosome not registered'
         continue
      else:
         ichrm = transcript_chrm[ pkchrm ]
         transcript = TRANSCRIPT[ichrm]
         transcriptLeft    = TRANSCRIPTLeft[ichrm]
         transcriptRight   = TRANSCRIPTRight[ichrm]
         transcriptRightID = TRANSCRIPTRightID[ichrm]
         
      # A set holding everything that has been marked.
      myNeighbor = set()
      
      # Find the search range.
      iMin = bisect_left(transcriptRight, peakLeft)

      setRight = set()
      for i in xrange(iMin, len(transcriptRight)):
         setRight.add(transcriptRightID[i]) 

      iMax = bisect_right(transcriptLeft, peakRight, lo=iMin+1)
      setLeft = set(range(iMax))

      # Search the range.
      found_exon_protein = False
      found_intron_protein = False

      #search exon region
      for transcriptID in setRight.intersection(setLeft):

         record = transcript[transcriptID]
         a, b = int(record[3]), int(record[4])

         if n < a or m > b:
            continue
         else: 
            geneid_B=transcript[transcriptID][8].strip().split(';')[0][9:-1]            
            if len(exonband[genemap[geneid_B]])>0:
               for i in xrange(len(exonband[genemap[geneid_B]])):
                  exona=exonband[genemap[geneid_B]][i][0]
                  exonb=exonband[genemap[geneid_B]][i][1]
                  if n < exona or m > exonb:
                     pass
                  else:
                     fout.write(out_string_exon(pkhd, middle, record, 'Exon\tNA')+ '\n')

                     #Check when exon finds, make stop label for SNP and adjacent=True
                     if record[1] == 'protein_coding':
                        if adjacent == True or SNP == True:
                           found_exon_protein = True
         myNeighbor.add(transcriptID)
      #search intron region
      if found_exon_protein == False:
         for transcriptID in setRight.intersection(setLeft):
            record = transcript[transcriptID]
            a, b = int(record[3]), int(record[4])
            if n < a or m > b:
               continue
            else:
               geneid_B=transcript[transcriptID][8].strip().split(';')[0][9:-1]
               if len(intronband[genemap[geneid_B]])>0:
                  for i in xrange(len(intronband[genemap[geneid_B]])):
                     introna=intronband[genemap[geneid_B]][i][0]
                     intronb=intronband[genemap[geneid_B]][i][1]
                     if n < introna or m > intronb:
                        pass
                     else:
                        fout.write(out_string_intron(pkhd, middle, record, 'Intron\tNA')+ '\n')

                        #Check when intron finds, make stop label for adjacent=True
                        if record[1] == 'protein_coding':
                           if adjacent == True:
                              found_intron_protein = True                  
            myNeighbor.add(transcriptID)
      #search cds region
      if found_exon_protein == False and found_intron_protein == False:
         for transcriptID in setRight.intersection(setLeft):
            record = transcript[transcriptID]
            a, b = int(record[3]), int(record[4])
            if n < a or m > b:
               continue
            else:
               geneid_B=transcript[transcriptID][8].strip().split(';')[0][9:-1]
               if len(cdsband[genemap[geneid_B]])>0:
                  for i in xrange(len(cdsband[genemap[geneid_B]])):
                     cdsa=cdsband[genemap[geneid_B]][i][0]
                     cdsb=cdsband[genemap[geneid_B]][i][1]
                     if n < cdsa or m > cdsb:
                        pass
                     else:
                        fout.write(out_string_cds(pkhd, middle, record, 'cds\tNA')+ '\n')
            myNeighbor.add(transcriptID)            
      #search utr region
      if found_exon_protein == False and found_intron_protein == False:
         for transcriptID in setRight.intersection(setLeft):
            record = transcript[transcriptID]
            a, b = int(record[3]), int(record[4])
            if n < a or m > b:
               continue
            else:
               geneid_B=transcript[transcriptID][8].strip().split(';')[0][9:-1]
               if len(utrband[genemap[geneid_B]])>0:
                  for i in xrange(len(utrband[genemap[geneid_B]])):
                     utra=utrband[genemap[geneid_B]][i][0]
                     utrb=utrband[genemap[geneid_B]][i][1]
                     if n < utra or m > utrb:
                        pass
                     else:
                        fout.write(out_string_utr(pkhd, middle, record, 'utr\tNA')+ '\n')                                       
            myNeighbor.add(transcriptID)
#############################
      #search internal-promoter region
      if found_exon_protein == False and found_intron_protein == False:
         for transcriptID in setRight.intersection(setLeft):
            record = transcript[transcriptID]
            a, b = int(record[3]), int(record[4])
            if n < a or m > b:
               continue
            else:
               if record[6] == '+':
                  if m < a + promoter_radius2:
                     fout.write(out_string(pkhd, middle, record, 'Promoter_internal\tNA')+ '\n')
                     myNeighbor.add(transcriptID)
               elif record[6] == '-':
                  if n > b - promoter_radius2:
                     fout.write(out_string(pkhd, middle, record, 'Promoter_internal\tNA')+ '\n')
                     myNeighbor.add(transcriptID)                                                  
#####################################            



      # If found protein_coding. Annotate the next peak.
      if found_exon_protein or found_intron_protein:
         continue
        
      # --------------------------------------------
      # Search promoter in protein coding genes.
      # --------------------------------------------
      pcbeginkey = PCBEGINKEY[ichrm]
      pclineid   = PCLINEID[ichrm]

      pcid = bisect(pcbeginkey, middle)

      # Find the left promoter.
      i = pcid - 1
      found_left_promoter, left_d = False, 0
      left_promoter = []
      while i >= 0 and not found_left_promoter and left_d < promoter_radius:
         record = transcript[ pclineid[i] ]
         a, b = int(record[3]), int(record[4])
         left_d = m - b

         if left_d < promoter_radius:
            if record[6] == '-':
               found_left_promoter = True
               left_promoter.append( record )
            i = i - 1

      # Print the non-nearest left promoters.
      for record in left_promoter[1:]:
         a, b = int(record[3]), int(record[4])
         #print out_string(pkhd, middle, record, 'Promotor_L\tN')
         fout.write(out_string(pkhd, middle, record, 'Promoter_L\tN')+ '\n')
         
      # Find the right promoter.
      i = pcid
      found_right_promoter, right_d = False, 0
      right_promoter = []
      while i < len(pcbeginkey) and not found_right_promoter and right_d < promoter_radius:
         record = transcript[ pclineid[i] ]
         a, b = int(record[3]), int(record[4])
         right_d = a - n

         if right_d < promoter_radius:
            if record[6] == '+':
               found_right_promoter = True
               right_promoter.append( record )
            i = i + 1
      # Print the non-nearest right promoters.
      for record in right_promoter[1:]:
         a, b = int(record[3]), int(record[4])
         #print out_string(pkhd, middle, record, 'Promotor_R\tN')
         fout.write(out_string(pkhd, middle, record, 'Promoter_R\tN') + '\n')

      # If promoters are found on both sides.
      if found_left_promoter and found_right_promoter:
         record = left_promoter[0]
         a, b = int(record[3]), int(record[4])
         #print out_string(pkhd, middle, record, 'Promotor_L\tY')
         fout.write(out_string(pkhd, middle, record, 'Promoter_L\tY') + '\n')

         record = right_promoter[0]
         a, b = int(record[3]), int(record[4])
         #print out_string(pkhd, middle, record, 'Promotor_R\tY')
         fout.write(out_string(pkhd, middle, record, 'Promoter_R\tY')+ '\n')

        
      # Find the right nearest neighbor if no right promoter found.
      if found_left_promoter and not found_right_promoter:
         is_right_bidirectional = False
         i = pcid
         found_right = False
         while i < len(pcbeginkey) and not found_right:
            record = transcript[ pclineid[i] ]
            a, b = int(record[3]), int(record[4])
            if min(a, b) > n:
               found_right = True
               if record[6] == '+':
                  is_right_bidirectional = True
                  #print out_string(pkhd, middle, record, 'Nearest_R\tY')
                  fout.write(out_string(pkhd, middle, record, 'Nearest_R\tY') + '\n')
               else:
                  #print out_string(pkhd, middle, record, 'Nearest_R\tN')
                  fout.write(out_string(pkhd, middle, record, 'Nearest_R\tN') + '\n')
            else:
               i = i + 1
               
         record = left_promoter[0]
         a, b = int(record[3]), int(record[4])
         if is_right_bidirectional:
            #print out_string(pkhd, middle, record, 'Promotor_L\tY')
            fout.write(out_string(pkhd, middle, record, 'Promoter_L\tY') + '\n')
         else:
            #print out_string(pkhd, middle, record, 'Promotor_L\tN')
            fout.write(out_string(pkhd, middle, record, 'Promoter_L\tN') + '\n')
         
      # Find the left nearest neighbor if no left promoter found.
      if not found_left_promoter and found_right_promoter:
         is_left_bidirectional = False
         i = pcid - 1
         found_left = False
         while i >= 0 and not found_left:
            record = transcript[ pclineid[i] ]
            a, b = int(record[3]), int(record[4])
            if max(a, b) < m:
               found_left = True
               if record[6] == '-':
                  is_left_bidirectional = True
                  #print out_string(pkhd, middle, record, 'Nearest_L\tY')
                  fout.write(out_string(pkhd, middle, record, 'Nearest_L\tY') + '\n')
               else:
                  #print out_string(pkhd, middle, record, 'Nearest_L\tN')
                  fout.write(out_string(pkhd, middle, record, 'Nearest_L\tN') + '\n')
            else:
               i = i - 1

         record = right_promoter[0]
         a, b = int(record[3]), int(record[4])
         if is_left_bidirectional:
            #print out_string(pkhd, middle, record, 'Promotor_R\tY')
            fout.write(out_string(pkhd, middle, record, 'Promoter_R\tY') + '\n')
         else:
            #print out_string(pkhd, middle, record, 'Promotor_R\tN')
            fout.write(out_string(pkhd, middle, record, 'Promoter_R\tN') + '\n')

      # Stop here if any promoter is found and if PromoterStop index is True, else if no promoter is found/or PromoterStop index is False, contitue for further
      # search in search radius.
      if PromoterStop == False:
         found_left_promoter = PromoterStop
         found_right_promoter = PromoterStop
      elif PromoterStop == True:
         pass
         
      if found_left_promoter or found_right_promoter:
         continue

      # ----------------------------------------------------
      # Search nearest neighbor in protein coding genes.
      # ----------------------------------------------------

      # Find the left nearest transcript.
      i = pcid - 1
      found_left = False
      while i >= 0 and not found_left:
         lineL = pclineid[i]
         record = transcript[lineL]
         a, b = int(record[3]), int(record[4])
         if max(a, b) < m:
            found_left = True
         else:
            i = i - 1

      # Find the right nearest transcript.
      i = pcid
      found_right = False
      while i < len(pcbeginkey) and not found_right:
         lineR = pclineid[i]
         record = transcript[lineR]
         a, b = int(record[3]), int(record[4])
         if min(a, b) > n:
            found_right = True
         else:
            i = i + 1

      # Check if is bidirectional.
      if found_left and found_right:
         recordL = transcript[lineL]
         recordR = transcript[lineR]

         if recordL[6] == '-' and recordR[6] == '+':

            a, b = int(recordL[3]), int(recordL[4])
            myNeighbor.add(lineL)
            #print out_string(pkhd, middle, recordL, 'Nearest_L\tY')
            fout.write(out_string(pkhd, middle, recordL, 'Nearest_L\tY') + '\n')

            a, b = int(recordR[3]), int(recordR[4])
            myNeighbor.add(lineR)
            #print out_string(pkhd, middle, recordR, 'Nearest_R\tY')
            fout.write(out_string(pkhd, middle, recordR, 'Nearest_R\tY') + '\n')

         else:

            a, b = int(recordL[3]), int(recordL[4])
            if recordL[6] == '+':
               dL = a - middle 
            else:
               dL = b - middle

            a, b = int(recordR[3]), int(recordR[4])
            if recordR[6] == '+':
               dR = a - middle 
            else:
               dR = b - middle

            if NearestTwoDirection == False:
               if abs(dL) < abs(dR):
                  #print out_string_v2(pkhd, dL, recordL, 'Nearest_L\tN')
                  fout.write(out_string_v2(pkhd, dL, recordL, 'Nearest\tN') + '\n')
                  myNeighbor.add(lineL)
               else:
                  #print out_string_v2(pkhd, dR, recordR, 'Nearest_R\tN')
                  fout.write(out_string_v2(pkhd, dR, recordR, 'Nearest\tN') + '\n')
                  myNeighbor.add(lineR)

            if NearestTwoDirection == True:
               fout.write(out_string_v2(pkhd, dL, recordL, 'Nearest_L\tN') + '\n')
               fout.write(out_string_v2(pkhd, dR, recordR, 'Nearest_R\tN') + '\n')
               myNeighbor.add(lineL)
               myNeighbor.add(lineR)

      elif found_left:

         myNeighbor.add(lineL)
         record = transcript[lineL]
         a, b = int(record[3]), int(record[4])
         myNeighbor.add(lineL)
         #print out_string(pkhd, middle, record, 'Nearest_L\tN')
         fout.write(out_string(pkhd, middle, record, 'Nearest\tN') + '\n')

      elif found_right:

         myNeighbor.add(lineR)
         record = transcript[lineR]
         a, b = int(record[3]), int(record[4])
         myNeighbor.add(lineR)
         #print out_string(pkhd, middle, record, 'Nearest_R\tN')
         fout.write(out_string(pkhd, middle, record, 'Nearest\tN')+ '\n')


      # -------------------------------------------
      # Print everything within searching radius.
      # -------------------------------------------
      tsbeginkey  = TSBEGINKEY[ichrm]
      tslineid    = TSLINEID[ichrm]

      lower_bound = max(0, middle - search_radius)
      upper_bound = min(tsbeginkey[-1], middle + search_radius)

      lower_id = bisect(tsbeginkey, lower_bound)
      upper_id = bisect(tsbeginkey, upper_bound, lo = lower_id)

      for key_id in xrange(lower_id, upper_id):
         line_id = tslineid[key_id]
         if line_id not in myNeighbor:
            distance = middle - tsbeginkey[key_id]
            record = transcript[line_id]
            #print out_string_v2(pkhd, distance, record, 'Neighbor\tN')
            fout.write(out_string_v2(pkhd, distance, record, 'Neighbor\tN') + '\n')

   fout.close()

   print 'Finish Annotation'
   print 'python process end:', datetime.now()



