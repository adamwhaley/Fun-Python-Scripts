#usr/bin
import sys
import os

def patterns(default=1):
	"""Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """

	question="Please select a comparison:\n 1:FL-FH   2:FL-ML   3:ML-MH   4:FH-MH\n"

	prompt=" [1,2,3,4] "
	valid = {"1":['402_404_diff','406_408_diff','412_410_diff','414_416_diff'],   
			"2":['402_403_diff','406_407_diff','412_413_diff','414_415_diff'],
			"3":['403_405_diff','407_409_diff','413_411_diff','415_417_diff'],
            "4":['FH-MH1-test','FH2-MH2-diff','FH3-MH3-diff','FH-MH4n']}

	
	sys.stdout.write(question + prompt)
	choice = raw_input().lower()
	if default is not None and choice == '':
		return valid[default]
	elif choice in valid:
		### Open output file for specified choice
		if choice=='1':
			fout = open('FL-FH_timing_analysis.txt','w')

		elif choice=='2':
			fout = open('FL-ML_timing_analysis.txt','w')
		elif choice=='3':
			fout = open('ML-MH_timing_analysis.txt','w')
		elif choice=='4':
			fout = open('FH-MH_timing_analysis.txt','w')
			fout1 = open('FH-MH_p1.txt','w')
			fout2 = open('FH-MH_p2.txt','w')
			fout3 = open('FH-MH_p3.txt','w')
			fout4 = open('FH-MH_p4.txt','w')
			fout5 = open('FH-MH_p5.txt','w')
			fout6 = open('FH-MH_p6.txt','w')

		fout.write("#Pattern codes:\n# 1: Starts with 0 expression in both samples, faster response in sample 2\n")
		fout.write("# 2: Starts with 0 expression in both samples, faster response in sample 1\n")
		fout.write("# 3: Starts with expression in both samples, sample 1 turns off while sample 2 remains expressed\n")
		fout.write("# 4: Starts with expression in both samples, sample 2 turns off while sample 1 remains expressed\n")
		fout.write("# 5: Sample 2 has no expression, high expression in sample 1.  Sample 2 begins to be expressed\n")
		fout.write("# 6: Sample 1 has no expression, high expression in sample 2.  Sample 1 begins to be expressed\n")

		### Open gene_exp.diff in directories for chosen comparison
		genes = {}
		dir = valid[choice]
		for file in dir:
			fh = open(file + '/gene_exp.diff','r')
			for line in fh.readlines():
				if line.startswith('test_id'):
					continue
				### Read line, looking for 'yes' to denote significant differential expression
				### as well as fold change value relative to 0
				line=line.rstrip('\n').split('\t')
				name = line[2]
				if line[13]=='yes':
					
					if float(line[7])==0 and float(line[8])!=0:
						fc = 'P'
						#return line
					elif float(line[7])!=0 and float(line[8])==0:
						fc = 'N'
						#return line
					

					else:
						if float(line[9])>0:
							fc = 'p'
						else:
							fc = 'n'
					### Add relevant results to dictionary
					if genes.has_key(name):
						genes[name].append((file,fc))
					else:
						genes[name] = [(file,fc)]
				elif float(line[7])>0 and float(line[8])>0:
					if genes.has_key(name):
						genes[name].append((file,'B'))
					else:
						genes[name] = [(file,'B')]
				else:
					if genes.has_key(name):
						genes[name].append((file,"U"))
					else:
						genes[name] = [(file,"U")]
						
			fh.close()
		return_lst = []
		for gene in genes.keys():
			#eliminate entries with a common pattern through all samples in series
			if genes[gene][0][1]==genes[gene][1][1] and genes[gene][0][1]==genes[gene][2][1] and genes[gene][0][1]==genes[gene][3][1]:
				continue
			else:
				lst = genes[gene]
				pattern = lst[0][1]+lst[1][1]+lst[2][1]+lst[3][1]
				return_lst.append(pattern)
				
				if pattern[0]=='U':
					valid1 = ["UUUP","UUPP","UPPP","UPPp","UPpp","UPpB","UPBB","UppB","UpBB"]
					valid2 = ["UUUN","UUNN","UNNN","UNNn","UNnn","UNnB","UNBB","UnnB","UnBB"]
					if pattern in valid1:
						fout.write(gene + "\t1\n")
						if ',' in gene:
							gene = gene.split(',')
							for g in gene:
								fout1.write(g + '\n')
						else:
							fout1.write(gene + '\n')

					elif pattern in valid2:
						fout.write(gene + "\t2\n")
						if ',' in gene:
							gene = gene.split(',')
							for g in gene:
								fout2.write(g + '\n')
						else:
							fout2.write(gene + '\n')

				elif pattern[0]=='B':
					valid = ["BBBp","BBpp","Bppp","BppP","BpPP","BpPU","BpUU","BPPU","BPUU"]
					valid2 = ["BBBn","BBnn","Bnnn","BnnN","BnNN","BnNU","BnUU","BNNU","BNUU"]
					if pattern in valid:
						fout.write(gene + "\t3\n")
						if ',' in gene:
							gene = gene.split(',')
							for g in gene:
								fout3.write(g + '\n')
						else:
							fout3.write(gene + '\n')

					elif pattern in valid2:
						fout.write(gene + "\t4\n")
						if ',' in gene:
							gene = gene.split(',')
							for g in gene:
								fout4.write(g + '\n')
						else:
							fout4.write(gene + '\n')

				elif pattern[0]=='N':
					valid = ["NNNn","NNnn","Nnnn","NnnB","NnBB","NBBB"]
					if pattern in valid:
						fout.write(gene + "\t5\n")
						if ',' in gene:
							gene = gene.split(',')
							for g in gene:
								fout5.write(g + '\n')
						else:
							fout5.write(gene + '\n')

				elif pattern[0]=='P':
					##print pattern
					valid2 = ["PPPp","PPpp","Pppp","PppB","PpBB","PBBB"]
					if pattern in valid2:
						fout.write(gene + "\t6\n")
						if ',' in gene:
							gene = gene.split(',')
							for g in gene:
								fout6.write(g + '\n')
						else:
							fout6.write(gene + '\n')
						



	else:
		sys.stdout.write("Not valid int! Try again\n\n")
