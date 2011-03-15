""" 
This library offers functions to parse different types of SNPs data from multiple formats into a snpsData objects.

Bjarni Vilhjalmsson, bvilhjal@usc.edu
"""
import time, sys, random
import os
import dbutils
import cPickle
import pdb
from snpsdata import *
from env import *

# this should be fixed


#Standard missing value is NA
missing_val = 'NA'
#Standard nt decoding, is using the IUPAC alphabet
nt_decoder = {'A':'A',
	      'C':'C',
	      'G':'G',
	      'T':'T',
	      'AG':'R',
	      'AC':'M',
	      'GT':'K',
	      'CT':'Y',
	      'AT':'W',
	      'CG':'S',
	      'Y':'Y',
	      'R':'R',
	      'W':'W',
	      'S':'S',
	      'K':'K',
	      'M':'M',
	      'D':'D',
	      'H':'H',
	      'V':'V',
	      'B':'B',
	      'X':'N', #Unknown base(s)
	      'N':'N', #Unknown base(s)
	      '-':'-', #Indel 
	      '|':missing_val}	#05/12/08 yh. add '-':'-' (deletion) and '|':'NA' (untouched)



# A couple of useful dictionaries:
#accessionName2Id = {'Ms-0': 92, 'Mt-0': 79, 'HR-10': 34, 'Knox-18': 4, 'Spr1-6': 20, 'Knox-10': 3, 'Spr1-2': 19, 'Nok-3': 80, 'Wa-1': 81, 'RRS-10': 2, 'C24': 57, 'Wt-5': 74, 'Ra-0': 69, 'Gu-0': 54, 'Mz-0': 73, 'Tsu-1': 78, 'Fei-0': 82, 'Bur-0': 93, 'Omo2-3': 22, 'Pu2-23': 30, 'Rmx-A180': 6, 'Kondara': 88, 'Tamm-2': 41, 'CS22491': 58, 'Zdr-6': 26, 'Ren-11': 48, 'Zdr-1': 25, 'Ren-1': 47, 'Ler-1': 55, 'Fab-2': 13, 'Yo-0': 61, 'Wei-0': 59, 'Got-7': 45, 'Tamm-27': 42, 'Kas-2': 75, 'Ts-5': 85, 'Ts-1': 84, 'Pu2-7': 29, 'Mr-0': 77, 'Ei-2': 53, 'Mrk-0': 72, 'Lz-0': 52, 'Bil-7': 16, 'Bil-5': 15, 'Sq-8': 38, 'Fab-4': 14, 'Sq-1': 37, 'Omo2-1': 21, 'Var2-1': 17, 'Var2-6': 18, 'Shahdara': 89, 'Uod-7': 50, 'Uod-1': 49, 'Lov-5': 12, 'Lov-1': 11, 'Gy-0': 68, 'Col-0': 62, 'Kin-0': 91, 'NFA-8': 35, 'Nd-1': 56, 'Got-22': 46, 'Br-0': 65, 'HR-5': 33, 'Ull2-3': 24, 'Ull2-5': 23, 'Est-1': 66, 'CIBC-17': 40, 'Ct-1': 76, 'Cvi-0': 51, 'Oy-0': 95, 'LL-0': 87, 'Bor-4': 28, 'Bor-1': 27, 'Pna-10': 8, 'Pna-17': 7, 'Ga-0': 71, 'Bay-0': 70, 'Eden-2': 10, 'Eden-1': 9, 'Pro-0': 86, 'Kz-1': 43, 'RRS-7': 1, 'Kz-9': 44, 'Edi-0': 94, 'An-1': 63, 'CIBC-5': 39, 'Ws-0': 60, 'Ws-2': 96, 'Van-0': 64, 'Rmx-A02': 5, 'Se-0': 83, 'Lp2-2': 31, 'Lp2-6': 32, 'NFA-10': 36, 'Ag-0': 67, 'Sorbo': 90}

#accessions250To2010 = { "Ag0":"Ag-0", "An1":"An-1", "Bay0B":"Bay-0", "Bay0":"Bay-0", "Bil-5":"Bil-5", "Bil-7":"Bil-7", "Bil5":"Bil-5", "Bil7":"Bil-7", "Bor1":"Bor-1", "Bor4":"Bor-4", "Bor4B":"Bor-4", "Br0":"Br-0", "BroA":"Br-0", "BurOB":"Bur-0", "Bur0":"Bur-0", "BUR0A":"Bur-0", "C24A":"C24", "C24B":"C24", "CIBC17":"CIBC-17", "Cibc17":"CIBC-17", "CIBC5":"CIBC-5", "Cibc5":"CIBC-5", "CS2249":"CS22491", "Col0":"Col-0", "LY1":"Col-0", "Ct1":"Ct-1", "Cvi0":"Cvi-0", "CviOA":"Cvi-0", "Eden-1":"Eden-1", "Eden-2":"Eden-2", "Eden1":"Eden-1", "Eden2":"Eden-2", "Edi0":"Edi-0", "El2":"Ei-2", "Ei2":"Ei-2", "Est1B":"Est-1", "Est1":"Est-1", "EST1A":"Est-1", "Fab-2":"Fab-2", "Fab-4":"Fab-4", "Fab2":"Fab-2", "Fab4":"Fab-4", "Fei0":"Fei-0", "FeiOB":"Fei-0", "Ga0":"Ga-0", "Got22":"Got-22", "Got7":"Got-7", "Got7A":"Got-7", "Gu0":"Gu-0", "Gy0":"Gy-0", "HR10":"HR-10", "HR5":"HR-5", "Kas2":"Kas-2", "Kas1":"Kas-2", "Kin0":"Kin-0", "Knox10":"Knox-10", "Knox-10":"Knox-10", "Kno10":"Knox-10", "Knox-18":"Knox-18", "Knox18":"Knox-18", "Kno18":"Knox-18", "Kondara":"Kondara", "Kz1":"Kz-1", "Kz9":"Kz-9", "Kz_9":"Kz-9", "LL0":"LL-0", "Ler1A":"Ler-1", "Ler1":"Ler-1","LY3":"Ler-1", "Lov-1":"Lov-1", "Lov1":"Lov-1", "Lov_1":"Lov-1", "Lov-5":"Lov-5", "Lov5":"Lov-5", "Lov5B":"Lov-5", "LOV5A":"Lov-5", "Lp22":"Lp2-2", "LP22":"Lp2-2", "Lp26":"Lp2-6", "LP26":"Lp2-6", "Lz0":"Lz-0", "Mr0":"Mr-0", "Mrk0":"Mrk-0", "Ms0":"Ms-0", "Ms_0":"Ms-0", "Mt0":"Mt-0", "Mz0":"Mz-0", "NFA10":"NFA-10", "Nfa10":"NFA-10", "NFA8":"NFA-8", "Nafa8B":"NFA-8", "Nd1":"Nd-1", "Nok3":"Nok-3", "Nok_3":"Nok-3", "Omo21":"Omo2-1", "Omo23":"Omo2-3", "Omo2_3":"Omo2-3", "Oy0":"Oy-0", "Pna-10":"Pna-10", "Pna10":"Pna-10", "Pna_10":"Pna-10", "Pna-17":"Pna-17", "Pna17":"Pna-17", "Pna_17":"Pna-17", "Pro0":"Pro-0", "Pu223":"Pu2-23", "Pu27":"Pu2-7", "Pu2_7":"Pu2-7", "RRS-10":"RRS-10", "RRS10":"RRS-10", "Rrs10B":"RRS-10", "RRS-7":"RRS-7", "RRS7":"RRS-7", "Rrs7B":"RRS-7", "Ra0":"Ra-0", "Ren1":"Ren-1", "Ren11":"Ren-11", "Ren_11":"Ren-11", "Rmx-A02":"Rmx-A02", "RmxA02":"Rmx-A02", "Rmx-A180":"Rmx-A180", "RmxA180":"Rmx-A180", "Se0":"Se-0", "Se0_2":"Se-0", "ShaA":"Shahdara", "Sorbo":"Sorbo", "Spr12":"Spr1-2", "Spr1_2":"Spr1-2", "Spr16":"Spr1-6", "Spr1_6":"Spr1-6", "Sq1":"Sq-1", "Sq_1":"Sq-1", "Sq8":"Sq-8", "Tamm2":"Tamm-2", "Tamm2B":"Tamm-2", "Tamm27":"Tamm-27", "Tamm_27":"Tamm-27", "Ts1":"Ts-1", "Ts1A":"Ts-1", "Ts5":"Ts-5", "Ts_5":"Ts-5", "Tsu1":"Tsu-1", "Tsu1B":"Tsu-1", "Ts1B":"Tsu-1", "Ull23":"Ull2-3", "Ull25":"Ull2-5", "Uod1":"Uod-1", "Uod_1":"Uod-1", "Uod7":"Uod-7", "Van0":"Van-0", "VanOA":"Van-0", "Var21":"Var2-1", "Var26":"Var2-6", "Var2_6":"Var2-6", "Wa1":"Wa-1","Wa_1":"Wa-1", "Wei0":"Wei-0", "Ws0":"Ws-0", "Ws2":"Ws-2", "Ws_2":"Ws-2", "Wt5":"Wt-5", "Yo0":"Yo-0", "Zdr1":"Zdr-1", "Zdr_1":"Zdr-1", "Zdr6":"Zdr-6"}


#accessionName2EcotypeId = {'NFA-8': '8346', 'Ei2': '8289', 'Pu27': '8362', 'Ms-0': '8340', 'Mt-0': '8341', 'Ull25': '8397', 'HR-10': '8308', 'Nfa10': '8345', 'Fei0': '8294', 'Ws0': '8405', 'Omo2_3': '8350', 'Est1B': '8291', 'Lp26': '8333', 'Var26': '8402', 'Zdr-1': '8409', 'Lp22': '8332', 'Ms_0': '8340', 'Tsu-1': '8394', 'Knox-18': '8318', 'HR5': '8309', 'Spr1-6': '8383', 'Knox-10': '8317', 'Spr1-2': '8382', 'Nok-3': '8347', 'Ra0': '8364', 'LP26': '8333', 'Bay0B': '8260', 'LP22': '8332', 'Mt0': '8341', 'Eden-2': '8287', 'Uod-1': '8398', 'LY1': '8279', 'FeiOB': '8294', 'LY3': '8324', 'Mrk0': '8339', 'C24B': '8273', 'C24A': '8273', 'Pro-0': '8360', 'Ra-0': '8364', 'Kin-0': '8316', 'Ts1A': '8392', 'Ts1B': '8394', 'Gu-0': '8301', 'LOV5A': '6046', 'Knox10': '8317', 'Mz-0': '8342', 'Ts5': '8393', 'Knox18': '8318', 'Sq_1': '8384', 'Ler1': '8324', 'Lz0': '8336', 'Ts1': '8392', 'Wa_1': '8403', 'Spr16': '8383', 'Fei-0': '8294', 'Nok_3': '8347', 'Bur-0': '8272', 'Spr1_6': '8383', 'Bur0': '8272', 'Gy0': '8302', 'Rrs10B': '8372', 'HR10': '8308', 'Tsu1': '8394', 'Ws_2': '8406', 'Ull2-3': '8396', 'Omo2-3': '8350', 'Omo2-1': '8349', 'Rmx-A180': '8371', 'Kondara': '8319', 'NFA10': '8345', 'Bil7': '8263', 'Tamm-2': '8390', 'Bil5': '8262', 'Kin0': '8316', 'CS22491': '8429', 'Nd1': '8344', 'Ga0': '8295', 'Cibc5': '8277', 'Ren11': '8368', 'Zdr-6': '8410', 'Bor4': '8268', 'Ren-11': '8368', 'Bor1': '5837', 'Var21': '8401', 'Est1': '8291', 'Lov_1': '6043', 'CIBC-17': '8276', 'Pu2-23': '8361', 'An-1': '8253', 'Zdr_1': '8409', 'Kas1': '8315', 'LL0': '8328', 'Ren-1': '8367', 'Ler-1': '8324', 'RRS10': '8372', 'Var2_6': '8402', 'Sq-8': '8385', 'Zdr6': '8410', 'EST1A': '8291', 'Rrs7B': '8373', 'Got-7': '8299', 'BroA': '8269', 'Uod_1': '8398', 'Tamm_27': '8391', 'Kz1': '8320', 'Wt5': '8407', 'Zdr1': '8409', 'Tamm-27': '8391', 'Kz9': '8322', 'RmxA02': '8370', 'Got7A': '8299', 'Ts_5': '8393', 'Eden2': '8287', 'Eden1': '6009', 'Pna10': '8358', 'Pu2-7': '8362', 'Tamm27': '8391', 'Mr-0': '8338', 'Ei-2': '8289', 'Mrk-0': '8339', 'Uod1': '8398', 'Pu2_7': '8362', 'CIBC17': '8276', 'Uod7': '8399', 'Lz-0': '8336', 'Gu0': '8301', 'Kz_9': '8322', 'Bil-7': '8263', 'Bil-5': '8262', 'Lp2-2': '8332', 'Fab-2': '8292', 'Fab-4': '8293', 'Se0': '8379', 'Sq-1': '8384', 'Oy-0': '8352', 'LL-0': '8328', 'Col0': '8279', 'Lov1': '6043', 'ShaA': '8248', 'Lov5': '6046', 'An1': '8253', 'Var2-1': '8401', 'Var2-6': '8402', 'CIBC5': '8277', 'Pna_17': '8359', 'Fab4': '8293', 'Got22': '8298', 'Wa-1': '8403', 'Fab2': '8292', 'Lov-5': '6046', 'Edi0': '8288', 'Lov-1': '6043', 'CS2249': '8429', 'RRS7': '8373', 'Ms0': '8340', 'Col-0': '8279', 'Gy-0': '8302', 'BUR0A': '8272', 'El2': '8289', 'Nok3': '8347', 'Wa1': '8403', 'Nd-1': '8344', 'Got-22': '8298', 'Br-0': '8269', 'HR-5': '8309', 'Cibc17': '8276', 'Ull2-5': '8397', 'Se0_2': '8379', 'Est-1': '8291', 'Tamm2B': '8390', 'Se-0': '8379', 'Pu223': '8361', 'Ct-1': '8280', 'Got7': '8299', 'Tamm2': '8390', 'Shahdara': '8248', 'Kas-2': '8315', 'Ct1': '8280', 'Cvi-0': '8281', 'Spr12': '8382', 'Omo21': '8349', 'Bor4B': '8268', 'Omo23': '8350', 'Mr0': '8338', 'Uod-7': '8399', 'Kno10': '8317', 'RmxA180': '8371', 'Bor-4': '8268', 'Bor-1': '5837', 'Kno18': '8318', 'Pna-10': '8358', 'Spr1_2': '8382', 'CviOA': '8281', 'Pna-17': '8359', 'C24': '8273', 'Ts-5': '8393', 'Nafa8B': '8346', 'Pna_10': '8358', 'Pna17': '8359', 'Bay-0': '8260', 'Br0': '8269', 'Ren1': '8367', 'NFA8': '8346', 'RRS-10': '8372', 'Eden-1': '6009', 'Yo-0': '8408', 'Ag0': '8251', 'Ren_11': '8368', 'Ts-1': '8392', 'Kz-1': '8320', 'RRS-7': '8373', 'Wt-5': '8407', 'Ws-2': '8406', 'Cvi0': '8281', 'Kz-9': '8322', 'Edi-0': '8288', 'Wei0': '8404', 'Wei-0': '8404', 'Oy0': '8352', 'Bay0': '8260', 'CIBC-5': '8277', 'Lov5B': '6046', 'Ws-0': '8405', 'Pro0': '8360', 'Van-0': '8400', 'Rmx-A02': '8370', 'Yo0': '8408', 'Van0': '8400', 'Ler1A': '8324', 'Sq8': '8385', 'Mz0': '8342', 'Lp2-6': '8333', 'Ull23': '8396', 'Ws2': '8406', 'Sq1': '8384', 'Ga-0': '8295', 'Tsu1B': '8394', 'BurOB': '8272', 'Kas2': '8315', 'NFA-10': '8345', 'Ag-0': '8251', 'VanOA': '8400', 'Sorbo': '8381'}


#def getEcotypeToAccessionDictionary(host="papaya.usc.edu", user=None, passwd=None, defaultValue=None):
#	class _ecotypeDict_(dict):
#		def __missing__(self, key):
#			return (defaultValue, defaultValue)
#
#	import MySQLdb
#	print "Connecting to db, host=" + host
#	if not user:
#		import sys
#		sys.stdout.write("Username: ")
#		user = sys.stdin.readline().rstrip()
#	if not passwd:
#		import getpass
#		passwd = getpass.getpass()
#	try:
#		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="at")
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#	#Retrieve the filenames
#	print "Fetching data"
#	numRows = int(cursor.execute("select distinct ei.tg_ecotypeid, e2a.accession_id, ei.nativename from stock.ecotypeid2tg_ecotypeid ei, at.complete_2010_strains_in_stock2tg_ecotypeid e2a where e2a.tg_ecotypeid=ei.tg_ecotypeid"))
#
#	ecotDict = _ecotypeDict_()
#	while(1):
#		row = cursor.fetchone()
#		if not row:
#			break;
#		ecotypeID = str(int(row[0]))
#		ecotDict[ecotypeID] = (str(int(row[1])), str(row[2]))
#	cursor.close ()
#	conn.close ()
#	return ecotDict
#
#def getEcotypeToNameDictionary(host="papaya.usc.edu", user=None, passwd=None, defaultValue=None):
#	class _ecotypeDict_(dict):
#		def __missing__(self, key):
#			return defaultValue
#
#	import MySQLdb
#	print "Connecting to db, host="+host
#	if not user:
#		import sys
#		sys.stdout.write("Username: ")
#		user = sys.stdin.readline().rstrip()
#	if not passwd:
#		import getpass
#		passwd = getpass.getpass()
#	try:
#		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#	#Retrieve the filenames
#	print "Fetching data"
#	numRows = int(cursor.execute("select distinct ei.id, ei.nativename from stock.ecotype ei order by ei.name"))
#	
#	ecotDict = _ecotypeDict_()
#	while(1):
#		row = cursor.fetchone()
#		if not row:
#			break;
#		ecotypeID = str(int(row[0]))
#		ecotDict[ecotypeID]=str(row[1])
#	cursor.close ()
#	conn.close ()
#	return ecotDict



def get250KDataFromDb(host="banyan.usc.edu", chromosomes=[1, 2, 3, 4, 5], db="stock_250k", withArrayIds=False,
			methodId=1, user=None, passwd=None, callProb=False, newBatch=False):
	"""
	Retrieve 2010 data from DB.  Returns a list of RawSnpsData objects. 
	
	withArrayInfo: if True, then a row of array used to genotype is added before the accession row.
	
	methodId: Specify the (base-calling/imputation) method.

	"""
	rt = time.time()
	decoder = RawDecoder()  #Other unused informative letters are ['R','Y','S','M','K','W']!!

	conn = dbutils.connect_to_default_lookup()
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"
	if newBatch:
		numRows = int(cursor.execute("select distinct ai.maternal_ecotype_id, ai.paternal_ecotype_id, ci.filename, ai.id from call_info ci, array_info ai where ci.array_id=ai.id and ci.method_id=" + str(methodId) + " and ai.maternal_ecotype_id is not null and ai.paternal_ecotype_id is not null and ai.id>492 order by ai.id"))  #An ugly hack..
	else:
		numRows = int(cursor.execute("select distinct ai.maternal_ecotype_id, ai.paternal_ecotype_id, ci.filename, ai.id from call_info ci, array_info ai where ci.array_id=ai.id and ci.method_id=" + str(methodId) + " and ai.maternal_ecotype_id is not null and ai.paternal_ecotype_id is not null order by ai.id"))
	print "Num rows:", numRows

	dict = {}
	accessions = []
	arrayIds = []
	dataFiles = []
	i = 0

	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		if int(row[0]) == int(row[1]):
			accession = str(int(row[0]))
		else:
			accession = str(int(row[0])) + "_" + str(int(row[1]))
		dict[accession] = i
		accessions.append(accession)
		dataFiles.append(row[2])
		arrayIds.append(str(int(row[3])))
		i = i + 1


	f = open(dataFiles[0], "r")
	lines = f.readlines()
	f.close()
	numSnps = len(lines)

	line = lines[1].split()

	linesList = []
	for i in range(0, len(accessions)): #Checking files
		f = open(dataFiles[i], "r")
		newLines = f.readlines()
		f.close()
		if len(newLines) != numSnps:
			raise Exception("Data files are not equal in length.")
		del newLines




	print "Reading and processing (cutting & merging) raw call files"
	import tempfile
	tempfile.tempdir = '/tmp'
	tmpFiles = []
	#Cutting files
	for i in range(0, len(dataFiles)): #Cutting the SNPs fields
		tmpFile = tempfile.mkstemp()
		os.close(tmpFile[0])
		st = "cut -f 2 " + str(dataFiles[i]) + " > " + tmpFile[1]
		#print st
		os.system(st)
		tmpFiles.append(tmpFile[1])


	def _mergeFiles_(fileList): #A recursive function that merges the files in the list
		l = len(fileList)
		if l > 2:
			tmpFile1 = _mergeFiles_(fileList[:l / 2])
			tmpFile2 = _mergeFiles_(fileList[l / 2:])
		elif l < 2:
			return fileList[0]
		elif l == 2:
			tmpFile1 = fileList[0]
			tmpFile2 = fileList[1]

		tmpFile = tempfile.mkstemp()
		os.close(tmpFile[0])
		st = "paste -d , " + str(tmpFile1) + " " + str(tmpFile2) + " > " + tmpFile[1] + ";\n"
		st += "rm " + str(tmpFile1) + " " + str(tmpFile2)
		#print st
		os.system(st)
		return tmpFile[1]


	print "Merging files"
	fileName = _mergeFiles_(tmpFiles)  #Merge all the fields
	print fileName


	chromosomes = []
	positionsList = []
	snpsList = []
	i = 1
	newChr = int(line[0].split("_")[0])

	if callProb:
		probTmpFiles = []
		#Cut the call probability fields
		print "Cutting call prob. files"
		for i in range(0, len(dataFiles)): #Cutting the SNPs fields
			tmpFile = tempfile.mkstemp()
			os.close(tmpFile[0])
			st = "cut -f 3 " + str(dataFiles[i]) + " > " + tmpFile[1]
			#print st
			os.system(st)
			probTmpFiles.append(tmpFile[1])

		print "Merging call prob. files"
		probFileName = _mergeFiles_(probTmpFiles)  #Merge all the fields
		print probFileName
		probFile = open(probFileName, "r")
		probFile.readline()
		callProbabilityList = []

		#Read file with call probabilities
		print "Starting to read file (snps and probabilities)"
		f = open(fileName, "r")
		f.readline()
		i = 1  #Skipping first line!
		while i < len(lines):
			chromosomes.append(int(newChr))
			oldChr = newChr
			positions = []
			snps = []
			callProbs = []
			while i < len(lines) and newChr == oldChr:
				line = lines[i].split()
				chrPos = line[0].split("_")
				oldChr = int(chrPos[0])
				positions.append(int(chrPos[1]))
				snp = []
				snpProbs = []
				newLine = f.readline().split(",")
				newProbLine = probFile.readline().split(",")
				for j in range(0, len(accessions)):
					nt = newLine[j].rstrip()
					prob = float(newProbLine[j].rstrip())
					snp.append(nt)
					snpProbs.append(prob)
				snps.append(snp)
				callProbs.append(snpProbs)
				i += 1
				if i < len(lines):
					line = lines[i].split()
					newChr = int(line[0].split("_")[0])
				else:
					break

			print "Fetched ", i, "SNPs of", len(lines)
			positionsList.append(positions)
			snpsList.append(snps)
			callProbabilityList.append(callProbs)
		probFile.close()
		os.system("rm " + probFileName)


	else: #Without call probabilities
		print "Starting to read file"
		f = open(fileName, "r")
		f.readline()
		while i < len(lines):
			chromosomes.append(int(newChr))
			oldChr = newChr
			positions = []
			snps = []
			while i < len(lines) and newChr == oldChr:
				line = lines[i].split()
				chrPos = line[0].split("_")
				oldChr = int(chrPos[0])
				positions.append(int(chrPos[1]))
				snp = []
				newLine = f.readline().split(",")
				for j in range(0, len(accessions)):
					nt = newLine[j].rstrip()
					snp.append(nt)
				snps.append(snp)
				i += 1
				if i < len(lines):
					line = lines[i].split()
					newChr = int(line[0].split("_")[0])
				else:
					break

			print "Fetched ", i, "SNPs of", len(lines)
			positionsList.append(positions)
			snpsList.append(snps)

	f.close()
	os.system("rm " + fileName)


	snpsds = []
	if callProb:
		for i in chromosomes:
			snpsds.append(RawSnpsData(snpsList[i - 1], positionsList[i - 1], accessions=accessions, arrayIds=arrayIds, callProbabilities=callProbabilityList[i - 1]))
			print "callProbabilityList length=", len(callProbabilityList[i - 1])
	else:
		for i in chromosomes:
			snpsds.append(RawSnpsData(snpsList[i - 1], positionsList[i - 1], accessions=accessions, arrayIds=arrayIds))
	dif = int(time.time() - rt)
	print "It took " + str(dif / 60) + " min. and " + str(dif % 60) + " sec. to fetch data."
	cursor.close ()
	conn.close ()
	return snpsds




def get149DataFromDb(host="papaya.usc.edu", chromosomes=[1, 2, 3, 4, 5], db="at", only96accessions=False, user=None, passwd=None):
	import MySQLdb
	"""
	Retrieve 149 data from DB.  Returns a list of RawSnpsData objects. 

	"""
	rt = time.time()
	decoder = RawDecoder()  #Other unused informative letters are ['R','Y','S','M','K','W']:

	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db=db)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	if only96accessions:
		numRows = int(cursor.execute("select distinct eva.ecotype_id, eva.accession_id, eva.nativename from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version=" + dataVersion + " and l.offset=0 and g.accession=eva.accession_id and eva.accession_id<97 and l.chromosome=1 order by eva.nativename"))
	else:
		numRows = int(cursor.execute("select distinct eva.ecotype_id, eva.accession_id, eva.nativename from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version=" + dataVersion + " and l.offset=0 and g.accession=eva.accession_id and l.chromosome=1 order by eva.nativename"))

	dict = {}
	accessions = []
	i = 0
	print "numRows", numRows
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		dict[int(row[1])] = i
		accessions.append(str(int(row[0])))
		i = i + 1

	print len(accessions), " accessions found."

	print "Fetching 2010 data (version " + dataVersion + "):"
	snpsds = []
	for chromosome in chromosomes:
		print "	Chromosome", chromosome
		if only96accessions:
			numRows = int(cursor.execute("select distinct l.position, g.accession, eva.nativename, al.base from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version>=" + dataVersion + " and l.offset=0 and g.accession=eva.accession_id and eva.accession_id<97 and l.chromosome=" + str(chromosome) + " order by l.position, eva.nativename"))
		else:
			numRows = int(cursor.execute("select distinct l.position, g.accession, eva.nativename, al.base from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version>=" + dataVersion + " and l.offset=0 and g.accession=eva.accession_id and l.chromosome=" + str(chromosome) + " order by l.position, eva.nativename"))
		print "	", numRows, "rows retrieved."
		positions = []
		snps = []
		if numRows > 0:
			row = cursor.fetchone()
			newPosition = int(row[0])
			while(1):
				if not row:
					break;
				print row
				positions.append(newPosition)
				oldPosition = newPosition
				snp = ['NA'] * len(accessions)  #Initialize to missing data.
				while(oldPosition == newPosition):
					snp[dict[int(row[1])]] = decoder[row[3]]
					row = cursor.fetchone()
					if not row:
						break;
					newPosition = int(row[0])
				snps.append(snp)

		snpsd = RawSnpsData(snps, positions, accessions=accessions)
		snpsds.append(snpsd)

	cursor.close ()
	conn.close ()
	dif = int(time.time() - rt)
	print "It took " + str(dif / 60) + " min. and " + str(dif % 60) + " sec. to fetch data."

	return snpsds



def get2010DataFromDb(host="papaya.usc.edu", chromosomes=[1, 2, 3, 4, 5], db="at", dataVersion="3", only96accessions=False, user=None, passwd=None):
	import MySQLdb, sys
	"""
	Retrieve 2010 data from DB.  Returns a list of RawSnpsData objects. 
	"""
	rt = time.time()
	decoder = RawDecoder({'-':'-'})  #Other unused informative letters are ['R','Y','S','M','K','W']:

	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db=db)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()

	locStr = " and l.offset=0 "
	if int(dataVersion) == 4:		   #If using version=4, then allow any locus.offset.
		locStr = ""

	#Get distinct accessions and their id.
	if only96accessions:
		print 'starting'
		numRows = int(cursor.execute("select distinct eva.ecotype_id, eva.accession_id, eva.nativename from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version=" + dataVersion + locStr + " and g.accession=eva.accession_id and eva.accession_id<97 order by eva.nativename"))
		print 'done'
	else:
		numRows = int(cursor.execute("select distinct eva.ecotype_id, eva.accession_id, eva.nativename from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version=" + dataVersion + locStr + " and g.accession=eva.accession_id order by eva.nativename"))

	dict = {}
	accessions = []
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		dict[int(row[1])] = i
		accessions.append(str(int(row[0])))
		i = i + 1

	print len(accessions), " accessions found."

	print "Fetching 2010 data (version " + dataVersion + "):"
	snpsds = []
	for chromosome in chromosomes:
		print "	Chromosome", chromosome
		if only96accessions:
			numRows = int(cursor.execute("select distinct l.position, g.accession, eva.nativename, al.base from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version>=" + dataVersion + locStr + " and g.accession=eva.accession_id and eva.accession_id<97 and l.chromosome=" + str(chromosome) + " order by l.position, eva.nativename"))
		else:
			numRows = int(cursor.execute("select distinct l.position, g.accession, eva.nativename, al.base from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version>=" + dataVersion + locStr + " and g.accession=eva.accession_id and l.chromosome=" + str(chromosome) + " order by l.position, eva.nativename"))
		print "	", numRows, "rows retrieved."
		positions = []
		snps = []
		if numRows > 0:
			row = cursor.fetchone()
			newPosition = int(row[0])
			while(1):
				if not row:
					break;
				positions.append(newPosition)
				oldPosition = newPosition
				snp = ['NA'] * len(accessions)  #Initialize to missing data.
				while(oldPosition == newPosition):
					snp[dict[int(row[1])]] = decoder[row[3]]
					row = cursor.fetchone()
					if not row:
						break;
					newPosition = int(row[0])
				snps.append(snp)

		snpsd = RawSnpsData(snps, positions, accessions=accessions)
		snpsd.alphabet.append('-')
		snpsds.append(snpsd)

	cursor.close ()
	conn.close ()
	dif = int(time.time() - rt)
	print "It took " + str(dif / 60) + " min. and " + str(dif % 60) + " sec. to fetch data."

	return snpsds



def get_2010_sequences_from_db(host="papaya.usc.edu", user="bvilhjal", passwd="*rri_bjarni@usc",
				chromosomes=[1, 2, 3, 4, 5], dataVersion="3"):
	"""
	Returns a list of AlignmentData objects.
	"""
	rt = time.time()
	decoder = RawDecoder({'-':'-'})  #Other unused informative letters are ['R','Y','S','M','K','W']:

	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()

	locStr = " and l.offset=0 "
	if int(dataVersion) == 4:		   #If using version=4, then allow any locus.offset.
		locStr = ""


	#Get distinct accessions and their id.

	numRows = int(cursor.execute("select distinct eva.accession_id, eva.nativename, eva.ecotype_id, e2tge.tg_ecotypeid from at.accession2tg_ecotypeid eva, stock.ecotypeid2tg_ecotypeid e2tge where eva.ecotype_id = e2tge.ecotypeid order by eva.accession_id"))
	aid_dict = {}
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		aid_dict[int(row[0])] = (row[1], row[2], row[3])

	numRows = int(cursor.execute("select distinct an.id, an.chromosome, an.start, an.end, an.version from at.alignment an where an.version>=" + dataVersion + locStr + " order by an.chromosome, an.start"))
	alignment_datas = []
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;

		a_id = int(row[0])
		a_chr = int(row[1])
		a_start = int(row[2])
		a_end = int(row[3])
		sequences = []
		positions = []
		offsets = []
		numRows = int(cursor.execute("select distinct l.id, l.position, l.offset, from at.locus l where l.alignment=" + str(a_id) + " order by l.position, l.offset"))




	print "Fetching 2010 data (version " + dataVersion + "):"
	snpsds = []
	for chromosome in chromosomes:
		print "	Chromosome", chromosome
		numRows = int(cursor.execute("select distinct l.position, g.accession, eva.nativename, al.base from at.genotype g, at.allele al, at.locus l, at.alignment an, at.accession2tg_ecotypeid eva where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version>=" + dataVersion + locStr + " and g.accession=eva.accession_id and l.chromosome=" + str(chromosome) + " order by l.position, eva.nativename"))
		print "	", numRows, "rows retrieved."
		positions = []
		snps = []
		if numRows > 0:
			row = cursor.fetchone()
			newPosition = int(row[0])
			while(1):
				if not row:
					break;
				positions.append(newPosition)
				oldPosition = newPosition
				snp = ['NA'] * len(accessions)  #Initialize to missing data.
				while(oldPosition == newPosition):
					snp[dict[int(row[1])]] = decoder[row[3]]
					row = cursor.fetchone()
					if not row:
						break;
					newPosition = int(row[0])
				snps.append(snp)

		snpsd = RawSnpsData(snps, positions, accessions=accessions)
		snpsd.alphabet.append('-')
		snpsds.append(snpsd)

	cursor.close ()
	conn.close ()
	dif = int(time.time() - rt)
	print "It took " + str(dif / 60) + " min. and " + str(dif % 60) + " sec. to fetch data."

	return snpsds




def getPerlgenDataFromDb(host="papaya.usc.edu", db="chip", chromosomes=[1, 2, 3, 4, 5], user=None, passwd=None):
	import MySQLdb
	"""
	Retrieve Perlgen data from DB.  Returns a list of RawSnpsData objects. 

	"""

	rt = time.time()
	perlgenTo2010 = {'bay-0':'Bay-0', 'bor-4':'Bor-4', 'br-0':'Br-0', 'bur-0':'Bur-0',
					 'c24':'C24', 'col-0':'Col-0', 'cvi-0':'Cvi-0', 'est-1':'Est-1',
					 'fei-0':'Fei-0', 'got-7':'Got-7', 'ler-1':'Ler-1', 'lov-5':'L\xf6v-5',
					 'nfa-8':'NFA-8', 'rrs-10':'RRS-10', 'rrs-7':'RRS-7', 'sha':'Sha',
					 'tamm-2':'Tamm-2', 'ts-1':'Ts-1', 'tsu-1':'Tsu-1', 'van-0':'Van-0'}

	perlgenAccessionToId = {}  #Get accession id's (or ecotype id?)  from DB
	e2n = getEcotypeToNameDictionary(host=host, user=user, passwd=passwd)
	n2e = {}
	#accessions = []
	for e in e2n:
		n2e[e2n[e]] = e
		#accessions.append(e2n[e])
	#print n2e
	#accessions.sort()
	#print accessions


	decoder = {'N':'NA', 'A':'A', 'C':'C', 'G':'G', 'T':'T'}
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db=db)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)

	cursor = conn.cursor()
	numRows = int(cursor.execute("select distinct d.ecotype from snp_combined_may_9_06_no_van d order by d.ecotype"))
	dict = {}
	accessions = []
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		dict[row[0]] = i
		accessions.append(n2e[perlgenTo2010[row[0]]])
		i = i + 1

	print "Fetching Perlgen data:"
	snpsds = []
	for chromosome in chromosomes:
		print "	Chromosome", chromosome
		numRows = int(cursor.execute("select distinct d.position, d.ecotype, d.mbml98 from snp_combined_may_9_06_no_van d where d.chromosome=" + str(chromosome) + " and d.mbml98 is not null order by d.position, d.ecotype;"))
		print "	", numRows, "rows retrieved."
		positions = []
		snps = []
		if numRows > 0:
			row = cursor.fetchone()
			newPosition = int(row[0])
			while(1):
				if not row:
					break;
				positions.append(newPosition)
				oldPosition = newPosition
				snp = ['NA'] * 20  #Initialize to missing data.
				while(oldPosition == newPosition):
					snp[dict[row[1]]] = decoder[row[2]]
					row = cursor.fetchone()
					if not row:
						break;
					newPosition = int(row[0])
				snps.append(snp)
		snpsd = RawSnpsData(snps, positions, accessions=accessions)
		snpsds.append(snpsd)

	cursor.close ()
	conn.close ()
	dif = int(time.time() - rt)
	print "It took " + str(dif / 60) + " min. and " + str(dif % 60) + " sec. to fetch data."
	return snpsds


def parseCSVDataAccessions(datafile, format=1, deliminator=",", missingVal='NA'):
	accessions = []
	f = open(datafile, 'r')
	line = f.readline()
 	withArrayIds = line == "Chromosome"
	if withArrayIds:
		line = line.split(deliminator)
		arrayIds = []
		for arrayId in line[2:]:
			arrayIds.append(arrayId.strip())
		line = f.readline()
	f.close()
	line = line.split(deliminator)
	for acc in line[2:]:
		accessions.append(acc.strip())
	print "Found", len(accessions), "arrays/strains."
	if withArrayIds:
		return accessions, arrayIds
	else:
		return accessions

def parseCSVData(datafile, format=1, deliminator=",", missingVal='NA', use_nt2number=0,
		returnChromosomes=False, useDecoder=True, filter=1, id=None, marker_type=None):
	"""
	05/12/08 yh. add argument use_nt2number, to turn nucleotide into numbers. default is not
	05/12/08 yh. add '-':'-' (deletion) and '|':'NA' (untouched) check nt2number in 'variation.common' to decoder
	05/12/2008 yh. more correct reporting of snp loading counts
	05/11/2008 yh. add chromosome. use RawSnpsData directly. ...
	06/25/2008 bjv. withArrayIds is now detected, i.e. the argument doesn't matter.	

	Parses raw CSV SNPs data files into a RawSnpsData.

	format=1: the function return a RawSnpsData object list
	format=0: the function return a SnpsData object list
	"""

	sys.stderr.write("Loading file: %s ... \n" % datafile)
	decoder = nt_decoder
	decoder[missingVal] = missing_val
	positions = [] #list[chr][position_index]
	genotypes = [] #list[chr][position_index][acces]
	accessions = []

	#Reading column data
	f = open(datafile, 'r')
	lines = f.readlines()
	f.close()

	chromosomes = []
	positionsList = []
	snpsList = []
	arrayIds = None
	line = lines[1].split(deliminator)
	i = 0
	if line[0] == "Chromosome":
		line = lines[i].split(deliminator)
		arrayIds = []
		for arrayId in line[2:]:
			arrayIds.append(arrayId.strip())
		i += 1
	line = lines[i].split(deliminator)
	for acc in line[2:]:
		accessions.append(acc.strip())
	print "Found", len(accessions), "arrays/strains."
	i += 1
	line = lines[i].split(deliminator)
	newChr = int(line[0])

	no_of_headers = i
	snpsd_ls = []
	while i < len(lines):
		chromosomes.append(int(newChr))
		oldChr = newChr
		rawSnpsData = RawSnpsData(accessions=accessions, arrayIds=arrayIds, id=id)	#05/11/2008 yh. use rawSnpsData
		rawSnpsData.snps = []
		rawSnpsData.positions = []
		rawSnpsData.chromosome = oldChr
		while i < len(lines) and newChr == oldChr:
			if filter < 1 and random.random() > filter:
	  			i += 1
				if i < len(lines):
					line = lines[i].split(deliminator)
					newChr = int(line[0])
			       	else:
					break
 				continue

			line = lines[i].split(deliminator)
			#print i,":",lines[i]
			oldChr = int(line[0])
			rawSnpsData.positions.append(int(line[1]))
			snp = []
			if useDecoder:
				for nt in line[2:]:
					snp.append(decoder[(nt.strip())])
			else:
				for nt in line[2:]:
					snp.append(nt.strip())
			rawSnpsData.snps.append(snp)
			i += 1
			if i < len(lines):
				line = lines[i].split(deliminator)
				newChr = int(line[0])
			else:
				break

		sys.stderr.write("Loaded %s of %s SNPs.\n" % (i - no_of_headers, len(lines) - no_of_headers))
		#Adding marker
		if marker_type:
			marker_types = [marker_type] * len(rawSnpsData.snps)
			rawSnpsData.marker_types = marker_types
		snpsd_ls.append(rawSnpsData)
		del rawSnpsData
	if format == 0:
		print "Converting raw SNPs data to binary SNPs."
		for i in range(0, len(chromosomes)):
			snpsd_ls[i] = snpsd_ls[i].getSnpsData()
	sys.stderr.write("\n")
	if returnChromosomes:
		return (snpsd_ls, chromosomes)
	else:
		return snpsd_ls






def parseCSVDataWithCallProb(datafile, callProbFile, format=1, deliminator=",", missingVal='NA', withArrayIds=False):
	"""
	Parses raw CSV SNPs data files into a RawSnpsData.

	format=1: the function return a RawSnpsData object list
	format=0: the function return a SnpsData object list
	"""
	sys.stderr.write("Loading file: %s ...\n" % datafile)
	decoder = {missingVal:'NA', 'A':'A', 'C':'C', 'G':'G', 'T':'T'}

	positions = [] #list[chr][position_index]
	genotypes = [] #list[chr][position_index][acces]
	accessions = []

	#Reading column data
	f = open(datafile, 'r')
	lines = f.readlines()
	f.close()
	numLines = len(lines)
	f = open(datafile, 'r')
	probFile = open(callProbFile, 'r')

	chromosomes = []
	positionsList = []
	snpsList = []
	accessions = []
	callProbList = []
	arrayIds = None
	i = 0
	if withArrayIds:
		line = (f.readline()).split(deliminator)
		probFile.readline()
		arrayIds = []
		for arrayId in line[2:]:
			arrayIds.append(arrayId.strip())
		i += 1
	line = (f.readline()).split(deliminator)
	probFile.readline()
	for acc in line[2:]:
		accessions.append(acc.strip())
	i += 1

	newChr = lines[3].split(deliminator)[0]
	while i < numLines:
		chromosomes.append(int(newChr))
		oldChr = newChr
		positions = []
		snps = []
		probsList = []
		while i < numLines and newChr == oldChr:
			line = (f.readline()).split(deliminator)
			probLine = (probFile.readline()).split(deliminator)
			oldChr = int(line[0])
			positions.append(int(line[1]))
			snp = []
			probs = []
			for j in range(2, len(line)):
				nt = line[j]
				snp.append(decoder[nt.strip()])
				probs.append(float(probLine[j]))
			snps.append(snp)
			probsList.append(probs)
			i += 1
			if i < numLines:
				line = lines[i].split(deliminator)
				newChr = int(line[0])
			else:
				break
		sys.stderr.write("Loaded %s of %s SNPs.\n" % (i, len(lines)))
		positionsList.append(positions)
		snpsList.append(snps)
		callProbList.append(probsList)
	f.close()
	probFile.close()

	snpsds = []
	for i in range(0, len(chromosomes)):
		snpsds.append(RawSnpsData(snpsList[i], positionsList[i], accessions=accessions, arrayIds=arrayIds, callProbabilities=callProbList[i]))
	if format == 0:
		for i in range(0, len(chromosomes)):
			snpsds[i] = snpsds[i].getSnpsData()
	sys.stderr.write("\n")
	return(snpsds)



def parse2010Data(datafile=None):
	"""
	Returns 2010 Data.
	
	Loads it from a file.
	"""

	#Imputed data
	datadir = homedir + "Projects/data/2010-Zhao_et_al/"
	if not datafile:
		datafile = datadir + "SNPs.csv"

	#Initialization
	accessions = []
	positions = [] #list[chr][position_index]
	genotypes = [] #list[chr][position_index][acces]
	for i in range(0, 5):
		positions.append([])


	chromosomes = []
	f = open(datafile, 'r')
	lines = f.readlines()
	f.close()
	for chr in (lines[0].strip()).split(",")[1:]:
		chromosomes.append(int(chr) - 1)

	line = (lines[1].strip()).split(",")[1:]
	for i in range(0, len(chromosomes)):
		positions[chromosomes[i]].append(int(line[i]))

	for j in range(0, 5):
		l = []
		for k in range(0, len(positions[j])):
			l.append([])
		genotypes.append(l)


	accessions = []
	for i in range(2, len(lines)):
		line = (lines[i].strip()).split(",")
		accessions.append(line[0])
		line = line[1:]
		for j in range(0, 5):
			for k in range(0, len(positions[j])):
				genotypes[j][k].append(line[i])

	print accessions

	import random
	#Converting genotype and filtering.
	countAll = 0
	countGood = 0
	countStupid = 0
	decoder = {'.':'NA'}
	newgenotypes = [[], [], [], [], []]
	newpositions = [[], [], [], [], []]
	for i in range(0, len(genotypes)):
		for j in range(0, len(positions[i])):
			countAll = countAll + 1
			k = 0
			ntl = [] #list of observed nucleotides.
			for nt in ['0', '1', '2', '3']:
				if nt in genotypes[i][j]:
					decoder[nt] = k
					ntl.append(nt)
					k = k + 1
			if k == 2:
				countGood = countGood + 1
				l = []
				for nt in genotypes[i][j]:
					l.append(decoder[nt])
				newgenotypes[i].append(l)
				newpositions[i].append(positions[i][j])
			else:
				if k == 1:
					countStupid = countStupid + 1
	print countAll, " SNPs in all"
	print countGood, " SNPs used"
	print countStupid, " Stupid SNPs thrown away"

	for i in range(0, 5):
		print newpositions[i][-1]
	del positions
	del genotypes
	positions = newpositions
	genotypes = newgenotypes


	chromasomes = []
	for i in range(0, 5):
		chromasomes.append(SnpsData(genotypes[i], positions[i], accessions=accessions))

	#print positions[4]
	return(chromasomes)


#def parse250DataRaw(imputed = True):
#	"""
#	WARNING: OUTDATED
#	Returns 250K Data in a list of RawSnpsData objects (one for each chromosome).
#
#	Set imputed to False, if un-imputed data is preferred.
#
#	"""
#
#	if imputed:
#		datadir = homedir+"Projects/data/NPUTE_results/"
#		datafile = datadir+"250K_snp2acc_NPUTE.csv"
#		datafilecol = datadir+"250K_snp2acc_NPUTE.colnames.csv"
#		datafilerow = datadir+"250K_snp2acc_NPUTE.rownames.csv"
#	else:
#		datadir = homedir+"Projects/data/020808dump/"
#		datafile = datadir+"250K_snp2acc_GOOD.csv"
#		datafilecol = datadir+"250K_snp2acc_GOOD.colnames.csv"
#		datafilerow = datadir+"250K_snp2acc_GOOD.rownames.csv"
#		
#	positions = [] #list[chr][position_index]
#	genotypes = [] #list[chr][position_index][acces]
#	accessions = []
#	
#	phenotype = None
#	filteredPhenotypes = [] #list[acces]
#	individuals = [] #list[acces][chr][position_index]  Filtered for specific phenotypes.
#	
#	#Reading column data
#	f = open(datafilecol, 'r')
#	lines = f.readlines()
#	f.close()
#	for line in lines:
#		if imputed:
#			acc = line.strip()
#		else:
#			acc = (line.strip().split("."))[1]			
#		accessions.append(accessionName2Id[accessions250To2010[acc]])
#
#	#Reading row data
#	f = open(datafilerow, 'r')
#	lines = f.readlines()
#	f.close()
#	positions = [[],[],[],[],[]] #1 list per chromasome
#	for line in lines:
#		if imputed:
#			l = line.strip().split(".")
#			chr = int(l[0])
#			pos = int(l[1])
#		else:
#			pos = int(line.strip())
#			chr = pos/100000000
#			pos = pos%100000000
#		positions[chr-1].append(pos)
#
#	#Reading genotype data
#	rawgenotypes = [[],[],[],[],[]] #1 list per chromasome
#	f = open(datafile, 'r')
#	for i in range(0,5):
#		for p in positions[i]:
#			line = f.readline()
#			l = (line.strip()).split(",")
#			if l == []:
#				raise Exception("Data problem") 
#			rawgenotypes[i].append(l)
#	f.close()
#	print "raw genotype read"
#
#	import random
#	#Converting genotype and filtering.
#	decoder = RawDecoder()
#	genotypes = [[],[],[],[],[]]
#	newpositions = [[],[],[],[],[]]
#	for i in range(0,len(rawgenotypes)):
#		for j in range(0,len(positions[i])):
#			l = []
#			for nt in rawgenotypes[i][j]:
#				l.append(decoder[nt])
#			genotypes[i].append(l)
#			newpositions[i].append(positions[i][j])
#		print "Chromosome",i+1,":",len(newpositions[i]),"SNPs."
#
#	"""
#	for i in range(0,5):
#		print newpositions[i][-1]
#	"""
#	del positions
#	del rawgenotypes
#	positions = newpositions
#	
#	chromasomes = []
#	for i in range(0,5):
#		chromasomes.append(RawSnpsData(genotypes[i],positions[i],accessions=accessions))
#			
#	return(chromasomes)



#def parse250KDataFiles(imputed = True):
#	"""
#	WARNING: OUTDATED
#	Returns 250K Data as a list of SnpsData objects (not RawSnpsData).
#	
#	Set imputed to False, if un-imputed data is preferred.
#
#	"""
#	
#	if imputed:
#		datadir = homedir+"Projects/data/NPUTE_results/"
#		datafile = datadir+"250K_snp2acc_NPUTE.csv"
#		datafilecol = datadir+"250K_snp2acc_NPUTE.colnames.csv"
#		datafilerow = datadir+"250K_snp2acc_NPUTE.rownames.csv"
#	else:
#		datadir = homedir+"Projects/data/020808dump/"
#		datafile = datadir+"250K_snp2acc_GOOD.csv"
#		datafilecol = datadir+"250K_snp2acc_GOOD.colnames.csv"
#		datafilerow = datadir+"250K_snp2acc_GOOD.rownames.csv"
#		
#	positions = [] #list[chr][position_index]
#	genotypes = [] #list[chr][position_index][acces]
#	accessions = []
#	
#	phenotype = None
#	filteredPhenotypes = [] #list[acces]
#	individuals = [] #list[acces][chr][position_index]  Filtered for specific phenotypes.
#	
#	#Reading column data
#	f = open(datafilecol, 'r')
#	lines = f.readlines()
#	f.close()
#	for line in lines:
#		if imputed:
#			acc = line.rstrip()
#		else:
#			acc = (line.rstrip().split("."))[1]			
#		accessions.append(acc)
#
#	#Reading row data
#	f = open(datafilerow, 'r')
#	lines = f.readlines()
#	f.close()
#	positions = [[],[],[],[],[]] #1 list per chromasome
#	for line in lines:
#		if imputed:
#			l = line.rstrip().split(".")
#			chr = int(l[0])
#			pos = int(l[1])
#		else:
#			pos = int(line.rstrip())
#			chr = pos/100000000
#			pos = pos%100000000
#		positions[chr-1].append(pos)
#
#	#Reading genotype data
#	rawgenotypes = [[],[],[],[],[]] #1 list per chromasome
#	f = open(datafile, 'r')
#	for i in range(0,5):
#		for p in positions[i]:
#			line = f.readline()
#			l = (line.rstrip()).split(",")
#			if l == []:
#				raise Exception("Data problem") 
#			rawgenotypes[i].append(l)
#	f.close()
#	print "raw genotype read"
#
#	import random
#	#Converting genotype and filtering.
#	countAll = 0
#	countGood = 0 
#	countStupid = 0
#	decoder = {'NA':'NA'}
#	genotypes = [[],[],[],[],[]]
#	newpositions = [[],[],[],[],[]]
#	for i in range(0,len(rawgenotypes)):
#		for j in range(0,len(positions[i])):
#			countAll = countAll+1
#			k = 0
#			ntl = [] #list of observed nucleotides.
#			for nt in ['A','C','G','T']:
#				if nt in rawgenotypes[i][j]:
#					decoder[nt]=k
#					ntl.append(nt)
#					k = k+1
#			if k==2:
#				countGood = countGood + 1
#				l = []
#				for nt in rawgenotypes[i][j]:
#					l.append(decoder[nt])
#				genotypes[i].append(l)
#				newpositions[i].append(positions[i][j])
#			else:
#				if k==1:
#					countStupid = countStupid+1
#	print countAll," SNPs in all"
#	print countGood," SNPs used"
#	print countStupid," Stupid SNPs thrown away"
#
#	for i in range(0,5):
#		print newpositions[i][-1]
#	del positions
#	del rawgenotypes
#	positions = newpositions
#	
#	chromasomes = []
#	for i in range(0,5):
#		chromasomes.append(SnpsData(genotypes[i],positions[i],accessions=accessions))
#			
#	return(chromasomes)


#def parseMSFile(filename):
#	"""
#	Parses a Hudson's ms file.
#
#	Returns a list of snpsdata.
#	"""
#	import random
#	f = open(filename, 'r')
#	lines = f.readlines()
#	f.close()
#	i = 0
#	data = []
#	while i	 < len(lines):
#		line = lines[i]
#		if line.startswith("//"):
#			num = 0
#			positions = []
#			snps = []
#			i = i + 1
#			while i < len(lines) and not lines[i].startswith("//"):
#				line = lines[i]
#				if line.startswith("segsites:"):
#					num = int(line[9:].rstrip())
#				if line.startswith("positions:"):
#					l1 = line[10:].rstrip().split()
#					l2 = [0.0] * len(l1)
#					for j in range(0, len(l1)):
#						l2[j] = float(l1[j])
#						snps.append([])		#Initializing the snps.
#					positions = l2
#				if line[0].isdigit():
#					line = line.rstrip()
#					snps[0].append(int(line[0]))
#					for j in range(1, len(positions)):
#						snps[j].append(int(line[j]))
#
#				i = i + 1
#			newSnps = []
#			newPositions = []
#			if len(positions) > 0:
#				newSnps.append(snps[0])
#				newPositions.append(positions[0])
#			k = 0
#			for j in range(1, len(positions)):
#				newSnps.append(snps[j])
#				newPositions.append(positions[j])
#
#			newSnpsd = SnpsData(newSnps, newPositions)
#			#print newSnpsd.snps
#			data.append(newSnpsd)
#		else:
#			i = i + 1
#	#print len(data[0].snps)
#	return data
#
#
#
#def parseMSFileBasescale(filename, baseScale=100000):
#	"""
#	Parses a Hudson's ms file.
#
#	Returns a list of snpsdata.
#	"""
#	import random
#	f = open(filename, 'r')
#	lines = f.readlines()
#	f.close()
#	i = 0
#	data = []
#	invBaseScale = 1 / float(baseScale)
#	while i	 < len(lines):
#		line = lines[i]
#		if line.startswith("//"):
#			num = 0
#			positions = []
#			snps = []
#			i = i + 1
#			while i < len(lines) and not lines[i].startswith("//"):
#				line = lines[i]
#				if line.startswith("segsites:"):
#					num = int(line[9:].rstrip())
#				if line.startswith("positions:"):
#					l1 = line[10:].rstrip().split()
#					l2 = [0.0] * len(l1)
#					for j in range(0, len(l1)):
#						l2[j] = float(l1[j])
#						snps.append([])		#Initializing the snps.
#					positions = l2
#				if line[0].isdigit():
#					line = line.rstrip()
#					snps[0].append(int(line[0]))
#					for j in range(1, len(positions)):
#						snps[j].append(int(line[j]))
#
#				i = i + 1
#			newSnps = []
#			newPositions = []
#			if len(positions) > 0:
#				newSnps.append(snps[0])
#				newPositions.append(positions[0])
#			k = 0
#			for j in range(1, len(positions)):
#				if (positions[j] - positions[j - 1]) >= invBaseScale:
#					newSnps.append(snps[j])
#					newPositions.append(positions[j])
#				else:
#					if baseScale > 10000 and k <= j:  #The precision of ms is 1/10000 (sometimes we have many more snps)
#						k = j + 1
#						while(k < len(positions) and (positions[k] - positions[j - 1]) < invBaseScale):
#							k = k + 1
#						last = 0.0
#						for h in range(0, k - j): # Order statistic 
#							last = random.betavariate(1, k - j - h) * (1 - last) + last
#							positions[j + h] = positions[j + h] + last / 10000
#						if (positions[j] - positions[j - 1]) >= invBaseScale:
#							newSnps.append(snps[j])
#							newPositions.append(positions[j])
#
#			newSnpsd = SnpsData(newSnps, newPositions, baseScale)
#			#print newSnpsd.snps
#			data.append(newSnpsd)
#		else:
#			i = i + 1
#	#print len(data[0].snps)
#	return data
#
#
#def parseMSFileFilter(filename, baseScale=1000000, fixPos=True, filterProb=1.0):
#	f = open(filename, 'r')
#	lines = 	f.readlines()
#	f.close()
#	i = 0
#	data = []
#	invBaseScale = 1 / float(baseScale)
#	while i	 < len(lines):
#		line = lines[i]
#		if line.startswith("//"):
#			num = 0
#			positions = []
#			snps = []
#			i = i + 1
#			while i	 < len(lines) and not lines[i].startswith("//"):
#				line = lines[i]
#				if line.startswith("segsites:"):
#					num = int(line[9:].rstrip())
#				if line.startswith("positions:"):
#					l1 = line[10:].rstrip().split()
#					l2 = [0.0] * len(l1)
#					for j in range(0, len(l1)):
#						l2[j] = float(l1[j])
#						snps.append([])		#Initializing the snps.
#					positions = l2
#				if line[0].isdigit():
#					line = line.rstrip()
#					snps[0].append(int(line[0]))
#					for j in range(1, len(positions)):
#						snps[j].append(int(line[j]))
#
#				i = i + 1
#			newSnps = []
#			newPositions = []
#			l = 0
#			pos = -1
#			for j in range(0, len(positions)):
#				if random.random() <= filterProb:
#					npos = positions[j]
#					if (npos - pos) >= invBaseScale:
#						newSnps.append(snps[j])
#						newPositions.append(npos)
#					else:
#						if fixPos and baseScale > 10000 and l <= j :  #The precision of ms
#							l = j + 1
#							while(l < len(positions) and (positions[l] - pos) < invBaseScale):
#								l = l + 1
#							last = 0.0
#							for h in range(0, l - j): # Order statistic 
#								last = random.betavariate(1, l - j - h) * (1 - last) + last
#								positions[j + h] = positions[j + h] + last / 10000
#
#							if (npos - pos) >= invBaseScale:
#								newSnps.append(snps[j])
#								newPositions.append(npos)
#					pos = positions[j]
#			data.append(SnpsData(newSnps, newPositions, baseScale))
#		else:
#			i = i + 1
#	return data
#
#
#
#def parseMSData(filename, baseScale=1000000, sampleNum=None, fixPos=True):
#	""" 
#	Parses randomly chosen ms outputs from previously calculated files.  Depends on the datafiles!
#	"""
#	f = open(filename, 'r')
#	lines = 	f.readlines()
#	totSampleNum = int(lines[0].split()[2])
#	f.close()
#	i = 0
#	klist = []
#	if sampleNum:
#		k = 0
#		while k < sampleNum:
#			val = int(random.random()*totSampleNum + 1)
#			if not val in klist:
#				klist.append(val)
#				k = k + 1
#	else:
#		klist = range(0, totSampleNum)
#	k = -1
#	data = []
#	invBaseScale = 1 / float(baseScale)
#	while i < len(lines):
#		line = lines[i]
#		if line.startswith("//"):
#			k = k + 1
#			i = i + 1
#			if k in klist:
#				num = 0
#				positions = []
#				snps = []
#				while i	 < len(lines) and not lines[i].startswith("//"):
#					line = lines[i]
#					if line.startswith("segsites:"):
#						num = int(line[9:].rstrip())
#					if line.startswith("positions:"):
#						l1 = line[10:].rstrip().split()
#						l2 = [0.0] * len(l1)
#						for j in range(0, len(l1)):
#							l2[j] = float(l1[j])
#							snps.append([])		#Initializing the snps.
#						positions = l2
#					if line[0].isdigit():
#						line = line.rstrip()
#						snps[0].append(int(line[0]))
#						for j in xrange(1, len(positions)):
#							snps[j].append(int(line[j]))
#
#					i = i + 1
#				newSnps = []
#				newPositions = []
#				if len(positions) > 0:
#					newSnps.append(snps[0])
#					newPositions.append(positions[0])
#					pos = positions[0]
#				l = 0
#				for j in range(1, len(positions)):
#					npos = positions[j]
#					if (npos - pos) >= invBaseScale:
#						newSnps.append(snps[j])
#						newPositions.append(npos)
#					else:
#						if fixPos and baseScale > 10000 and l <= j :  #The precision of ms
#							l = j + 1
#							while(l < len(positions) and (positions[l] - pos) < invBaseScale):
#								l = l + 1
#							last = 0.0
#							for h in range(0, l - j): # Order statistic 
#								last = random.betavariate(1, l - j - h) * (1 - last) + last
#								positions[j + h] = positions[j + h] + last / 10000
#
#							if (npos - pos) >= invBaseScale:
#								newSnps.append(snps[j])
#								newPositions.append(npos)
#					pos = positions[j]
#				data.append(SnpsData(newSnps, newPositions, baseScale))
#		else:
#			i = i + 1
#	return data
#
#
#def parseMSDataFilter(filename, baseScale=1000000, sampleNum=None, fixPos=True, filterProb=1.0):
#	""" 
#	Parses randomly chosen ms outputs from previously calculated files.  Depends on the datafiles!
#	"""
#	#print filename
#	f = open(filename, 'r')
#	lines = 	f.readlines()
#	totSampleNum = int(lines[0].split()[2])
#	f.close()
#	i = 0
#	klist = []
#	if sampleNum:
#		k = 0
#		while k < sampleNum:
#			val = int(random.random()*totSampleNum + 1)
#			if not val in klist:
#				klist.append(val)
#				k = k + 1
#	else:
#		klist = range(0, totSampleNum)
#	k = -1
#	data = []
#	invBaseScale = 1 / float(baseScale)
#	while i	 < len(lines):
#		line = lines[i]
#		if line.startswith("//"):  #New data
#			k = k + 1
#			i = i + 1
#			if k in klist:
#				num = 0
#				positions = []
#				snps = []
#				while i	 < len(lines) and not lines[i].startswith("//"):
#					line = lines[i]
#					if line.startswith("segsites:"):
#						num = int(line[9:].rstrip())
#					if line.startswith("positions:"):
#						l1 = line[10:].rstrip().split()
#						l2 = [0.0] * len(l1)
#						for j in range(0, len(l1)):
#							l2[j] = float(l1[j])
#							snps.append([])		#Initializing the snps.
#						positions = l2
#					if line[0].isdigit():
#						if len(positions) == len(line.rstrip()):  #Hack added to fix ms generated errors in files!
#							line = line.rstrip()
#								#snps[0].append(int(line[0]))
#							for j in xrange(0, len(line)):
#								snps[j].append(int(line[j]))
#						else:
#							print "line nr.", i
#							print "positions:", positions
#							print "line:", line.rstrip()
#							print "old line:", line
#							raise Exception()
#					i = i + 1
#				newSnps = []
#				newPositions = []
#				l = 0
#				pos = -1
#				debug1 = positions[:]
#				debug2 = snps[:]
#				r = 0
#				for j in range(0, len(positions)):
#					if random.random() <= filterProb:
#						r = r + 1
#						npos = positions[j]
#						if (npos - pos) >= invBaseScale:
#							newSnps.append(snps[j])
#							newPositions.append(npos)
#						else:
#							if fixPos and baseScale > 10000 and l <= j :  #The precision of ms
#								l = j + 1
#								while(l < len(positions) and (positions[l] - pos) < invBaseScale):
#									l = l + 1
#								last = 0.0
#								for h in range(0, l - j): # Order statistic 
#									last = random.betavariate(1, l - j - h) * (1 - last) + last
#									positions[j + h] = positions[j + h] + last / 10000
#
#								if (npos - pos) >= invBaseScale:
#									newSnps.append(snps[j])
#									newPositions.append(npos)
#						pos = positions[j]
#				if r == 0 :
#					print "No luck"
#					print len(positions)
#				if len(newSnps) == 0 or len(newSnps[0]) == 0:
#					print "data is empty!!!"
#					print "It succeded ", r, "number of times."
#					print debug1
#					print newPositions
#					print newSnps
#					print debug2
#				snpsd = SnpsData(newSnps, newPositions, baseScale)
#				if len(snpsd.snps) == 0 or len(snpsd.snps[0]) == 0:
#					print snpsd.positions
#					print snpsd.snps
#				data.append(snpsd)
#		else:
#			i = i + 1
#	if data[0] == []:
#		print
#	return data


def parse_numerical_snp_data(data_file, delimiter=",", missing_val='NA', filter=1, filter_accessions=None,
			use_pickle=True, dtype='int8', data_format='binary'):
	"""
	A sped-up version, to load a int (e.g. binary) file.
	
	If pickle is used then more memory is required, but it speeds up.
	"""
	import numpy as np
	import random
	import time
	pickle_file = data_file + '.pickled'
	if use_pickle:
		if os.path.isfile(pickle_file):
			t_ = time.time()
			f = open(pickle_file, 'rb')
			sd = cPickle.load(f)
			f.close()
			if filter_accessions and set(filter_accessions) != set(sd.accessions):
				print 'Filtering accessions.'
				sd.filter_accessions(filter_accessions)
			if filter < 1:
				sd.sample_snps(filter)
			print 'Loading genotype data took %.2f s...' % (time.time() - t_)
			return sd
		else:
			filter_accessions_ = filter_accessions
			filter_accessions = None
	sys.stderr.write("Loading binary SNPs data file: %s \n" % data_file)
	accessions = []
	snpsd_ls = []
	chromosomes = []

	#Reading accession data
	f = open(data_file, 'rb')
	line = f.next().split(delimiter)
	for acc in line[2:]:
		accessions.append(acc.strip())
	print "Found", len(accessions), "arrays/strains."

	if filter_accessions:
		indices_to_include = []
		new_accessions = []
		for i, acc in enumerate(accessions):
			if acc in filter_accessions:
				indices_to_include.append(i)
				new_accessions.append(acc)
		accessions = new_accessions
		print "Loading only", len(accessions), "arrays/strains."

	num_accessions = len(accessions)

	positions = []
	snps = []
	line = f.next().split(delimiter)
	old_chromosome = int(line[0])
	positions.append(int(line[1]))
	snps.append(np.array(line[2:], dtype=dtype))
	#snps.append(map(int,line[2:]))
	i = 0

	if filter < 1.0:
		if filter_accessions:
			for i, line in enumerate(f):
				if random.random() > filter: continue
				line_list = line.split(delimiter)
				chromosome = int(line_list[0])
				if chromosome != old_chromosome: #Then save snps_data
					snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
					snpsd_ls.append(snpsd)
					chromosomes.append(old_chromosome)
					positions = []
					snps = []
					old_chromosome = chromosome
					sys.stderr.write("Loaded %s SNPs.\n" % i)
				positions.append(int(line_list[1]))
				l_list = line_list[2:]
				#snp = [int(l_list[j]) for j in indices_to_include]
				snp = np.array([l_list[j] for j in indices_to_include], dtype=dtype)
				snps.append(snp)
		else:
			for i, line in enumerate(f):
				if random.random() > filter: continue
				line_list = line.split(delimiter)
				chromosome = int(line_list[0])
				if chromosome != old_chromosome: #Then save snps_data
					snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
					snpsd_ls.append(snpsd)
					chromosomes.append(old_chromosome)
					positions = []
					snps = []
					old_chromosome = chromosome
					sys.stderr.write("Loaded %s SNPs.\n" % i)
				positions.append(int(line_list[1]))
				#snps.append(map(int,line_list[2:]))
				snps.append(np.array(line_list[2:], dtype=dtype))
	else:
		if filter_accessions:
			for i, line in enumerate(f):
				line_list = line.split(delimiter)
				chromosome = int(line_list[0])
				if chromosome != old_chromosome: #Then save snps_data
					snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
					snpsd_ls.append(snpsd)
					chromosomes.append(old_chromosome)
					positions = []
					snps = []
					old_chromosome = chromosome
					sys.stderr.write("Loaded %s SNPs.\n" % i)
				positions.append(int(line_list[1]))
				l_list = line_list[2:]
				#snp = [int(l_list[j]) for j in indices_to_include]
				snp = np.array([l_list[j] for j in indices_to_include], dtype=dtype)
				snps.append(snp)
		else:
			for i, line in enumerate(f):
				line_list = line.split(delimiter)
				chromosome = int(line_list[0])
				if chromosome != old_chromosome: #Then save snps_data
					snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
					snpsd_ls.append(snpsd)
					chromosomes.append(old_chromosome)
					positions = []
					snps = []
					old_chromosome = chromosome
					sys.stderr.write("Loaded %s SNPs.\n" % i)
				positions.append(int(line_list[1]))
				#snps.append(map(int,line_list[2:]))
				snps.append(np.array(line_list[2:], dtype=dtype))



	f.close()
	snpsd = SNPsData(snps, positions, accessions=accessions, chromosome=old_chromosome)
	snpsd_ls.append(snpsd)
	chromosomes.append(old_chromosome)
	sys.stderr.write("Loaded %s SNPs.\n" % i)
	sd = SNPsDataSet(snpsd_ls, chromosomes, data_format=data_format)
	if use_pickle:
		print 'Saving a pickled version of genotypes.'
		f = open(pickle_file, 'wb')
		cPickle.dump(sd, f, protocol=2)
		f.close()
	if filter_accessions_:
		sd.filter_accessions(filter_accessions_)
	return sd



def parse_snp_data(data_file, delimiter=",", missingVal='NA', format='nucleotides', filter=1,
		id=None, useDecoder=True, look_for_binary=True, filter_accessions=None,
		use_pickle=True):
	"""
	Load snps data..
	"""
	data_format_dict = {'binary':0, 'nucleotides':1, 'int':2, 'float':3}
	if format == 'binary' and look_for_binary:
		print 'Looking for binary SNPs'
		sd_binary_file = data_file + '.binary'
		if os.path.isfile(sd_binary_file):
			sd = parse_numerical_snp_data(sd_binary_file, delimiter=delimiter, missing_val=missingVal,
						filter=filter, filter_accessions=filter_accessions,
						use_pickle=use_pickle, dtype='int8', data_format=format)
		else: #Try nucleotide format
			sd = parse_snp_data(data_file , format='nucleotides', delimiter=delimiter,
					      missingVal=missingVal, filter=filter, look_for_binary=False,
					      filter_accessions=filter_accessions)
			sd.convert_2_binary()
			print 'Save a binary snps data file:', sd_binary_file
			sd.writeToFile(sd_binary_file, binary_format=True)
			if use_pickle:
				pickle_file = sd_binary_file + '.pickled'
				print 'Saving a pickled version of genotypes.'
				f = open(pickle_file, 'wb')
				cPickle.dump(sd, f, protocol=2)
				f.close()
	elif format == 'int':
		sd = parse_numerical_snp_data(data_file, delimiter=delimiter, missing_val=missingVal,
					filter=filter, filter_accessions=filter_accessions,
					use_pickle=use_pickle, dtype='int8', data_format=format)
	elif format == 'float':
		sd = parse_numerical_snp_data(data_file, delimiter=delimiter, missing_val=missingVal,
					filter=filter, filter_accessions=filter_accessions,
					use_pickle=use_pickle, dtype='float32', data_format=format)
	elif format == 'nucleotides':
		print 'Looking for nucleotide SNPs'
		(snpsds, chromosomes) = parseCSVData(data_file, deliminator=delimiter, missingVal=missingVal,
					format=data_format_dict[format], filter=filter, id=id, returnChromosomes=True)
		sd = SNPsDataSet(snpsds, chromosomes, data_format=format)
	else:
		print "Unknown file format!"
		raise Exception()
 	return sd




def parse_snp_data_region(datafile, chromosome, start_pos, end_pos, delimiter=",",
			missingVal='NA', format=1, filter=1, id=None):
	"""
	Return a region of snpsd.
	"""
	snpsds = parseCSVData(datafile, deliminator=delimiter, missingVal=missingVal, format=format, filter=filter, id=id)
	snpsd = SNPsDataSet(snpsds, [1, 2, 3, 4, 5])
	return snpsd.get_region_snpsd(chromosome, start_pos, end_pos)



def parse_chr_pos_list(datafile, delimiter=",", min_marf=0.0):
	"""
	Return a chr_pos list without loading all data...
	"""

	sys.stderr.write("Loading file: %s ... \n" % datafile)

	chr_pos_list = []

	#Reading column data
	f = open(datafile, 'r')
	lines = f.readlines()
	f.close()

	line = lines[1].split(delimiter)
	withArrayIds = line[0] == "Chromosome"
	i = 1
	if withArrayIds:
		i += 1
	while i < len(lines):
		line = lines[i].split(delimiter)
		chr_pos_list.append((int(line[0]), int(line[1])))
		i += 1
	sys.stderr.write("Chromosomes and positions read\n")
	return chr_pos_list



#--------------------------------------------------------------------------------#

#def _testDBParser_():
#
#	snpsDatas250 = parse250DataRaw()
#
#
#	#snpsDatasPerlgen = get2010DataFromDb()
#	snpsDatasPerlgen = getPerlgenDataFromDb()
#	print len(snpsDatasPerlgen)
#	for snpsd in snpsDatasPerlgen:
#		print len(snpsd.positions)
#
#
#	snpsDatas250[0].compareWith(snpsDatasPerlgen[0])
#
#	merged250nPerlgen = []
#	for i in range(0, len(snpsDatas250)):
#		merged250nPerlgen.append(snpsDatasPerlgen[i].mergeData(snpsDatas250[i]))
#
#	writeRawSnpsDatasToFile("Perlgen_250.out", merged250nPerlgen)
#
#	snpsDatas2010 = get2010DataFromDb()
#	print len(snpsDatas2010)
#	for snpsd in snpsDatas2010:
#		print len(snpsd.positions)
#
#
#	snpsDatas2010[0].compareWith(snpsDatasPerlgen[0])
#
#	merged250n2010nPerlgen = []
#	for i in range(0, len(merged250nPerlgen)):
#		merged250n2010nPerlgen.append(snpsDatas2010[i].mergeData(merged250nPerlgen[i]))
#
#	writeRawSnpsDatasToFile("2010_Perlgen_250.out", merged250n2010nPerlgen)
#
#
#	"""  
#	snpErrorRate = []
#	for chr in [0,1,2,3,4]:
#		snpErrorRate +=snpsds2[chr].compareWith(snpsds1[chr])
#	print "Average Genome Wide Snp Error:",sum(snpErrorRate)/float(len(snpErrorRate))
#	"""


def parse_full_snp_data(data_file_prefix, chromosome, data_format='nucleotides'):
	"""
	Handles large dataset (e.g. full SNPs data) in a managable way... even for a 4GB memory mac.
	"""
	data_format_dict = {'binary':0, 'nucleotides':1, 'int':2, 'float':3}
	data_file = data_file_prefix + '_Chr%i.csv' % chromosome
	if data_format == 'nucleotides':
		print 'Looking for nucleotide SNPs'
		(snpsds, chromosomes) = parseCSVData(data_file, deliminator=',', missingVal='N',
					format=data_format_dict[data_format], filter=filter, id=id, returnChromosomes=True)
		sd = SNPsDataSet(snpsds, chromosomes, data_format=data_format)
	elif data_format == 'int': #E.g. binary data.
		pass
	return sd



def parse_1001genomes_data(file_name, convert_to_binary=True, remove_non_binary=True, mac_filter=5, delimiter=',',
		reference_ecotype='6909'):
	"""
	Parse SNPs...
	"""
	if convert_to_binary:
	 	if not remove_non_binary:
			 raise Exception("Can't convert SNPs to binary values when non-binary values are in the data.  Set remove_non_binary option!")

	with open(file_name) as f:
		header = map(str.strip, f.next().split(delimiter))
		print header[0]
		if header[0] != 'Chromosome':
			array_ids = map(int, header[2:])
			header = map(str.strip, f.next().split(delimiter))
		ecotypes = header[2:]
		if convert_to_binary:
			if '6909' in ecotypes:
				ref_i = ecotypes.index('6909')
			else:
				ref_i = 0
				print "Couldn't find reference ecotype!"

		snps = []
		positions = []
		num_mac_filtered = 0
		num_non_binary_filtered = 0
		for i, line in enumerate(f):
			if i % 10000 == 0:
				print i, num_non_binary_filtered, num_mac_filtered, num_non_binary_filtered + num_mac_filtered
			l = map(str.strip, line.split(delimiter))
			chrom = int(l[0])
			pos = int(l[1])
			snp = l[2:]
			nts = set(snp)
			if remove_non_binary:
				if len(nts) != 2:
					num_non_binary_filtered += 1
					continue

			if mac_filter > 0:
				counts = [snp.count(nt) for nt in nts]
				if min(counts) < mac_filter:
					num_mac_filtered += 1
					continue

			if convert_to_binary:
				allele_0 = snp[ref_i]
				f = lambda x: 0 if x == allele_0 else 1
				snp = sp.array(map(f, snp), dtype='int8')

			snps.append(snp)
			positions.append(pos)
	if num_mac_filtered:
		print '%d SNPs were filtered due to small MACs' % num_mac_filtered
	if num_non_binary_filtered:
		print '%d non-binary SNPs were filtered.' % num_non_binary_filtered
	return ecotypes, snps, positions





def _test_full_seq_parser_():
	file_prefix = '/Users/bjarni.vilhjalmsson/Projects/Data/1001genomes/Imputed_2_Chr'
	chromosomes = [4, 1, 2, 3, 5]
	for mac in [0, 15]:
		for chrom in chromosomes:
			snpsd_list = []
			file_name = file_prefix + '%d_mac%d.binary.csv' % (chrom, mac)
			sd = parse_numerical_snp_data(file_name, use_pickle=True)


def load_250K_snps(call_method_id=75, data_format='binary', debug_filter=1.0):
	data_file = '%s250K_t%d.csv' % (env['data_dir'], call_method_id)
	return parse_snp_data(data_file , format=data_format, filter=debug_filter)

def load_1001_full_snps(mac=0, chromosomes=[1, 2, 3, 4, 5], debug_filter=1):
	"""
	Parses the binary SNPs
	"""
	import env
	data_1001_dir = env.env['data_1001_dir']
	print 'Loading full genotype data.'
	t = time.time()
	snpsds = []
	num_snps = 0
	for chrom in chromosomes:
		print 'Loading pickled chromosome %d.' % chrom
		file_name = data_1001_dir + 'Imputed_2_Chr%d_mac%d.binary.csv.pickled' % (chrom, mac)
		with open(file_name) as f:
			sd = cPickle.load(f)
		sd.sample_snps(debug_filter)
		num_snps += len(sd.snpsDataList[0].snps)
		snpsds.append(sd.snpsDataList[0])
	t = time.time() - t
	print 'Loaded %d SNPs.' % num_snps
	print 'It took %d minutes and %0.2f seconds to load the SNPs' % (t / 60, t % 60)
	return SNPsDataSet(snpsds, chromosomes, data_format='binary')


def _bug_test_():
	sd_1001 = load_1001_full_snps()
	sd = parse_numerical_snp_data('/home/GMI/bjarni.vilhjalmsson/Projects/data/250K_t72.csv.binary', use_pickle=True)
	accs1 = set(sd_1001.accessions)
	accs2 = set(sd.accessions)
	print len(accs1)
	print len(accs2)
	print len(accs2 - accs1)
	print accs2 - accs1
	print len(accs1 - accs2)
	print accs1 - accs2

if __name__ == "__main__":

#	snpsds = get2010DataFromDb(host="papaya.usc.edu",chromosomes=[1,2,3,4,5], db = "at", dataVersion="3", user = "bvilhjal",passwd = "bamboo123")
#	print len(snpsds)
#	for i in range(0,len(snpsds)):
#		print len(snpsds[i].snps)
#	sd = parse_snp_data('/Users/bjarni.vilhjalmsson/Projects/Data/250k/250K_t54.csv',look_for_binary=False)
#	sd = parse_binary_snp_data('/Users/bjarni.vilhjalmsson/Projects/Data/250k/250K_t54.csv.binary')
#	import time
#	import random
#	t_ = time.time()
#	f = open('/tmp/test.pickled','wb')
#	cPickle.dump(sd, f, protocol=2)
#	print 'Took %.2f s...'%(time.time()-t_)
#	f.close()

#	f = open('/tmp/test.pickled', 'rb')
#	t_ = time.time()
#	sd = cPickle.load(f)
#	f_accessions = random.sample(sd.accessions, 900)
#	sd.filter_accessions(f_accessions)
#	print 'Took %.2f s...' % (time.time() - t_)
#	f.close()

	#get250KDataFromDb(user="bvilhjal", passwd="bamboo123")

	#pass
	#_testDBParser_()
	#load_1001_full_snps(chromosomes=[1, 2, 3, 4, 5])
	_bug_test_()
	pdb.set_trace()
