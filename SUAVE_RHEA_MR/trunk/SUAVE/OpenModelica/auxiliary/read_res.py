"""
	Wenn read_mat eine Methode einer Klasse ist, in __init__
		self.state = 0
	hinzufuegen.
	
	Der bottle neck war die whileSchleife. Die braucht man jedoch nur beim ersten Mal.
	Deshalb wird sie ueberprungen, wenn self.state = 1 ist. Das Auslesen ist nun viel schneller.
"""

def read_mat(self, filename, filepath=''):
	''' 
	Auswertungsfunktion für *.mat-Files, die mit OMPython generiert wurden
	
	Eingangsvariablen:
		filename := Name des *.mat-File (String)
		filepath := Pfad zum *.mat-File (String) 
	Ausgangsvariablen:
		res := Dictonary, key := Variablenname des Modells (Parameter bekommen
															nur zwei Werte 
															übergeben)
	
	'''
	###########################################################################
	# Einlesen der *.mat-Datei
	if len(filepath)>0 and filepath[-1] != '/':
		filename='/'+filename
	try:
		mat = scipy.io.loadmat(filepath+filename)
	except:
		print('This file does not exsist: ', filepath+filename)
		
	###########################################################################
	# Auswerten der *.mat-Datei
	num_var = len(mat['name'][0])
	max_len = len(mat['name'])
	res = {}
	
	if self.state == 0:
		for k in range(num_var):
			i = 0
			# Einlesen der Variablennamen
			name = ''
			while i < max_len-1:
				try:
					if mat['name'][i][k] != '\x00':
						name +=  mat['name'][i][k]
						i += 1
					else:
						break
				except:
					break
			else:
				pass
			self.name.append(name)
		self.state = 1
	else:
		pass
		
	for n in range(num_var):
		# Parameter:= data_1 / Variable := data_2
		if mat['dataInfo'][0][n] == 0:
			dataSet = 'data_2'
		else:
			dataSet = 'data_' + str(mat['dataInfo'][0][n])
			
		# Identifizieren von Alias
		pointer = mat['dataInfo'][1][n] 
		if pointer < 0:
			pointer *= -1
			pre = -1
		else:
			pre = 1
		res.update({self.name[n]: self.get_res(mat, pre, dataSet, pointer)})
	
	return res    	