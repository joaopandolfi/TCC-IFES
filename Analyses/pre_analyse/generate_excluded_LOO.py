import pandas as pd
#leitura de arquivo
def read_arq(file_in):
	baseFile = "user/input/SZ_KATO/seeds/seeds_schizophrenia.txt"
	conteudo =""
	lista= []
	i = 0 

	s0 = pd.read_table(baseFile, header=None)
	arq = open('./'+file_in , 'rt')
	conteudo = arq.readline()
	while conteudo != '':
		if(conteudo[-1] =='\n'):
			conteudo = conteudo[:-1]
			
		path = "user/input/"+conteudo+"/seeds/seeds_schizophrenia.txt"
		s1 = pd.read_table(path, header=None)

		m = set(s0[1]) - set(s1[1])

		elem = m.pop()
		elem = s0.ix[s0[1] == elem]
		#print(elem.values[0])
	
		path +="_EXCLUDED.txt"
		arq2 = open('./'+path , 'w')
		print(path)
		arq2.write(str(elem.values[0][0])+"\t"+elem.values[0][1]+"\t"+elem.values[0][2])
		arq2.close()

		conteudo = arq.readline()
	arq.close()
	return lista

read_arq("mapLooIn.txt")
