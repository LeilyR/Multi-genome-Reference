

fIn = open("/ebio/abt6_projects7/small_projects/mdubarry/Documents/data/before_rr.fastg","r")
fOut = open("/ebio/abt6_projects7/small_projects/mdubarry/Documents/data/before_rr.dot","w")

fOut.write("digraph adj {\n")
for line in fIn:
	if line[0] == ">":
		line = line.rstrip() #delete end line "/n"
		line =  line[1:-1] #delete charater ">" and last charater ";"
		node1 = line.split(":")[0]
		if node1.find("'") != -1:
			newNode1 = node1[:-1]+"-"
			#fOut.write(newNode1)
		else:
			newNode1 = node1+"+"
			#fOut.write(newNode1)

		if line.find(":") != -1:
			nodes = line.split(":")[1]
			#fOut.write(" -> ")
			#if nodes.find(",") != -1: #more than one node
			for n in nodes.split(","):
				fOut.write(newNode1+" -> ")
				if n.find("'") != -1:
					fOut.write(n[:-1]+"- ;\n")
				else:
					fOut.write(n+"+ ;\n")
		else:
			fOut.write(newNode1+" ;\n")

			
fOut.write("}")
fIn.close()
fOut.close()
