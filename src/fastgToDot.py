

fIn = open("/ebio/abt6_projects7/small_projects/mdubarry/Documents/pbsim/before_rr.fastg","r")
fOut = open("/ebio/abt6_projects7/small_projects/mdubarry/Documents/pbsim/before_rr.dot","w")

fOut.write("digraph adj {\n")
for line in fIn:
	if line[0] == ">":
		line = line.rstrip() #delete end line "/n"
		line =  line[1:-1] #delete charater ">" and last charater ";"
		node1 = line.split(":")[0]
		if node1.find("'") != -1:
			fOut.write(node1[:-1]+"-")
		else:
			fOut.write(node1+"+")
		if line.find(":") != -1:
			node2 = line.split(":")[1]
			fOut.write(" -> ")
			if node2.find("'") != -1:
				fOut.write(node2[:-1]+"- ;\n")
			else:
				fOut.write(node2+"+ ;\n")
		else:
			fOut.write(" ;\n")
			
fOut.write("}")
fIn.close()
fOut.close()
