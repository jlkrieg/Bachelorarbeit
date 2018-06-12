infile  = open("run199datanames.txt","r" )
outfile = open("dataNames.C","w")
outfile.write("void dataNames(){\n  TNtuple* nt = new TNtuple(\"nt\",\"nt\",\"names\");\n")
for line in infile:
    textline=line.rstrip()
    if ".fits" in textline:
        print textline
        str = "  nt->Fill(\""+textline+"\");\n"
        outfile.write(str)
outfile.write("}")
infile.close()
outfile.close()
