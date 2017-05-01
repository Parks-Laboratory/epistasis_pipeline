import os
import argparse

gwas = set()


p = argparse.ArgumentParser()
p.add_argument('-p', action = 'store', default = '.')
args = p.parse_args()
path = args.p
print("path: " + path)
for file in os.listdir(path):
	if file.endswith(".gwas"):
		id = file.split('_')[-1]
		id = id.split('.')[0]
		gwas.add(id)
output = set()
for file in os.listdir('.'):
        if 'output' in file:
 #               print(file)
		id = file.split('.')[-1]
#                print(id)
		output.add(id)
diff = output.difference(gwas)
print((list(diff)))
# write list to file 
with open("diff.txt",'w') as f:
	f.write("\n".join(diff))
print("write numbers in diff.txt")

