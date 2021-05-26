import sys

if(len(sys.argv) < 2): #if no file location is provided then show error message and exit
    print('Provide File Location')
    sys.exit()

fil = sys.argv[1]
csvfilename = open(fil, 'r').readlines()
#store header values
header = csvfilename[0] 
#remove header from list
csvfilename.pop(0) 
file = 1
#Number of lines to be written in new file
record_per_file = 50

for j in range(len(csvfilename)):
    if j % record_per_file == 0:
	write_file = csvfilename[j:j+record_per_file]
        #adding header at the start of the write_file
        write_file.insert(0, header)
 	#write in file
 	open(str(fil)+ str(file) + '.csv', 'w+').writelines(write_file)
        file += 1
