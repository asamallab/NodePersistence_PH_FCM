# This is a Jython script calculated the persistent homology of a clique dictionary. 
# It has been developed to work with a weighted clique filtration, but it 
# can adapted to be used to any series of persistent homology dictionaries
#
# NOTE: THIS CODE RUNS IN JYTHON, NOT PYTHON!!! 

import pickle, sys, os 
from collections import defaultdict
import Holes 

def list2simplexes(list,dim):
	num=dim+1;
	simplexes=[];
	for i in range(0,len(list),num):
		simplexes.append(list[i:i+num]);
	return simplexes
			
if len(sys.argv)>=4:
	clique_dic_file=str(sys.argv[1]);
	dimension=int(sys.argv[2]);                                               
	Dir=str(sys.argv[3]);
	stringie=str(sys.argv[4]);
	javaplex_path=str(sys.argv[5]);
	save_generators = bool(sys.argv[6])
else:
	print('This code needs as input:');
	print('1) full filtration file name');
	print('2) the maximum homology dimension to calculate');
	print('3) the directory name for output');
	print('4) the tag name for the output files');
	print('6) the full path to your javaplex directory');
	sys.exit();

print('Opening filtration file...');
file=open(clique_dic_file,'rb');
Clique_dictionary=pickle.load(file);                                                                                        


# My modification of the code!
## NOTE: you need to put here the path to the javaPlex distribution on your system
libs = [                                                                                                 
	 os.path.join(javaplex_path,'javaplex.jar')
 	]                    
print(libs)                                                  

for s in libs:                                                                                           
	sys.path.append(s)                                                                                      
print(sys.path)

# I just add it brutally - this is for my local machine  
if os.path.exists(' C:/Users/hp/Documents/GitHub/NodePersistence_PH_FCM/Code/Holes/jython_utils/lib/javaplex.jar'):
   sys.path.append('C:/Users/hp/Documents/GitHub/NodePersistence_PH_FCM/Code/Holes/jython_utils/lib/javaplex.jar')
    


import edu.stanford.math.plex4                                                                              
import edu.stanford.math.plex4.api
#import edu.stanford.math.plex.Persistence
Complex=edu.stanford.math.plex4.api.Plex4.createExplicitSimplexStream();
from java.io import BufferedWriter, FileWriter


max_index=0;
print("Clique dictionary: parsing started.");
for key in Clique_dictionary:
	original_key=key;
	key=str(key);    
	key=key.strip('[]');
	key=key.split(', ');
	key_buona=[];
	for n in range(len(key)):
		key_buona.append(int(float(eval(key[n]))));
	if len(key_buona)==1:
		Complex.addVertex(key_buona[0],0);
	else:	
		Complex.addElement(key_buona, int(Clique_dictionary[original_key][0]));
		if int(Clique_dictionary[original_key][0])>max_index:
			max_index=int(Clique_dictionary[original_key][0]);

print("Parsing over. Closing now.")
Complex.finalizeStream();
print("Complex is valid? ", Complex.validateVerbose());
print("Size of complex filtration:" , Complex.getSize());
max_filtration_value=max_index;
pH=edu.stanford.math.plex4.api.Plex4.getModularSimplicialAlgorithm(dimension+1,2);
print("Starting pH calculation...")
complex_computation=pH.computeIntervals(Complex);
print("Done!")
print("Results incoming:")
infinite_barcodes = complex_computation.getInfiniteIntervals()
annotated_intervals = pH.computeAnnotatedIntervals(Complex);
betti_numbers_string = infinite_barcodes.getBettiNumbers()
print('The betti numbers are:', betti_numbers_string);
print('while the annotated intervals are: \n', annotated_intervals);

# betti_file_name = Dir + 'Betti_number_'+str(stringie)+'.txt'
# output_file = betti_file_name
# writer1 = BufferedWriter(FileWriter(output_file))
# writer1.write(betti_numbers_string)
# print("Annotated intervals saved to", output_file)
# writer1.close()


# generators = annotated_intervals.getGeneratorsAtDimension(1)
# intervals_file = Dir + 'intervals_1d_'+str(stringie)+'.txt'
# try:
    # writer = BufferedWriter(FileWriter(intervals_file))
    # writer.write(str(annotated_intervals.getIntervalsAtDimension(1)))
# #    print("Persistence intervals saved to", intervals_file)
# finally:
    # writer.close();

# # Save generators to file
# generators_file = Dir + 'generators_1d_'+str(stringie)+'.txt'
# try:
    # writer = BufferedWriter(FileWriter(generators_file))
    # for gen in generators:
        # writer.write(str(gen) + "\n")
        
# finally:
    # writer.close()
    

# Here we save the full generator dictionary and save the interval files in order to be
# able to reopen them later for other purposes, for example comparison of random and null 
# models.. 


gendir=Dir+'gen'
if not os.path.exists(gendir):
    os.makedirs(gendir)

Generator_dictionary={};
import re
import string
for h in range(dimension+1):
	Generator_dictionary[h]=[];
	list_gen=list(annotated_intervals.getGeneratorsAtDimension(h))
    #list_gen_strings = [str(gen) for gen in list_gen]
	#if save_generators == True:
        #genfilename=gendir + '/details_generators_'+str(h)+'_'+str(stringie)+'.txt';
        #writer = BufferedWriter(FileWriter(genfilename));
        #writer.close();
    
    
	list_intervals=list(annotated_intervals.getIntervalsAtDimension(h))
	for n,key in enumerate(list_gen):
		test=str(list_intervals[n]).split(',');
		test[0]=test[0].strip(' [')
		test[1]=test[1].strip(' )')
		if test[1]=='infinity':
			test[1]=str(max_filtration_value);
		line=str(key);
        #writer.write(line + "\n");
		line = line.translate(string.maketrans('', ''), '-[]')
		line = re.sub('[,+ ]', ' ', line)
		line=line.split();
		tempcycle=Holes.Cycle(h,list2simplexes(line,h),test[0],test[1]);
		Generator_dictionary[h].append(tempcycle);
		del tempcycle;	
	for cycle in Generator_dictionary[h]:
		cycle.summary();

filename=os.path.join(gendir,'generators_'+str(stringie)+'.pck')
generator_dict_file=open(filename,'wb');
pickle.dump(Generator_dictionary,generator_dict_file);
print('Generator dictionary dumped to '+filename)
generator_dict_file.close();



	
	
