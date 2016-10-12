#!/usr/bin/python

'''
The MIT License (MIT)

Copyright (c) 2016 Charles Y. Lin, Rocio Dominguez-Vidana

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys
print "using python version %s" % sys.version



#================================================================================
#=============================LOG INTO RESDK=====================================
#================================================================================

import resdk
res = resdk.Resolwe('admin', 'admin', 'https://torta.bcm.genialis.com')
resdk.start_logging()


#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#for testing
#add locations of files and global parameters in this section
collection_slug = 'primary_chordoma'
genome = 'hg19'
projectFolder = '/grail/projects/chordoma/'


#================================================================================
#===========================DEFINING THE CLASSES=================================
#================================================================================




#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

class ResCollection:

    #this __init__section is called from the get go
    def __init__(self,collection_slug,relationship_file=''):
        
        print('Loading collection %s' % (collection_slug))
        collection = res.collection.get(collection_slug)        
        self._collection = collection

        sample_dict = {}
        for data_obj in collection.data:

            d = res.data.get(data_obj)
            if d.process_type.startswith('data:alignment:bam:'):
                sample = res.sample.filter(data=d.id)[0] #only one bam per id
                #new_line = '{}\t{}\t{}\t{}\t{}\n'.format(sample.name,sample.slug, sample.id, '', '') 
                #data_dict[sample.name] = new_line
                #make sample_dict a nested dictionary
                sample_dict[sample.name]={}
                sample_dict[sample.name]['unique_id']=sample.id #returns unique id
                sample_dict[sample.name]['slug']=sample.slug #returns slug

        sample_names = sample_dict.keys()
        sample_names.sort()
        self._names = sample_names

        for name in sample_names:
            
            sample=res.sample.get(sample_dict[name]['slug'])            
            for data_obj in sample.data:
                d = res.data.get(data_obj)
                if d.process_type.startswith('data:alignment:bam'):

                    sample_dict[name]['bam'] = d.id
                
        self._sample_dict = sample_dict


        #now establish a relationship dictionary for groups and sample relationships
        self._background_dict = {}
        for name in sample_names:
            self._background_dict[name]='NONE'

        #now establish a relationship dictionary for groups and sample relationships
        self._group_dict = {}
        for name in sample_names:
            self._group_dict[name]='NONE'

        #if a relationship table is provided, import it.
        #if not, create a basic table
        
        if len(relationship_file) > 0: # something provided by user
            #first step is to see if it exists
            try:
                f = open(relationship_file,'r')
                f.close()
                self.importRelationships(relationship_file)
            except IOError: #nothing present
                self.exportRelationships(relationship_file)

    def exportRelationships(self,output=''):
        print('Creating relationship table')
        if len(output) >0:
            print('Writing relationships to %s' % (output))
        rel_table = [['SAMPLE_NAME','SAMPLE_SLUG','U_ID','BACKGROUND_NAME','GROUP']]

        for name in self._names:            
            new_line = [name,self._sample_dict[name]['slug'],self._sample_dict[name]['unique_id'],self._background_dict[name],'']
            rel_table.append(new_line)

        if len(output) == 0:
            return rel_table
        else:
            outfile = open(output,'w')
            for line in rel_table:
                outfile.write('\t'.join([str(x) for x in line])+'\n')

            outfile.close()


    def importRelationships(self,input_table):
        print('Importing relationship data from %s' % (input_table))

        rel_table = [['SAMPLE_NAME','SAMPLE_SLUG','U_ID','BACKGROUND_NAME','GROUP']]
        
        f = open(input_table,'r')
        lines = f.readlines()
        for line in lines:
            line = line.rstrip().split('\t')
            name = line[0]
            background = line[3]
            group = line[4]
            self._background_dict[name]=background
            self._group_dict[name]=group


    def setBackground(self,name,background_name):

        #take in a sample name and a background name and set the background dict
        
        #check that each exist
        if self._names.count(name) == 0:
            print('ERROR: %s not in collection' % (name))
        if self._names.count(background_name) == 0:
            print('ERROR: %s not in collection' % (background_name))

        self._background_dict[name] = background_name


    def setGroup(self,name,group_name):

        #take in a sample name and a background name and set the background dict
        
        #check that each exist
        if self._names.count(name) == 0:
            print('ERROR: %s not in collection' % (name))

        self._group_dict[name] = group_name

    def names(self):
        #returns all sample names
        return self._names

    def group(self,name):
        #returns all sample names
        return self._group_dict[name]


    def getBamID(self,name):

        return self._sample_dict[name]['bam']


    def getBackground(self,name):
        #return the sample name of the background!!!
        background_name = self._background_dict[name]
        if background_name == 'NONE' or background_name == '':
            return None
        else:
            return self._background_dict[name]


#class locus:
    #locus class once i get macs going
#    def __init__():




#================================================================================
#=================================FUNCTIONS======================================
#================================================================================

#retrieve information from using the slug



def make_sample_dict(collection):

    '''
    dictionary keyed by sample name that has some useful pointers
    '''
    sample_dict = {}
    for data_obj in collection.data:
        d = res.data.get(data_obj)
        if d.process_type.startswith('data:alignment:bam:'):
            sample = res.sample.filter(data=d.id)[0] #only one bam per id
            #new_line = '{}\t{}\t{}\t{}\t{}\n'.format(sample.name,sample.slug, sample.id, '', '') 
            #data_dict[sample.name] = new_line
            #make sample_dict a nested dictionary
	    sample_dict[sample.name]={}
	    sample_dict[sample.name]['unique_id']=sample.id #returns unique id
	    sample_dict[sample.name]['slug']=sample.slug #returns slug
    return sample_dict

def get_bam(sample_name,sample_dict):

    '''
    from the sample dict given the sample name, return the aligned bam data object
    '''
    sample=res.sample.get(sample_dict[sample_name]['slug'])
    for data_obj in sample.data:
        d = res.data.get(data_obj)
        if d.process_type.startswith('data:alignment:bam'):
            return d.id

def calculate_macs(case, background):
    '''
    given a sample and a background name, calculate macs
    '''
    #macs = res.run("macs", inputs=(case = case.id, control = background.id))
    macs = res.run("macs", inputs=(case = case.id))
    return macs

#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():

    '''
    this is the main run function for the script
    all of the work should occur here, but no functions should be defined here
    '''
    projectFolder = '/grail/projects/pipeline_resdk/'
    
    #ideal situation
    res_collection = ResCollection(collection_slug,'%sfoo.txt' % (projectFolder)) #if foo exists, load it, if not write it out to disk

    #all of the datasets that we have
    #names_list = res_collection.names()
    
    #only want k27ac datasets
    names_list = [name for name in res_collection.names() if res_collection.group(name) == 'H3K27AC']
    for name in names_list:
        background_name = res_collection.getBackground(name)
        if background_name:
            bamID = res_collection.getBamID(background_name)
            print('For dataset %s, the background is %s and the bamID for the background is %s' % (name,background_name,bamID))
        else:
            print('For dataset %s, No background was found' % (name))

    #run all MACS
    #using the collection list above...
    for case_name in names_list: 
        background_name = res_collection.getBackground(name)
        if background_name:
            caseID = res_collection.getBamID(case_name)
            backgroundID = res_collection.getBamID(background_name)
            macs_result = calculate_macs(case,bamID) #sending funny error of not access permission

        



main()

