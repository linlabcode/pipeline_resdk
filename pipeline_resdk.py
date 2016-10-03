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

collection = res.collection.get(collection_slug)


#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes her

class DataDict:

    #this __init__section is called from the get go
    def __init__(self,collection,relationship_file=''):
        print('Loading collection %s' % (collection.slug))
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
            self._background_dict[name]=''

        
    def exportRelationships(self,output=''):
        
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
        background_dict={}
        rel_table = [['SAMPLE_NAME','SAMPLE_SLUG','U_ID','BACKGROUND_NAME','GROUP']]
        filename=input_table
        with open(input_table) as f:
            for line in f:
                sample.name,sample.slug,sample.id,background.name,group=f.split("\t")
                background_dict[sample.name]={}
                background_dict[sample.name]['unique_id']=sample.id #returns unique id
                background_dict[sample.name]['slug']=sample.slug #returns slug
                background_dict[sample.name]['background']=background.name #returns background name
                background_dict[sample.name]['group']=group #returns group



    def setBackground(self,name,background_name):

        #take in a sample name and a background name and set the background dict
        
        #check that each exist
        if self._names.count(name) == 0:
            print('ERROR: %s not in collection' % (name))
        if self._names.count(background_name) == 0:
            print('ERROR: %s not in collection' % (background_name))

        self._background_dict[name] = background_name

    def names(self):
        #returns all sample names
        return self._names

    def getBamID(self,name):

        return self._sample_dict[name]['bam']


    def getBackground(self,name):
        #return the sample name of the background!!!
        return self._background_dict[name]









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


#macs = res.run("macs", inputs=(case=d.id))

#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():

    '''
    this is the main run function for the script
    all of the work should occur here, but no functions should be defined here
    '''

    # sample_dict = make_sample_dict(collection)

    # print(sample_dict)

    # bam_id = get_bam('PRIMARY_CHOR_142A2_H3K27AC',sample_dict)

    # print(bam_id)


    #my_name=sample_dict.keys()[1]
    
    #my_bam=get_bam(my_name,sample_dict)
    
    #load the dataDict object
    dataDict = DataDict(collection)

    dataDict.exportRelationships('/home/rociod/src/foo.txt')
    dataDict.importRelationships('/home/rociod/src/foo.2.txt')


    # print(dataDict.names())

    # print(dataDict.getBam('PRIMARY_CHOR_142A2A4_H3K27AC'))

    
main()

