#!/usr/bin/env python

'''
The MIT License (MIT)

Copyright (c) 2016 Charles Y. Lin, Rocio Dominguez-Vidana, Barbara Jenko

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

import os
import sys
import string
import time
import utils
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
#projectFolder = '/home/barbara/gen/linlab/data/chordoma/'
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
    def __init__(self,collection_slug,genome,relationship_file=''):


        print('Loading collection %s' % (collection_slug))
        collection = res.collection.get(collection_slug)
        self._collection = collection

        self._id = collection.id
        self._slug = collection_slug
        print("Using genome %s" % (genome))
        self._genome = genome

        sample_dict = {}
        samples = []
        for d in collection.data.filter(type='data:alignment:bam:'):
            sample = d.sample or d.presample  # get sample or presample of data object
            #make sample_dict a nested dictionary
            sample_dict[sample.name]={
                'unique_id': sample.id,  # returns unique id
                'slug': sample.slug,  # returns slug
            }
            samples.append(sample)

        sample_names = sample_dict.keys()
        sample_names.sort()
        self._names = sample_names

        for sample in samples:
            for d in sample.data.filter(type='data:alignment:bam:'):
                sample_dict[sample.name]['bam'] = d.id
            for d in sample.data.filter(type='data:chipseq:macs14:'):
                sample_dict[sample.name]['bed'] = d.id

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

    def id(self):
        '''
        returns the id
        '''
        return self._id

    def slug(self):
        '''
        returns the slug
        '''
        return self._slug

    def update(self):
        '''
        re-loads the collection
        '''
        self._collection = res.collection.get(self._slug)

    def exportRelationships(self,output=''):
        print('Creating relationship table')

        if len(output) > 0:
            print('Writing relationships to %s' % (output))
        rel_table = [['SAMPLE_NAME','SAMPLE_SLUG','U_ID','BACKGROUND_NAME','GROUP']]

        for name in self._names:
            new_line = [name,self._sample_dict[name]['slug'],self._sample_dict[name]['unique_id'],self._background_dict[name],self._group_dict[name]]
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

        f = open(input_table,'r')
        lines = f.readlines()
        for line in lines[1:]:
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
            sys.exit()
        if self._names.count(background_name) == 0:
            print('ERROR: Background name %s not in collection' % (background_name))
            sys.exit()
        self._background_dict[name] = background_name


    def setGroup(self,name,group_name):

        #take in a sample name and a background name and set the background dict

        #check that each exist
        if self._names.count(name) == 0:
            print('ERROR: %s not in collection' % (name))
            sys.exit()
        self._group_dict[name] = group_name

    def names(self):
        #returns all sample names
        return self._names

    def group(self,name):
        #returns all sample names
        return self._group_dict[name]

    def getBamID(self,name):
        return self._sample_dict[name]['bam']

    def getBedID(self,name):
        return self._sample_dict[name]['bed']

    def getBackground(self,name):
        #return the sample name of the background!!!
        background_name = self._background_dict[name]
        if background_name == 'NONE' or background_name == '':
            return None
        else:
            return self._background_dict[name]

    def run_macs14(self, sample_name, useBackground=True, p_value='1e-9', output=''):
        '''
        given a sample and a background name, calculate macs
        '''
        macs_slug = 'macs14' #macs processor slug
        #in order to run this processor we need the slug, the control, treat, genome, p-value

        #get the treat bam id
        treat_id = self.getBamID(sample_name)

        #figuring out genome string
        genome_string_dict = {'HG19':'hs'} #probably should make this dictionary bigger

        genome_string = genome_string_dict[string.upper(self._genome)]

        input_dict = {'t':treat_id,
                      'g':genome_string,
                      'pvalue': p_value,
                      }

        if useBackground:
            background_name = self.getBackground(sample_name)
            if background_name:
                control_id = self.getBamID(background_name)
                input_dict['c'] = control_id
            else:
                print('WARNING: no background dataset found for %s' % (sample_name))
                print('INFO: macs will be run without control')

        macs = res.get_or_run(slug='macs14', input = input_dict)
        print("Calculating macs...")

        while True:
            macs.update()
            if macs.status=='OK':
                break
            elif macs.status=='ER':
                print(macs.stdout())
                print('oh snap')
                sys.exit()

            time.sleep(1)
            # print('Still working ;)')

        self._sample_dict[sample_name]['macs14'] = macs.id

        if len(output) > 0:
            macs.download(download_dir = output)

        return macs

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
    # for data_obj in collection.data:
    #     d = res.data.get(data_obj)
    for d in collection.data.filter(type='data:alignment:bam:'):
        sample = d.sample or d.presample  # get sample or presample form data object
        #make sample_dict a nested dictionary
        sample_dict[sample.name] = {
            'un_ique_id': sample.id,  # returns unique id
            'slug': sample.slug,  # returns slug
        }
    return sample_dict

def get_bam(sample_name,sample_dict):

    '''
    from the sample dict given the sample name, return the aligned bam data object
    '''
    sample=res.sample.get(sample_dict[sample_name]['slug'])
    for d in sample.data.filter(type='data:alignment:bam:'):
        return d.id

def get_bed(sample_name,sample_dict):

     '''
     from the sample dict given the sample name, return the aligned bam data object
     '''
     sample=res.sample.get(sample_dict[sample_name]['slug'])
     for d in sample.data.filter(type='data:chipseq:macs14:'):
        return d.id



# in rose2 we add macs processor and we add  function get_or_run, and if the data object already exist then it take it otherwise it run macs2 with defined parameters.

def run_rose2(res_collection,sample_name, useBackground=True, t=0, output='', s='', macs_params={}):
    '''
    given a sample and a background name, calculate macs
    '''
    # processor check if is macs with this parameters was already run, if not it will run it again
    macs_params.update({
        'res_collection': res_collection,
        'sample_name': sample_name,
        'useBackground': useBackground,
    })

    res_collection = run_macs14(**macs_params)
    # run_macs14 sets sample_dict[sample_name]['bed'] to the last macs
    # object, so we will get the right one
    macs_bed = res_collection.getBedID(sample_name)

    # get the treat bam id
    treat_id = res_collection.getBamID(sample_name)

    #figuring out genome string
    #genome_string_dict = {'HG19':'hs'}

    genome_string = string.upper(res_collection._genome)

    res_collection.getBedID(sample_name)

    input_rose2 = {
        'g': genome_string,
        'i': macs_bed,
        'r': treat_id,
        't': t,
        's': s,
    }

    if useBackground:
        background_name = res_collection.getBackground(sample_name)
        if background_name:
            control_id = res_collection.getBamID(background_name)
            input_rose2['c'] = control_id
        else:
            print('WARNING: no background dataset found for %s' % (sample_name))
            print('INFO: rose-2 will be run without control')

    # Rose2 is automatically added to all collections that sample is in,
    # so we don't need to define collection
    # get_or_run function, checks if macs processor with this inputs has already run
    # if not it runs it otherwise it took the data object
    rose2 = res.get_or_run(slug='rose2', input=input_rose2)
    print("Running rose2...")

    while True:
        rose2.update()
        if rose2.status=='OK':
            break
        elif rose2.status=='ER':
            print(rose2.stdout())
            print('oh snap')
            sys.exit()

        time.sleep(1)

    if len(output) > 0:
        rose2.download(download_dir=output)

    return res_collection


def bamplot(res_collection, g='', gff_path=None, gff_source=None, i_region=None, b=[], stretch_input=None, c=None,
            s=None, e=None, r=None, y=None, n=None, p=None, t=None, scale=None, bed=None,
            multi_page=None, verbose=None, output=''):

    bam_ids = []

    for bam_path in b:
        bam = res.run('bam-upload', input={'src': bam_path})
        bam_ids.append(bam.id)

    genome_string = string.upper(res_collection._genome)


    input_bamplot = {
        'g': genome_string,
        'b': bam_ids,
    }

    if c is not None:
        input_bamplot['c'] = c

    if s is not None:
        input_bamplot['s'] = s

    if e is not None:
        input_bamplot['e'] = e

    if r is not None:
        input_bamplot['r'] = r

    if y is not None:
        input_bamplot['y'] = y

    if n is not None:
        input_bamplot['n'] = n

    if p is not None:
        input_bamplot['p'] = p

    if t is not None:
        input_bamplot['t'] = t

    if scale is not None:
        input_bamplot['scale'] = scale

    if multi_page is not None:
        input_bamplot['multi_page'] = multi_page

    if verbosee is not None:
        input_bamplot['verbose'] = verbose

    if gff_path is not None:
        i_gff = res.run('upload-gff3', input={'src':gff_path, 'source':gff_source})
        input_bamplot['i_gff'] = i_gff

    if i_region is not None:
        input_bamplot['i_region'] = i_region

    if bed is not None:
        input_bamplot['bed'] = bed

    bamplot = res.get_or_run(slug='bamplot', input=input_bamplot)
    print("Running bamplot...")

    while True:
        bamplot.update()
        if bamplot.status=='OK':
            break
        elif bamplot.status=='ER':
            print(bamplot.stdout())
            print('oh snap')
            sys.exit()

        time.sleep(1)

    if len(output) > 0:
        bamplot.download(download_dir=output)

    return res_collection
#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():


    '''
    this is where we are building the main run function for the script
    all of the work should occur here, but no functions should be defined here
    for now it is not executed allowing a check
    '''
    projectFolder = '/grail/projects/pipeline_resdk_2/'
    #projectFolder = '/home/barbara/gen/linlab/data/pipeline_resdk/' #pipeline_resdk needs write permissions

    # ideal situation
    res_collection = ResCollection(collection_slug,'hg19','%sCHORDOMA_TABLE.txt' % (projectFolder)) #if foo exists, load it, if not write it out to disk

    # all of the datasets that we have from the k27ac group
    h3k27ac_list = [name for name in res_collection.names() if res_collection.group(name) == 'H3K27AC']

    #=======================
    #schema for how we want to run, record and get back data for any processor
    sample_name = h3k27ac_list[0]
    print(sample_name)
    #test macs on this and make sure the collection sample dict is appropriately updated
    #this makes sure you can return the macs data id when you run it
    #also want to make sure that the collection gets updated appropriately
    macs = res_collection.run_macs14(sample_name,useBackground=True,p_value='1e-9')

    print(macs.id)
    print(res_collection._sample_dict[sample_name])

    #now we should be able to retrieve the macs output easily by doing
    macs_id = res_collection._sample_dict[sample_name]['macs14']
    macs_copy = res.data.get(macs_id)
    print(macs_copy.id)
    print(macs_copy.files())

    #========================
    # if we want to run only for one sample
    # h3k27ac_list = h3k27ac_list[:1]

    # #run macs on everybody w/ background at p of 1e-9 and download to a folder
    # macs_parent_folder = utils.formatFolder('%smacsFolder' % (projectFolder),True)

    # for sample_name in h3k27ac_list:
    #     #macs_folder = utils.formatFolder('%s%s_MACS14/' % (macs_parent_folder,sample_name),True)
    #     #res_collection = run_macs14(res_collection,sample_name,useBackground=True,p_value='1e-9',output=macs_folder)
    #     res_collection = run_macs14(res_collection,sample_name,useBackground=True,p_value='1e-9')

    for sample_name in h3k27ac_list:
        res_collection = run_rose2(res_collection,sample_name,useBackground=True, t=0, s='', macs_params={'p_value': '1e-9'})


    #retrieve an arbitrary macs output
    #macs_list = res_collection._analysis_dict[sample_name]['macs14']
    #filter for a given p-value or presence/absence of a control
    #macs =res_collection._analysis_dict[sample_name]['macs14'][0]
    #then i could get file paths or object ids necessary to run other stuff


if __name__=="__main__":
        main()
