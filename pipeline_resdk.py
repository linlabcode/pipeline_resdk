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
import time
print("using python version %s" % sys.version)


#================================================================================
#=============================LOG INTO RESDK=====================================
#================================================================================

import resdk
res = resdk.Resolwe('admin', 'admin', 'https://torta.bcm.genialis.com')
resdk.start_logging()


#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

# for testing
# add locations of files and global parameters in this section
collection_slug = 'primary_chordoma'
genome = 'hg19'
# projectFolder = '/home/barbara/gen/linlab/data/chordoma/'
projectFolder = '/grail/projects/chordoma/'


#================================================================================
#===========================DEFINING THE CLASSES=================================
#================================================================================




#================================================================================
#===================================CLASSES======================================
#================================================================================



class ResCollection(object):

    #: list of processing objects that need to be downloaded
    to_download = []

    # this __init__section is called from the get go
    def __init__(self, collection_slug, genome, relationship_file=None):
        print('Loading collection %s' % (collection_slug))
        collection = res.collection.get(collection_slug)
        self._collection = collection
        self._id = collection.id
        self._slug = collection_slug

        print("Using genome %s" % (genome))
        self._genome = genome

        sample_dict = {}
        for data in collection.data.filter(type='data:alignment:bam:'):
            sample = data.sample or data.presample  # get sample or presample of data object
            sample_dict[sample.name] = {
                'sample': sample,
                'unique_id': sample.id,  # returns unique id
                'slug': sample.slug,  # returns slug
                'background': 'NONE',
                'group': 'NONE',
            }
        self._sample_dict = sample_dict

        if relationship_file:  # something provided by user
            if os.path.exists(relationship_file):
                self.importRelationships(relationship_file)
            else:
                self.exportRelationships(relationship_file)

    def importRelationships(self, input_table):
        print('Importing relationship data from %s' % (input_table))

        with open(input_table, 'r') as fn:
            fn.readline()  # skip the header
            for line in fn.readlines():
                line = line.rstrip().split('\t')
                name = line[0]
                background = line[3]
                group = line[4]
                self._sample_dict[name]['background'] = background
                self._sample_dict[name]['group'] = group

    def exportRelationships(self, output=None):
        print('Creating relationship table')

        if output:
            print('Writing relationships to %s' % (output))

        rel_table = [['SAMPLE_NAME', 'SAMPLE_SLUG', 'U_ID', 'BACKGROUND_NAME', 'GROUP']]
        for name in self._sample_dict:
            sample_data = self._sample_dict[name]
            rel_table.append([
                name,
                sample_data['slug'],
                sample_data['unique_id'],
                sample_data['background'],
                sample_data['group'],
            ])

        if not output:
            return rel_table
        else:
            with open(output, 'w') as outfile:
                for line in rel_table:
                    outfile.write('\t'.join([str(x) for x in line]) + '\n')

    def names(self):
        # returns all sample names
        return sorted(self._sample_dict.keys())

    def getGroup(self,name):
        # returns all sample names
        return self._sample_dict[name]['group']

    def getBackground(self, name):
        # return the sample name of the background!!!
        background_name = self._sample_dict[name]['background']
        if background_name == 'NONE' or background_name == '':
            return None
        else:
            return background_name

    def getMacs(self, name):
        sample = self._sample_dict[name]['sample']
        return sample.get_macs()

    def getBam(self, name):
        sample = self._sample_dict[name]['sample']
        return sample.get_bam()

    def download(self, output=''):
        print("Waiting for analysis to finish...")
        while True:
            new_list = []
            for obj in self.to_download:
                obj.update()
                if obj.status == 'OK':
                    print('Downloading data for: {}'.format(obj.name))
                    obj.download(download_dir=output)
                elif obj.status == 'ER':
                    print('Error in object: {}'.format(obj.name))
                    print(obj.stdout())
                    print('----------')
                else:
                    new_list.append(obj)

            self.to_download = new_list

            if new_list == []:
                break

            time.sleep(1)

    def run_macs(self, sample_name, useBackground=True, p_value='1e-9', watch=False):
        sample = self._sample_dict[sample_name]['sample']

        if useBackground:
            background_name = self.getBackground(sample_name)
            if background_name:
                background_slug = self._collection.samples.get(name=background_name).slug
            else:
                print('WARNING: no background dataset found for %s' % (sample_name))
                print('INFO: macs will be run without control')
                background_slug = None
        else:
            background_slug = None

        macs = sample.run_macs(
            use_background=useBackground,
            background_slug=background_slug,
            p_value=p_value
        )
        if watch:
            self.to_download.append(macs)

        return macs

    def run_rose2(self, sample_name, useBackground=True, tss=0, stitch=None, watch=False, macs_params={}):
        sample = self._sample_dict[sample_name]['sample']

        if useBackground:
            background_name = self.getBackground(sample_name)
            if background_name:
                background_slug = self._collection.samples.get(name=background_name).slug
            else:
                print('WARNING: no background dataset found for %s' % (sample_name))
                print('INFO: rose2 will be run without control')
                background_slug = None
        else:
            background_slug = None

        genome_string = self._genome.upper()

        # get existing macs or create a new one
        macs = sample.run_macs(
            use_background=useBackground,
            background_slug=background_slug,
            **macs_params
        )

        rose = sample.run_rose2(
            use_background=useBackground,
            background_slug=background_slug,
            genome=genome_string,
            tss=tss,
            stitch=stitch,
            beds=macs  # run rose only on one macs peaks bed file
        )
        if watch:
            self.to_download.extend(rose)  # run_rose returns list

        return rose[0]  # we only run one rose at the time

    def run_bamplot(self, sample_names, input_gff=None, input_region=None, stretch_input=None,
                    color=None, sense=None, extension=None, rpm=None, yscale=None, names=None,
                    plot=None, title=None, scale=None, bed=None, multi_page=None, watch=False):

        bams = [self.getBam(name) for name in sample_names]
        genome_string = self._genome.upper()

        bamplot = res.run_bamplot(bams, genome_string, input_gff, input_region, stretch_input,
                color, sense, extension, rpm, yscale, names, plot, title, scale, bed, multi_page)

        if watch:
            self.to_download.append(bamplot)

        self._collection.add_data(bamplot)

        return bamplot


def main():


    '''
    this is where we are building the main run function for the script
    all of the work should occur here, but no functions should be defined here
    for now it is not executed allowing a check
    '''
    projectFolder = '/grail/projects/pipeline_resdk/'
    #projectFolder = '/home/barbara/gen/linlab/data/pipeline_resdk/' #pipeline_resdk needs write permissions

    # ideal situation
    res_collection = ResCollection(collection_slug,'hg19','%sCHORDOMA_TABLE.txt' % (projectFolder)) #if foo exists, load it, if not write it out to disk

    # all of the datasets that we have from the k27ac group
    h3k27ac_list = [name for name in res_collection.names() if res_collection.getGroup(name) == 'H3K27AC']

    #=======================
    #schema for how we want to run, record and get back data for any processor
    sample_name = h3k27ac_list[0]
    print(sample_name)
    #test macs on this and make sure the collection sample dict is appropriately updated
    #this makes sure you can return the macs data id when you run it
    #also want to make sure that the collection gets updated appropriately
    macs = res_collection.run_macs(sample_name,useBackground=True,p_value='1e-9', watch=True)
    rose = res_collection.run_rose2(sample_name, useBackground=True, tss=0, stitch=None, macs_params={'p_value': '1e-9'}, watch=True)
    res_collection.download(output='/grail/genialis/pipeline_resdk')

    #gff = res.run('upload-gtf', input={'src':'<path/to/gff'})  # upload gff file
    #gff = res.data.get()  # get gff file onece it is uploaded
    #bed = res.run('upload-bed', input={'src':'<path/to/bed>'})
    gff_region = 'chr17:+:41468594-41566948'
    bamplot = res_collection.run_bamplot(sample_names=h3k27ac_list, input_region=gff_region, watch=True)
    res_collection.download(output='/grail/genialis/pipeline_resdk')

    # if we would like to use bams that are not part of the res_collection
    #bam = res.run('upload-bam', input='src':'<path/to/bam>')
    #bam1 = res.run('upload-bam', input='src':'<path/to/bam1>')
    #bamplot = res.run_bamplot(bam=[bam.id, bam1.id], genome='',... )

    print(macs.id)
    print(rose.id)
    print(res_collection._sample_dict[sample_name])

    #now we should be able to retrieve the macs output easily by doing

    print(macs.files())

    #========================
    # #run macs on everybody w/ background at p of 1e-9 and download to a folder
    # macs_parent_folder = utils.formatFolder('%smacsFolder' % (projectFolder),True)

    # for sample_name in h3k27ac_list:
    #     #macs_folder = utils.formatFolder('%s%s_MACS14/' % (macs_parent_folder,sample_name),True)
    #     #res_collection = run_macs14(res_collection,sample_name,useBackground=True,p_value='1e-9',output=macs_folder)
    #     res_collection = run_macs14(res_collection,sample_name,useBackground=True,p_value='1e-9')

    for sample_name in h3k27ac_list:
        res_collection.run_rose2(sample_name, useBackground=True, tss=0, stitch=None, macs_params={'p_value': '1e-9'}, watch=True)
        res_collection.run_bamplot(sample_names=h3k27ac_list, input_region=gff_region, watch=True, title='h3k27ac_list')
    res_collection.download(output='/grail/genialis/pipeline_resdk')


    #retrieve an arbitrary macs output
    #macs_list = res_collection.get;acs(sample_name)
    #then i could get file paths or object ids necessary to run other stuff


if __name__=="__main__":
        main()
