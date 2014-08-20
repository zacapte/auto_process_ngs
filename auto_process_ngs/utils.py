#!/bin/env python
#
#     auto_process_utils: utility classes & funcs for auto_process module
#     Copyright (C) University of Manchester 2013 Peter Briggs
#
########################################################################
#
# auto_process_utils.py
#
#########################################################################

__version__ = "0.0.16"

"""auto_process_utils

Utility classes and functions to support auto_process module.

Ultimately these should be relocated in the main 'genomics' code
tree at some point.

"""

#######################################################################
# Imports
#######################################################################

import os
import logging
import IlluminaData
import TabFile
import JobRunner
import Pipeline
import qcreporter
import bcf_utils
from ConfigParser import ConfigParser,NoOptionError

#######################################################################
# Classes
#######################################################################

class AutoProcessConfigParser(ConfigParser):
    """Wraps ConfigParser to set defaults for missing options

    Implements a wrapper for ConfigParser:

    - 'get' and 'getint' methods take a 'default' argument, which
      is returned if the specified option is missing from the file
    - implements a 'getrunner' method that returns a JobRunner
      instance based on a specification string.

    """
    def __init__(self):
        ConfigParser.__init__(self)
    def get(self,section,option,default=None):
        try:
            return ConfigParser.get(self,section,option)
        except NoOptionError:
            return default
    def getint(self,section,option,default):
        try:
            return ConfigParser.getint(self,section,option)
        except NoOptionError:
            return default
    def getrunner(self,section,option,default):
        return fetch_runner(self.get(section,option,default))

class AnalysisFastq:
    """Class for extracting information about Fastq files

    Given the name of a Fastq file, extract data about the sample name,
    barcode sequence, lane number, read number and set number.

    The name format can be a 'full' Fastq name as generated by CASAVA,
    which follows the general form:

    <sample_name>_<barcode_sequence>_L<lane_number>_R<read_number>_<set_number>.fastq.gz

    e.g. for

    NA10831_ATCACG_L002_R1_001.fastq.gz

    sample_name = 'NA10831_ATCACG_L002_R1_001'
    barcode_sequence = 'ATCACG'
    lane_number = 2
    read_number = 1
    set_number = 1

    Alternatively it can be a 'reduced' version where one or more
    of the components has been omitted (typically because they are
    redundant in uniquely distinguishing a Fastq file within a
    set of Fastqs).

    The reduced formats are:

    <sample_name>
    <sample_name>_L<lane_number>
    <sample_name>_<barcode_sequence>
    <sample_name>_<barcode_sequence>_L<lane_number>

    with an optional suffix '_R<read_number>' for paired end sets.

    e.g.

    NA10831
    NA10831_L002
    NA10831_ATCACG
    NA10831_ATCACG_L002

    Provides the follow attributes:

    fastq:            the original fastq file name
    sample_name:      name of the sample (leading part of the name)
    barcode_sequence: barcode sequence (string or None)
    lane_number:      integer (or None if no lane number)
    read_number:      integer (or None if no read number)
    set_number:       integer (or None if no set number)

    """

    def __init__(self,fastq):
        """Create and populate a new AnalysisFastq object

        Arguments:
          fastq: name of the fastq.gz (optionally can include leading path)

        """
        # Store name
        self.fastq = fastq
        # Values derived from the name
        self.sample_name = None
        self.barcode_sequence = None
        self.lane_number = None
        self.read_number = None
        self.set_number = None
        # Base name for sample (no leading path or extension)
        fastq_base = os.path.basename(fastq)
        try:
            i = fastq_base.index('.')
            fastq_base = fastq_base[:i]
        except ValueError:
            pass
        # Identify which part of the name is which
        #
        # Full Illumina-style names are e.g.
        # NH1_ChIP-seq_Gli1_ACAGTG_L001_R1_001
        #
        # We have shorter name formats where redundant parts are
        # omitted, the patterns are:
        # NAME          e.g. NH1_ChIP-seq_Gli1
        # NAME+LANE     e.g. NH1_ChIP-seq_Gli1_L001
        # NAME+TAG      e.g. NH1_ChIP-seq_Gli1_ACAGTG
        # NAME+TAG+LANE e.g. NH1_ChIP-seq_Gli1_ACAGTG_L001
        #
        # Also read number (i.e. R1 or R2) is appended but only for
        # paired end samples
        #
        # The set number is never included, except for full names
        fields = fastq_base.split('_')
        # Deal with set number first e.g. 001
        field = fields[-1]
        ##logging.debug("Test for set number %s" % field)
        if len(field) == 3 and field.isdigit():
            self.set_number = int(field)
            fields = fields[:-1]
        # Deal with trailing read number e.g. R1
        field = fields[-1]
        ##logging.debug("Test for read number %s" % field)
        if len(field) == 2 and field.startswith('R'):
            self.read_number = int(field[1])
            fields = fields[:-1]
        # Deal with trailing lane number e.g. L001
        field = fields[-1]
        ##logging.debug("Test for lane number %s" % field)
        if len(field) == 4 and field.startswith('L') and field[1:].isdigit():
            self.lane_number = int(field[1:])
            fields = fields[:-1]
        # Deal with trailing index tag e.g. ATTGCT or ATTGCT-CCTAAG
        field = fields[-1]
        ##logging.debug("Test for barcode sequence %s" % field)
        if len(fields) > 1:
            # This mustn't be the last field: if it is then it's
            # not the tag - it's the name
            is_tag = True
            for f in field.split('-'):
                for c in f:
                    is_tag = is_tag and c in 'ACGTN'
            if is_tag:
                self.barcode_sequence = field
                fields = fields[:-1]
                ##logging.debug("Identified barcode sequence as %s" % self.barcode_sequence)
        # What's left is the name
        ##logging.debug("Remaining fields: %s" % fields)
        self.sample_name = '_'.join(fields)
        assert(self.sample_name != '')

    def __repr__(self):
        """Implement __repr__ built-in

        """
        if self.set_number is not None:
            # Return the full name
            fq = "%s_%s_L%03d_R%d_%03d" % (self.sample_name,
                                           'NoIndex' \
                                           if self.barcode_sequence is None \
                                           else self.barcode_sequence,
                                           self.lane_number,
                                           self.read_number,
                                           self.set_number)
        else:
            # Reconstruct partial name
            fq = "%s" % self.sample_name
            if self.barcode_sequence is not None:
                fq = "%s_%s" % (fq,self.barcode_sequence)
            if self.lane_number is not None:
                fq = "%s_L%03d" % (fq,self.lane_number)
            if self.read_number is not None:
                fq = "%s_R%d" % (fq,self.read_number)
        return fq

class AnalysisDir:
    """Class describing an analysis directory

    Conceptually an analysis directory maps onto a sequencing run.
    It consists of one or more sets of samples from that run,
    which are represented by subdirectories.

    It is also possible to have one or more subdirectories containing
    outputs from the CASAVA or bclToFastq processing software.

    """
    def __init__(self,analysis_dir):
        """Create a new AnalsysDir instance for a specified directory

        Arguments:
          analysis_dir: name (and path) to analysis directory

        """
        # Store location
        self._analysis_dir = os.path.abspath(analysis_dir)
        self._name = os.path.basename(analysis_dir)
        self._bcl2fastq_dirs = []
        self._project_dirs = []
        self.sequencing_data = []
        self.projects = []
        # Look for outputs from bclToFastq and analysis projects
        logging.debug("Examining subdirectories of %s" % self._analysis_dir)
        for dirn in list_dirs(self._analysis_dir):
            # Look for sequencing data
            try:
                data = IlluminaData.IlluminaData(self._analysis_dir,
                                                 unaligned_dir=dirn)
                logging.debug("- %s: sequencing data" % dirn)
                self._bcl2fastq_dirs.append(dirn)
                self.sequencing_data.append(data)
                continue
            except IlluminaData.IlluminaDataError, ex:
                pass
            # Look for analysis data
            data = AnalysisProject(self._name,
                                   os.path.join(self._analysis_dir,dirn))
            if data.is_analysis_dir:
                logging.debug("- %s: project directory" % dirn)
                self._project_dirs.append(dirn)
                self.projects.append(data)
                continue
            # Unidentified contents
            logging.debug("- %s: unknown" % dirn)

    @property
    def n_projects(self):
        """Return number of projects found

        """
        return len(self.projects)

    @property
    def n_sequencing_data(self):
        """Return number of sequencing data dirs found

        """
        return len(self.sequencing_data)
        
class AnalysisProject:
    """Class describing an analysis project

    Conceptually an analysis project consists of a set of samples
    from a single sequencing experiment, plus associated data e.g.
    QC results.

    Practically an analysis project is represented by a directory
    with a set of fastq files.

    Provides the following properties:

    name        : name of the project
    dirn        : associated directory (full path)
    samples     : list of AnalysisSample objects
    fastq_dir   : subdirectory holding fastq files

    
    multiple_fastqs: True if at least one sample has more than one fastq
                     file per read associated with it
    fastq_format: either 'fastqgz' or 'fastq'

    There is also an 'info' property with the following additional
    properties:

    run         : run name
    user        : user name
    PI          : PI name
    library_type: library type, either None or e.g. 'RNA-seq' etc
    organism    : organism, either None or e.g. 'Human' etc
    platform    : sequencing platform, either None or e.g. 'miseq' etc
    comments    : additional comments, either None or else string of text
    paired_end  : True if data is paired end, False if not

    """
    def __init__(self,name,dirn,user=None,PI=None,library_type=None,
                 organism=None,run=None,comments=None,platform=None):
        """Create a new AnalysisProject instance

        Arguments:
          name: name of the project
          dirn: project directory (can be full or relative path)
          user: optional, specify name of the user
          PI: optional, specify name of the principal investigator
          library_type: optional, specify library type e.g. 'RNA-seq',
            'miRNA' etc
          organism: optional, specify organism e.g. 'Human', 'Mouse'
            etc
          platform: optional, specify sequencing platform e.g 'miseq'
          run: optional, name of the run
          comments: optional, free text comments associated with the
            run

        """
        self.name = name
        self.dirn = os.path.abspath(dirn)
        self.fastq_dir = None
        self.fastq_format = None
        self.samples = []
        self.info = AnalysisProjectInfo()
        self.info_file = os.path.join(self.dirn,"README.info")
        self.populate()
        # (Re)set metadata
        if run is not None:
            self.info['run'] = run
        if user is not None:
            self.info['user'] = user
        if PI is not None:
            self.info['PI'] = PI
        if library_type is not None:
            self.info['library_type'] = library_type
        if organism is not None:
            self.info['organism'] = organism
        if platform is not None:
            self.info['platform'] = platform
        if comments is not None:
            self.info['comments'] = comments

    def populate(self):
        """Populate data structure from directory contents

        """
        if not os.path.exists(self.dirn):
            # Nothing to do, yet
            return
        # Locate fastq files
        self.fastq_dir = os.path.join(self.dirn,'fastqs')
        if not os.path.exists(self.fastq_dir):
            # If special 'fastqs' doesn't exist then
            # look in top level of project
            self.fastq_dir = self.dirn
        # Populate from fastq file names
        logging.debug("Acquiring fastqs...")
        fastqs = Pipeline.GetFastqGzFiles(self.fastq_dir)
        if not fastqs:
            logging.debug("No fastq.gz files found")
            fastqs = Pipeline.GetFastqFiles(self.fastq_dir)
            if not fastqs:
                logging.debug("No fastq files found")
            else:
                self.fastq_format = 'fastq'
        else:
            self.fastq_format = 'fastqgz'
        logging.debug("Assigning fastqs to samples...")
        for fastq in fastqs:
            # GetFastqGzFile returns a list of tuples
            for fq in fastq:
                name = AnalysisFastq(fq).sample_name
                try:
                    sample = self.get_sample(name)
                except KeyError:
                    sample = AnalysisSample(name)
                    self.samples.append(sample)
                sample.add_fastq(os.path.join(self.fastq_dir,fq))
        logging.debug("Listing samples and files:")
        for sample in self.samples:
            logging.debug("* %s: %s" % (sample.name,sample.fastq))
        # Get data from info file, if present
        if os.path.isfile(self.info_file):
            self.info.load(self.info_file)
        # Set paired_end flag for project
        paired_end = True
        for sample in self.samples:
            paired_end = (paired_end and sample.paired_end)
        self.info['paired_end'] = paired_end

    def create_directory(self,illumina_project=None,fastqs=None):
        """Create and populate analysis directory for an IlluminaProject

        Creates a new directory corresponding to the AnalysisProject
        object, and optionally also populates with links to FASTQ files
        from a supplied IlluminaProject object.

        The directory structure it creates is:

        dir/
           fastqs/
           logs/
           ScriptCode/

        It also creates an info file with metadata about the project.

        Arguments:
          illumina_project: (optional) populated IlluminaProject object
            from which the analysis directory will be populated
          fastqs: (optional) list of fastq files to import
    
        """
        logging.debug("Creating analysis directory for project '%s'" % self.name)
        # Check for & create directory
        if os.path.exists(self.dirn):
            logging.warning("Directory %s already exists" % self.dirn)
        else:
            logging.debug("Making analysis directory %s" % self.dirn)
            bcf_utils.mkdir(self.dirn,mode=0775)
        # Make a 'logs' directory
        log_dir = os.path.join(self.dirn,'logs')
        bcf_utils.mkdir(log_dir,mode=0775)
        # Make a 'ScriptCode' directory
        scriptcode_dir = os.path.join(self.dirn,"ScriptCode")
        bcf_utils.mkdir(scriptcode_dir,mode=0775)
        # Put a file in ScriptCode to make sure it's
        # not pruned on subsequent rsync operations
        fp = open(os.path.join(self.dirn,'ScriptCode','README.txt'),'w')
        fp.write("The ScriptCode directory is a place to put custom scripts and programs")
        fp.close()
        # Make a 'fastqs' directory
        fastqs_dir = os.path.join(self.dirn,"fastqs")
        bcf_utils.mkdir(fastqs_dir,mode=0775)
        # Check for & create links to fastq files
        if fastqs is None:
            # Make a list of fastqs to import from the supplied
            # IlluminaProject object
            fastqs = []
            if illumina_project is not None:
                for sample in illumina_project.samples:
                    for fastq in sample.fastq:
                        fastqs.append(os.path.join(sample.dirn,fastq))
        # Get mapping to unique names    
        fastq_names = IlluminaData.get_unique_fastq_names(fastqs)
        for fastq in fastqs:
            fastq_ln = os.path.join(fastqs_dir,fastq_names[fastq])
            if os.path.exists(fastq_ln):
                logging.warning("Link %s already exists" % fastq_ln)
            else:
                logging.debug("Linking to %s" % fastq)
                bcf_utils.mklink(fastq,fastq_ln,relative=True)
        # Populate
        self.populate()
        # Update metadata information summarising the samples
        n_samples = len(self.samples)
        if n_samples == 0:
            sample_description = "No samples"
        else:
            sample_description = "%s %s" % (n_samples,
                                            'sample' if n_samples == 1 else 'samples')
            sample_description += " (%s" % \
                                  bcf_utils.pretty_print_names(
                                      [s.name for s in self.samples])
            if self.multiple_fastqs:
                sample_description += ", multiple fastqs per sample"
            sample_description += ")"
        self.info['samples'] = sample_description
        # Save metadata
        self.info.save(self.info_file)

    @property
    def exists(self):
        """Check if analysis project directory already exists

        """
        return os.path.exists(self.dirn)

    @property
    def is_analysis_dir(self):
        """Determine if directory really is an analysis project

        """
        return len(self.samples) > 0

    @property
    def qc_dir(self):
        # Return path to qc dir, if present
        qc_dir = os.path.join(self.dirn,'qc')
        if os.path.exists(qc_dir):
            return qc_dir
        else:
            return None

    @property
    def qc(self):
        # Return IlluminaQCReporter object for this project
        if self.qc_dir is None:
            return None
        else:
            return qcreporter.IlluminaQCReporter(self.dirn,
                                                 data_format=self.fastq_format)

    @property
    def qc_report(self):
        # Create zipped QC report and return name of zip file
        qc_reporter = self.qc
        try:
            if self.verify_qc():
                return self.qc.zip()
        except AttributeError:
            logging.error("Failed to generate QC report")
        return None

    @property
    def multiple_fastqs(self):
        # Determine if there are multiple fastqs per sample
        if not len(self.samples):
            return False
        else:
            return reduce(lambda x,y: x and y,
                          [len(s.fastq_subset(read_number=1)) > 1 for s in self.samples])

    def verify_qc(self):
        # Verify if the QC was successful
        try:
            return self.qc.verify()
        except AttributeError:
            return False

    def get_sample(self,name):
        """Return sample that matches 'name'

        Arguments:
          name: name of a sample

        Returns:
          AnalysisSample object with the matching name; raises
          KeyError exception if no match is found.

        """
        for sample in self.samples:
            if sample.name == name: return sample
        raise KeyError, "No matching sample for '%s'" % name

    def get_samples(self,pattern):
        """Return list of sample matching pattern

        Arguments:
          pattern: simple 'glob' style pattern

        Returns:
          Python list of samples with names matching the supplied
          pattern (or an empty list if no names match).

        """
        samples = []
        for sample in self.samples:
            if bcf_utils.name_matches(sample.name,pattern):
                samples.append(sample)
        return samples

    def prettyPrintSamples(self):
        """Return a nicely formatted string describing the sample names

        Wraps a call to 'pretty_print_names' function.

        """
        return bcf_utils.pretty_print_names(self.samples)

class AnalysisSample:
    """Class describing an analysis sample

    An analysis sample consists of a set of fastqs file corresponding
    to single sample.

    AnalysisSample has the following properties:

    name      : name of the sample
    fastq     : list of fastq files associated with the sample
    paired_end: True if sample is paired end, False if not

    """

    def __init__(self,name):
        """Create a new AnalysisSample instance

        Arguments:
          name: sample name

        """
        self.name = name
        self.fastq = []
        self.paired_end = False

    def add_fastq(self,fastq):
        """Add a reference to a fastq file in the sample

        Arguments:
          fastq: name of the fastq file

        """
        self.fastq.append(fastq)
        # Sort fastq's into order
        self.fastq.sort()
        # Check paired-end status
        if not self.paired_end:
            fq = AnalysisFastq(fastq)
            if fq.read_number == 2:
                self.paired_end = True

    def fastq_subset(self,read_number=None,full_path=False):
        """Return a subset of fastq files from the sample

        Arguments:
          read_number: select subset based on read_number (1 or 2)
          full_path  : if True then fastq files will be returned
            with the full path, if False (default) then as file
            names only.

        Returns:
          List of fastq files matching the selection criteria.

        """
        # Build list of fastqs that match the selection criteria
        fastqs = []
        for fastq in self.fastq:
            fq = AnalysisFastq(fastq)
            if fq.read_number is None:
                logging.debug("Unable to determine read number for %s, assume R1" % fastq)
                fq_read_number = 1
            else:
                fq_read_number = fq.read_number
            if fq_read_number == read_number:
                if full_path:
                    fastqs.append(os.path.join(self.dirn,fastq))
                else:
                    fastqs.append(fastq)
        # Sort into dictionary order and return
        fastqs.sort()
        return fastqs

    def qc_sample(self,qc_dir,fastq):
        """Fetch IlluminaQCSample object for a fastq file

        Arguments:
          qc_dir: name of the QC directory
          fastq : fastq file to get the QC information for

        Returns:
          Populated IlluminaQCSample object.

        """
        name = bcf_utils.rootname(os.path.basename(fastq))
        return qcreporter.IlluminaQCSample(name,qc_dir)

    def verify_qc(self,qc_dir,fastq):
        """Check if QC completed for a fastq file

        Arguments:
          qc_dir: name of the QC directory
          fastq : fastq file to get the QC information for

        Returns:
          True if QC completed correctly, False otherwise.

        """
        return self.qc_sample(qc_dir,fastq).verify()

    def __repr__(self):
        """Implement __repr__ built-in

        Return string representation for the sample -
        i.e. the sample name.

        """
        return str(self.name)

class MetadataDict(bcf_utils.AttributeDictionary):
    """Class for storing metadata in an analysis project

    Provides storage for arbitrary data items in the form of
    key-value pairs, which can be saved to and loaded from
    an external file.

    The data items are defined on instantiation via a dictionary
    supplied to the 'attributes' argument. For example:

    Create a new metadata object:
    >>> metadata = MetadataDict(attributes={'salutation':'Salutation',
    ...                                     'valediction': 'Valediction'})
 
    The dictionary keys correspond to the keys in the MetadataDict
    object; the corresponding values are the keys that are used
    when saving and loading the data to and from a file.

    Set attributes:
    >>> metadata['salutation'] = 'hello'
    >>> metadata['valediction'] = 'goodbye'

    Retrieve values:
    >>> print "Salutation is %s" % metadata.salutation

    Save to file:
    >>> metadata.save('metadata.tsv')

    Load data from a file:
    >>> metadata = MetadataDict('metadata.tsv')
    or
    >>> metadata = MetadataDict()
    >>> metadata.load('metadata.tsv')

    The external file storage is intended to be readable by
    humans so longer names are used to describe the keys; also
    Python None values are stored as '.', and True and False
    values are stored as 'Y' and 'N' respectively. These values
    are automatically converted back to the Python equivalents
    on reload.

    """

    def __init__(self,attributes=dict(),order=None,filen=None):
        """Create a new MetadataDict object

        By default an empty metadata object is created
        i.e. all attributes will have be None.

        If an input file is specified then the attributes
        will be assigned values according to the key-value
        pairs in that file.

        Arguments:
          attributes: dictionary defining metadata items
          filen: (optional) name of the tab-delimited file
            with key-value pairs to load in.

        """
        bcf_utils.AttributeDictionary.__init__(self)
        self.__filen = filen
        # Set up empty metadata attributes
        self.__attributes = attributes
        for key in self.__attributes:
            self[key] = None
        if self.__filen:
            # Load data from external file
            load(self,self.__filen)
        # Set up order of keys for output
        if order is None:
            self.__key_order = self.__attributes.keys()
            self.__key_order.sort()
        else:
            # Use supplied key order
            self.__key_order = []
            for key in order:
                if key in self.__attributes:
                    self.__key_order.append(key)
                else:
                    raise KeyError,"Key '%s' not defined in attributes"
            # Append keys not explicitly listed in the order
            extra_keys = []
            for key in self.__attributes:
                if key not in self.__key_order:
                    extra_keys.append(key)
            if extra_keys:
                extra_keys.sort()
                self.__key_order.extend(extra_keys)

    def __setitem__(self,key,value):
        if key in self.__attributes:
            bcf_utils.AttributeDictionary.__setitem__(self,key,value)
        else:
            raise AttributeError,"Key '%s' not defined" % key

    def load(self,filen):
        """Load key-value pairs from a tab-delimited file
        
        Loads the key-value pairs from a previously created
        tab-delimited file written by the 'save' method.

        Note that this overwrites any existing values
        already assigned to keys within the metadata object.

        Arguments:
          filen: name of the tab-delimited file with key-value
            pairs

        """
        self.__filen = filen
        metadata = TabFile.TabFile(filen)
        for line in metadata:
            try:
                # Get data from file and convert special values
                # to Python equivalents
                attr,value = line[0],line[1]
                if value == '.' or value == 'None':
                    value = None
                elif value == 'Y' or value == 'True':
                    value = True
                elif value == 'N' or value == 'False':
                    value = False
                # Locate dictionary key matching file key
                found_key = False
                for key in self.__attributes:
                    if self.__attributes[key] == attr:
                        self[key] = value
                        found_key = True
                        break
                if not found_key:
                    logging.debug("Unrecognised key in %s: %s" % (filen,key))
            except IndexError:
                logging.warning("Bad line in %s: %s" % (filen,line))

    def save(self,filen=None):
        """Save metadata to tab-delimited file

        Writes key-value paires to a tab-delimited file.
        The data can be recovered using the 'load' method.
 
        Note that if the specified file already exists then
        it will be overwritten.

        Arguments:
          filen: name of the tab-delimited file with key-value
            pairs; if None then the file specified when the
            object was instantiated will be used instead.

        """
        metadata = TabFile.TabFile()
        for key in self.__key_order:
            # Retrieve value and convert to appropriate
            # format for persistent storage
            value = self[key]
            if value is None:
                value = '.'
            elif value is True:
                value = 'Y'
            elif value is False:
                value = 'N'
            # Get the equivalent file key
            attr = self.__attributes[key]
            # Store in the file
            metadata.append(data=(attr,value))
        # Write the file
        if filen is not None:
            self.__filen = filen
        metadata.write(self.__filen)

class AnalysisDirMetadata(MetadataDict):
    """Class for storing metadata in an analysis project

    Provides a set of metadata items which are loaded from
    and saved to an external file.

    The data items are:

    analysis_dir: path to the analysis directory
    data_dir: path to the directory holding the raw sequencing data
    platform: sequencing platform e.g. 'miseq'
    sample_sheet: path to the customised SampleSheet.csv file
    bases_mask: bases mask string
    project_metadata: name of the project metadata file
    primary_data_dir: directory used to hold copies of primary data
    unaligned_dir: output directory for bcl2fastq conversion
    stats_file: name of file with statistics about the run

    """
    def __init__(self,filen=None):
        """Create a new AnalysisDirMetadata object

        Arguments:
          filen: (optional) name of the tab-delimited file
            with key-value pairs to load in.

        """
        MetadataDict.__init__(self,
                              attributes = {
                                  'analysis_dir':'analysis_dir',
                                  'data_dir':'data_dir',
                                  'platform':'platform',
                                  'sample_sheet':'sample_sheet',
                                  'bases_mask':'bases_mask',
                                  'project_metadata':'project_metadata',
                                  'primary_data_dir':'primary_data_dir',
                                  'unaligned_dir':'unaligned_dir',
                                  'stats_file':'stats_file',
                                  'source': 'source',
                                  'run_number': 'run_number',
                                  'assay': 'assay'
                              },
                              filen=filen)

class AnalysisProjectInfo(MetadataDict):
    """Class for storing metadata in an analysis project

    Provides a set of metadata items which are loaded from
    and saved to an external file.

    The data items are:

    run: the name of the sequencing run
    user: the user associated with the project
    PI: the principal investigator associated with the project
    organism: the organism associated with the project
    library_type: the library type e.g. 'RNA-seq'
    platform: the platform name e.g. 'miseq'
    paired_end: True if the data is paired end, False if not
    samples: textual description of the samples in the project
    comments: free-text comments

    """
    def __init__(self,filen=None):
        """Create a new AnalysisProjectInfo object

        Arguments:
          filen: (optional) name of the tab-delimited file
            with key-value pairs to load in.

        """
        MetadataDict.__init__(self,
                              attributes = {
                                  'run':'Run',
                                  'platform':'Platform',
                                  'user':'User',
                                  'PI':'PI',
                                  'organism':'Organism',
                                  'library_type':'Library type',
                                  'paired_end':'Paired_end',
                                  'samples':'Samples',
                                  'comments':'Comments',
                              },
                              order = (
                                  'run',
                                  'platform',
                                  'user',
                                  'PI',
                                  'organism',
                                  'library_type',
                                  'paired_end',
                                  'samples',
                                  'comments',
                              ),
                              filen=filen)

class ProjectMetadataFile(TabFile.TabFile):
    """File containing metadata about multiple projects in analysis dir
 
    The file consists of a header line plus one line per project
    with the following tab-delimited fields:

    Project: name of the project
    Samples: list/description of sample names
    User: name(s) of the associated user(s)
    Library: the library type
    Organism: name(s) of the organism(s)
    PI: name(s) of the associated principal investigator(s)
    Comments: free text containing additional information
              about the project

    Any fields set to None will be written to file with a '.'
    placeholder.

    """
    def __init__(self,filen=None):
        """Create a new ProjectsMetadataFile instance

        Arguments:
          filen: (optional) name of an existing file to read
            projects in from.

        """
        self.__filen = filen
        TabFile.TabFile.__init__(self,filen=filen,
                                 column_names=('Project',
                                               'Samples',
                                               'User',
                                               'Library',
                                               'Organism',
                                               'PI',
                                               'Comments'),
                                 first_line_is_header=True)

    def add_project(self,project_name,sample_names,user=None,
                    library_type=None,organism=None,PI=None,
                    comments=None):
        """Add information about a project into the file

        Arguments:
          project_name: name of the project
          sample_names: Python list of sample names
          user: (optional) user name(s)
          library_type: (optional) library type
          organism: (optional) organism(s)
          PI: (optional) principal investigator name(s)
          comments: (optional) additional information about
            the project

        """
        # Add project info to the metadata file
        self.append(data=(project_name,
                          ','.join(sample_names),
                          '.' if user is None else user,
                          '.' if library_type is None else library_type,
                          '.' if organism is None else organism,
                          '.' if PI is None else PI,
                          '.' if comments is None else comments))

    def project(self,name):
        """Return AttributeDictionary for a project

        """
        raise NotImplementedError

    def save(self,filen=None):
        """Save the data back to file

        Arguments:
          filen: name of the file to save to (if not specified then
            defaults to the same file as data was read in from)

        """
        if filen is not None:
            self.__filen = filen
        self.write(filen=self.__filen,include_header=True)

#######################################################################
# Functions
#######################################################################

def fetch_runner(definition):
    """Return job runner instance based on a definition string

    Given a definition string, returns an appropriate runner
    instance.

    Definitions are of the form:

      RunnerName[(args)]

    RunnerName can be 'SimpleJobRunner' or 'GEJobRunner'.
    If '(args)' are also supplied then these are passed to
    the job runner on instantiation (only works for
    GE runners).

    """
    if definition.startswith('SimpleJobRunner'):
        return JobRunner.SimpleJobRunner(join_logs=True)
    elif definition.startswith('GEJobRunner'):
        if definition.startswith('GEJobRunner(') and definition.endswith(')'):
            ge_extra_args = definition[len('GEJobRunner('):len(definition)-1].split(' ')
            return JobRunner.GEJobRunner(ge_extra_args=ge_extra_args)
        else:
            return JobRunner.GEJobRunner()
    raise Exception,"Unrecognised runner definition: %s" % definition

def bases_mask_is_paired_end(bases_mask):
    # Determine if run is paired end based on bases mask string
    non_index_reads = []
    for read in bases_mask.split(','):
        try:
            read.index('I')
        except ValueError:
            non_index_reads.append(read)
    if len(non_index_reads) == 2:
        # Paired end
        return True
    elif len(non_index_reads) < 2:
        # Single end
        return False
    else:
        # An error?
        raise Exception, "Bad bases mask '%s'?" % bases_mask

def split_user_host_dir(location):
    # Split a location of the form [[user@]host:]dir into its
    # user, hostname and directory components
    try:
        location = location.strip()
    except AttributeError:
        # Not a string?
        logging.error("Bad input to split_user_host_dir: '%s'" % location)
        return (None,None,None)
    if not location:
        return (None,None,None)
    try:
        location.index(':')
        location,dirn = location.split(':')
        try:
            location.index('@')
            user,host = location.split('@')
        except ValueError:
            user = None
            host = location
    except ValueError:
        user = None
        host = None
        dirn = location
    return (user,host,dirn)

def list_dirs(parent,matches=None,startswith=None):
    # Return list of subdirectories relative to 'parent'
    # If startswith is set then return subset that starts with
    # the specified string
    # If matches is set then return an exact match only
    dirs = []
    for d in os.listdir(parent):
        if os.path.isdir(os.path.join(parent,d)):
            if startswith is None or d.startswith(startswith):
                if matches is None or d == matches:
                    dirs.append(d)
    dirs.sort()
    return dirs
