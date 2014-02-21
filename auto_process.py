#!/bin/env python
#
#     auto_process.py: automated processing of Illumina sequence data
#     Copyright (C) University of Manchester 2013 Peter Briggs
#
#########################################################################
#
# auto_process.py
#
#########################################################################

"""
First attempt at an automated data processing & QC pipeline in Python

Implements a program for automating stages of a standard protocol for
processing and QC'ing Illumina sequencing data.

The stages are:

    setup
    make_fastqs
    setup_analysis_dirs
    run_qc
    archive
    publish_qc

The 'setup' stage creates an analysis directory and acquires the basic
data about the sequencing run from a source directory. Subsequent stages
should be run in sequence to create fastq files, set up analysis
directories for each project, and run QC scripts for each sample in
each project.

"""

__version__ = "0.0.32"

#######################################################################
# Imports
#######################################################################

import sys
import os
import re
import subprocess
import logging
import optparse
import shutil
import time
import IlluminaData
import platforms
import TabFile
import FASTQFile
import JobRunner
import Pipeline
import bcf_utils
import analyse_illumina_run
import qcreporter
import bclToFastq
import applications
import auto_process_utils
import auto_process_settings

#######################################################################
# Classes
#######################################################################

class AutoProcess:
    # Class implementing an automatic fastq generation and QC
    # processing procedure

    def __init__(self,analysis_dir=None):
        # analysis_dir: name/path for existing analysis directory
        self.params = auto_process_utils.AnalysisDirMetadata()
        self.analysis_dir = analysis_dir
        if self.analysis_dir is not None:
            self.analysis_dir = os.path.abspath(analysis_dir)
            self.load_parameters()
            self.params['analysis_dir'] = self.analysis_dir

    def add_directory(self,sub_dir):
        # Add a directory to the AutoProcess object
        dirn = os.path.join(self.analysis_dir,sub_dir)
        self.create_directory(dirn)
        return dirn

    def create_directory(self,dirn):
        # Make the specified directory, and any leading directories
        # that don't already exist
        if not os.path.exists(dirn):
            dir_path = os.sep
            for sub_dir in dirn.split(os.sep):
                dir_path = os.path.join(dir_path,sub_dir)
                if not os.path.exists(dir_path):
                    print "Making %s" % dir_path
                    bcf_utils.mkdir(dir_path)

    def load_parameters(self):
        # Get parameters from info file
        #
        # check for info file
        info_file_name = os.path.join(self.analysis_dir,'auto_process.info')
        if not os.path.isfile(info_file_name):
            raise Exception, "No info file %s" % info_file_name
        # Read contents of info file and assign values
        logging.debug("Loading settings from %s" % info_file_name)
        self.params.load(info_file_name)

    def save_parameters(self):
        # Save parameters to info file
        #
        info_file_name = os.path.join(self.analysis_dir,'auto_process.info')
        self.params.save(info_file_name)

    def log_path(self,*args):
        # Return path appended to log directory
        # Use for getting paths of files under the logs directory
        return os.path.join(self.log_dir,*args)

    def __del__(self):
        logging.debug("Saving parameters to file")
        self.save_parameters()

    @property
    def log_dir(self):
        # Generate and return full path to log directory
        return self.add_directory('logs')

    @property
    def tmp_dir(self):
        # Generate and return full path to tmp directory
        return self.add_directory('tmp')

    def setup(self,data_dir,analysis_dir=None):
        # Set up the initial analysis directory
        #
        # This does all the initialisation of the analysis directory
        # and processing parameters
        #
        # Arguments:
        # data_dir: source data directory
        # analysis_dir: corresponding analysis dir
        data_dir = data_dir.rstrip(os.sep)
        if analysis_dir is None:
            self.analysis_dir = os.path.join(os.getcwd(),
                                             os.path.basename(data_dir))+'_analysis'
        else:
            self.analysis_dir = os.path.abspath(analysis_dir)
        # Create the analysis directory structure
        if not os.path.exists(self.analysis_dir):
            # Create directory structure
            self.create_directory(self.analysis_dir)
        else:
            # Directory already exists
            logging.warning("Analysis directory already exists")
            # check for info file
            info_file_name = os.path.join(self.analysis_dir,'auto_process.info')
            if os.path.isfile(info_file_name):
                self.load_parameters()
            else:
                logging.warning("No info file found in %s" % self.analysis_dir)
        # Identify missing data and attempt to acquire
        # Sequencing platform
        platform = self.params.platform
        if platform is None:
            print "Identifying platform from data directory name"
            platform = platforms.get_sequencer_platform(data_dir)
        print "Platform identified as '%s'" % platform
        # Custom SampleSheet.csv file
        custom_sample_sheet = self.params.sample_sheet
        if custom_sample_sheet is None:
            print "Acquiring sample sheet..."
            tmp_sample_sheet = os.path.join(self.tmp_dir,'SampleSheet.csv')
            rsync = applications.general.rsync(os.path.join(data_dir,
                                                            'Data/Intensities/BaseCalls/SampleSheet.csv'),
                                               self.tmp_dir)
            print "%s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.sample_sheet.log'))
            custom_sample_sheet = os.path.join(self.analysis_dir,'custom_SampleSheet.csv')
            sample_sheet = make_custom_sample_sheet(tmp_sample_sheet,
                                                    custom_sample_sheet)
            print "Keeping copy of original sample sheet"
            os.rename(tmp_sample_sheet,os.path.join(self.analysis_dir,'SampleSheet.orig.csv'))
        sample_sheet = IlluminaData.CasavaSampleSheet(custom_sample_sheet)
        print "Sample sheet '%s'" % custom_sample_sheet
        # Bases mask
        bases_mask = self.params.bases_mask
        if bases_mask is None:
            print "Acquiring RunInfo.xml to determine bases mask..."
            tmp_run_info = os.path.join(self.tmp_dir,'RunInfo.xml')
            rsync = applications.general.rsync(os.path.join(data_dir,'RunInfo.xml'),
                                               self.tmp_dir)
            status = rsync.run_subprocess(log=self.log_path('rsync.run_info.log'))
            bases_mask = get_bases_mask(tmp_run_info,custom_sample_sheet)
            os.remove(tmp_run_info)
        print "Corrected bases mask: %s" % bases_mask
        # Print the predicted ouputs
        projects = sample_sheet.predict_output()
        print "Predicted output from sample sheet:"
        print "Project\tSample\tFastq"
        for project in projects:
            project_name = project[8:]
            sample_names = []
            for sample in projects[project]:
                sample_name = sample[7:]
                for fastq_base in projects[project][sample]:
                    print "%s\t%s\t%s" % (project_name,sample_name,fastq_base)
                sample_names.append(sample_name)
        # Store the parameters
        self.params['data_dir'] = data_dir
        self.params['analysis_dir'] = self.analysis_dir
        self.params['platform'] = platform
        self.params['sample_sheet'] = custom_sample_sheet
        self.params['bases_mask'] = bases_mask

    def setup_from_fastq_dir(self,analysis_dir,fastq_dir):
        # Do setup for an existing directory containing fastq files
        # with the same structure as that produced by CASAVA and bcl2fastq
        #
        # Assumes that the files are in a subdirectory of the analysis
        # directory specified by the 'fastq_dir' argument, and
        # that within that they are arranged in the structure
        # 'Project_<name>/Sample_<name>/<fastq>'
        self.analysis_dir = os.path.abspath(analysis_dir)
        # Create directory structure
        self.create_directory(self.analysis_dir)
        # Get information
        print "Identifying platform from data directory name"
        platform = platforms.get_sequencer_platform(analysis_dir)
        # Store the parameters
        self.params['analysis_dir'] = self.analysis_dir
        self.params['unaligned_dir'] = fastq_dir
        self.params['platform'] = platform
        # Generate statistics
        self.generate_stats()
        # Make a 'projects.info' metadata file
        self.make_project_metadata_file()

    def make_project_metadata_file(self,project_metadata_file='projects.info'):
        # Generate a project metadata file based on the fastq
        # files and directory structure
        filen = os.path.join(self.params.analysis_dir,project_metadata_file)
        illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        print "Project metadata file: %s" % filen
        if not os.path.exists(filen):
            # Create a new metadata file
            project_metadata = ProjectMetadataFile()
            print "Project\tSample\tFastq"
            for project in illumina_data.projects:
                project_name = project.name
                sample_names = []
                for sample in project.samples:
                    sample_name = sample.name
                    for fastq in sample.fastq:
                        print "%s\t%s\t%s" % (project_name,sample_name,fastq)
                    sample_names.append(sample_name)
                project_metadata.add_project(project_name,sample_names)
        else:
            # Load existing file and check for consistency
            print "Metadata file already exists, checking"
            project_metadata = ProjectMetadataFile(filen)
            # Check that each project listed actually exists
            bad_projects = []
            for line in project_metadata:
                pname = line['Project']
                try:
                    illumina_data.get_project(pname)
                except IlluminaDataError:
                    # Project doesn't exist
                    logging.warning("Project '%s' listed in metadata file doesn't exist" \
                                    % pname)
                    bad_projects.append(line)
            # Remove bad projects
            for bad_project in bad_projects:
                del(bad_project)
            # Check that all actual projects are listed
            for project in illumina_data.projects:
                if len(project_metadata.lookup('Project',project.name)) == 0:
                    # Project not listed, add it
                    logging.warning("Project '%s' not listed in metadata file" % project.name)
                    sample_names = []
                    for sample in project.samples:
                        sample_name = sample.name
                        for fastq in sample.fastq:
                            print "%s\t%s\t%s" % (project_name,sample_name,fastq)
                        sample_names.append(sample_name)
                    project_metadata.add_project(project_name,sample_names)
        # Finish
        project_metadata.save(filen)
        self.params['project_metadata'] = project_metadata_file
        print "Project metadata in %s" % self.params.project_metadata

    def get_analysis_projects(self,pattern=None):
        # Return the analysis projects in a list
        #
        # By default returns all projects within the analysis
        #
        # If the 'pattern' is not None then it should be a simple pattern
        # used to match against available names to select a subset of
        # projects (see bcf_utils.name_matches).
        project_metadata = ProjectMetadataFile(os.path.join(self.analysis_dir,
                                                            self.params.project_metadata))
        projects = []
        for line in project_metadata:
            name = line['Project']
            if pattern is not None:
                if not bcf_utils.name_matches(name,pattern):
                    # Name failed to match, ignore
                    continue
            print "Acquiring data for project %s" % name
            project_dir = os.path.join(self.analysis_dir,name)
            if os.path.isdir(project_dir):
                projects.append(
                    auto_process_utils.AnalysisProject(
                        name,
                        os.path.join(self.analysis_dir,name)))
            else:
                logging.warning("No directory found for project '%s'" % name)
        return projects

    def get_primary_data(self):
        # Copy the primary sequencing data (bcl files etc) to a local area
        # using rsync
        data_dir = self.params.data_dir
        self.params["primary_data_dir"] = self.add_directory('primary_data')
        try:
            rsync = applications.general.rsync(data_dir,self.params.primary_data_dir,
                                               prune_empty_dirs=True,
                                               extra_options=('--include=*/',
                                                              '--include=Data/**',
                                                              '--include=RunInfo.xml',
                                                              '--include=SampleSheet.csv',
                                                              '--exclude=*'))
            print "Running %s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.primary_data.log'))
        except Exception, ex:
            logging.error("Exception getting primary data: %s" % ex)
            status = -1
        if status != 0:
            logging.error("Failed to acquire primary data (status %s)" % status)
        return status
        
    def bcl_to_fastq(self,ignore_missing_bcl=False,ignore_missing_stats=False,
                     skip_rsync=False,keep_primary_data=False,generate_stats=False):
        # Convert bcl files to fastq
        #
        # Arguments:
        # ignore_missing_bcl  : if True then run bcl2fastq with --ignore-missing-bcl
        # ignore_missing_stats: if True then run bcl2fastq with --ignore-missing-stats
        # skip_rsync          : if True then don't rsync primary data at the start of
        #                       bcl2fastq conversion
        # keep_primary_data   : if True then don't remove primary data at
        #                       the end of bcl2fastq conversion
        # generate_stats      : if True then (re)generate statistics file for fastqs
        #
        # Directories
        analysis_dir = self.params.analysis_dir
        self.params['unaligned_dir'] = 'bcl2fastq'
        bcl2fastq_dir = self.add_directory(self.params.unaligned_dir)
        sample_sheet = self.params.sample_sheet
        if self.verify_bcl_to_fastq():
            print "Bcl to fastq outputs already present"
            # Check for project metadata file
            self.make_project_metadata_file()
            # (Re)generate stats?
            if generate_stats:
                self.generate_stats()
            return
        # Fetch primary data
        if not skip_rsync:
            if self.get_primary_data() != 0:
                logging.error("Failed to acquire primary data")
                raise Exception, "Failed to acquire primary data"
        primary_data = os.path.join(self.params.primary_data_dir,
                                    os.path.basename(self.params.data_dir))
        # Get info about the run
        print "Primary data dir    : %s" % primary_data
        illumina_run = IlluminaData.IlluminaRun(primary_data)
        bases_mask = self.params.bases_mask
        nmismatches = bclToFastq.get_nmismatches(bases_mask)
        print "%s" % illumina_run.run_dir
        print "Platform            : %s" % illumina_run.platform
        print "Bcl format          : %s" % illumina_run.bcl_extension
        print "Sample sheet        : %s" % sample_sheet
        print "Bases mask          : %s" % bases_mask
        print "Nmismatches         : %d (determined from bases mask)" % nmismatches
        print "Ignore missing bcl  : %s" % ignore_missing_bcl
        print "Ignore missing stats: %s" % ignore_missing_stats
        # Set up runner
        runner = auto_process_settings.runners.bcl2fastq
        runner.log_dir(self.log_dir)
        # Run bcl2fastq
        bcl2fastq = applications.Command('bclToFastq.py',
                                         '--use-bases-mask',bases_mask,
                                         '--nmismatches',nmismatches,
                                         '--ignore-missing-control')
        if ignore_missing_bcl:
            bcl2fastq.add_args('--ignore-missing-bcl')
        if ignore_missing_stats:
            bcl2fastq.add_args('--ignore-missing-stats')
        bcl2fastq.add_args(primary_data,
                           bcl2fastq_dir,
                           sample_sheet)
        print "Running %s" % bcl2fastq
        bcl2fastq_job = Pipeline.Job(runner,
                                     'bclToFastq',
                                     os.getcwd(),
                                     bcl2fastq.command,
                                     bcl2fastq.args)
        bcl2fastq_job.start()
        bcl2fastq_job.wait()
        print "bcl2fastq completed"
        # Verify outputs
        illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        if not analyse_illumina_run.verify_run_against_sample_sheet(illumina_data,
                                                                    sample_sheet):
            logging.error("Failed to verify bcl to fastq outputs against sample sheet")
            raise Exception, "Failed to verify bcl to fastq outputs against sample sheet"
        # Remove primary data
        if not keep_primary_data:
            self.remove_primary_data()
        # Generate statistics
        self.generate_stats()
        # Make a 'projects.info' metadata file
        self.make_project_metadata_file()

    def generate_stats(self,stats_file='statistics.info'):
        # Generate statistics for initial fastq files from bcl2fastq
        # Set up runner
        runner = auto_process_settings.runners.stats
        runner.log_dir(self.log_dir)
        # Generate statistics
        fastq_statistics = applications.Command('fastq_statistics.py',
                                                '--unaligned',self.params.unaligned_dir,
                                                '--output',
                                                os.path.join(self.params.analysis_dir,
                                                             stats_file),
                                                self.params.analysis_dir)
        print "Generating statistics: running %s" % fastq_statistics
        fastq_statistics_job = Pipeline.Job(runner,
                                            'fastq_statistics',
                                            self.params.analysis_dir,
                                            fastq_statistics.command,
                                            fastq_statistics.args)
        fastq_statistics_job.start()
        fastq_statistics_job.wait()
        self.params['stats_file'] = stats_file
        print "Statistics generation completed: %s" % self.params.stats_file

    def remove_primary_data(self):
        # Remove primary data
        primary_data = os.path.join(self.params.primary_data_dir,
                                    os.path.basename(self.params.data_dir))
        if os.path.isdir(primary_data):
            print "Removing copy of primary data in %s" % primary_data
            shutil.rmtree(primary_data)

    def verify_bcl_to_fastq(self):
        # Check that bcl to fastq outputs match sample sheet predictions
        bcl_to_fastq_dir = os.path.join(self.analysis_dir,self.params.unaligned_dir)
        if not os.path.isdir(bcl_to_fastq_dir):
            # Directory doesn't exist
            return False
        # Try to create an IlluminaData object
        try:
            illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                      unaligned_dir=self.params.unaligned_dir)
        except IlluminaData.IlluminaDataError, ex:
            # Failed to initialise
            logging.debug("Failed to get information from %s: %s" % (bcl_to_fastq_dir,ex))
            return False
        # Do check
        return analyse_illumina_run.verify_run_against_sample_sheet(illumina_data,
                                                                    self.params.sample_sheet)

    def setup_analysis_dirs(self):
        # Construct and populate the analysis directories for each project
        illumina_data = IlluminaData.IlluminaData(self.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        project_metadata = TabFile.TabFile(os.path.join(analysis_dir,
                                                        self.params.project_metadata),
                                           first_line_is_header=True)
        # Sanity check that the project data file has been populated
        got_project_data = True
        for line in project_metadata:
            for item in ('User','PI','Organism','Library',):
                if line[item] == '.':
                    logging.warning("Missing data from %s for '%s': %s" %
                                    (self.params.project_metadata,
                                     line['Project'],item))
                    got_project_data = False
        if not got_project_data:
            logging.error("Missing project metadata")
            raise Exception, "Missing project metadata"
        # Create the projects
        n_projects = 0
        for line in project_metadata:
            # Acquire the run name
            if self.params.data_dir is not None:
                run_name = os.path.basename(self.params.data_dir)
            else:
                run_name = os.path.basename(self.params.analysis_dir)
            # Look up project data
            project_name = line['Project']
            user = line['User']
            PI = line['PI']
            organism = line['Organism']
            library_type = line['Library']
            comments = line['Comments']
            # Create the project
            project = auto_process_utils.AnalysisProject(project_name,
                                                         os.path.join(self.analysis_dir,
                                                                      project_name),
                                                         user=user,
                                                         PI=PI,
                                                         organism=organism,
                                                         library_type=library_type,
                                                         run=run_name,
                                                         comments=comments,
                                                         platform=self.params.platform)
            print "Creating project: '%s'" % project_name
            project.create_directory(illumina_data.get_project(project_name))
            n_projects += 1
        # Tell us how many were made
        print "Created %d projects" % n_projects

    def run_qc(self,projects=None,max_jobs=4):
        # Run QC pipeline for all projects
        #
        # Tests whether QC outputs already exist and only runs
        # QC for those files where the outputs are not all present
        #
        # projects: specify a pattern to match one or more projects to
        # run the QC for (default is to run QC for all projects)
        #
        # Process pattern matching
        if projects is None:
            project_pattern = '*'
            sample_pattern = '*'
        else:
            project_pattern = projects.split('/')[0]
            try:
                sample_pattern = projects.split('/')[1]
            except IndexError:
                sample_pattern = '*'
        # Setup a pipeline runner
        qc_runner = auto_process_settings.runners.qc
        pipeline = Pipeline.PipelineRunner(qc_runner,
                                           max_concurrent_jobs=max_jobs)
        # Get project dir data
        projects = self.get_analysis_projects(project_pattern)
        # Check we have projects
        if len(projects) == 0:
            logging.warning("No projects found for QC analysis")
            return
        # Look for samples with no/invalid QC outputs and populate
        # pipeline with the associated fastq.gz files
        for project in projects:
            print "*** Setting up QC for %s ***" % project.name
            # Make the qc directory if it doesn't exist
            qc_dir = os.path.join(project.dirn,'qc')
            if not os.path.exists(qc_dir):
                print "Making 'qc' subdirectory"
                bcf_utils.mkdir(qc_dir,mode=0775)
            # Loop over samples and queue up those where the QC
            # isn't validated
            samples = project.get_samples(sample_pattern)
            if len(samples) == 0:
                logging.warning("No samples found for QC analysis in project '%s'" %
                                project.name)
            for sample in samples:
                print "Examining files in sample %s" % sample.name
                for fq in sample.fastq:
                    if sample.verify_qc(qc_dir,fq):
                        logging.debug("\t%s: QC verified" % fq)
                    else:
                        print "\t%s: setting up QC run" % os.path.basename(fq)
                        group = "%s.%s" % (project.name,sample.name)
                        fastq = os.path.join(project.dirn,'fastqs',fq)
                        label = str(auto_process_utils.AnalysisFastq(fq))
                        pipeline.queueJob(project.dirn,'illumina_qc.sh',(fastq,),
                                          label=label,group=group)
        # Run the pipeline
        if pipeline.nWaiting() > 0:
            pipeline.run()
        # Verify the outputs
        for project in projects:
            if not project.verify_qc():
                logging.error("QC failed for one or more samples in %s" % project.name)
            else:
                print "QC okay, generating report for %s" % project.name
                project.qc_report

    def copy_to_archive(self,platform=None,year=None):
        # Copy the analysis directory and contents to an archive area
        archive_dir = auto_process_settings.archive.dirn
        if archive_dir is None:
            raise Exception, "No archive directory specified in settings"
        # Construct subdirectory structure i.e. platform and year
        if platform is None:
            platform = self.params.platform
        if platform is None:
            raise Exception, "No platform specified (use --platform option?)"
        if year is None:
            year = time.strftime("%Y")
        archive_dir = os.path.join(archive_dir,year,platform)
        print "Copying to archive directory: %s" % archive_dir
        try:
            rsync = applications.general.rsync(self.analysis_dir,archive_dir,
                                               prune_empty_dirs=True)
            print "Running %s" % rsync
            status = rsync.run_subprocess(log=self.log_path('rsync.archive.log'))
        except Exception, ex:
            logging.error("Exception rsyncing to archive: %s" % ex)
            status = -1
        if status != 0:
            logging.error("Failed to rsync to archive (returned status %d)" % status)

    def log_analysis(self):
        # Add a record of the analysis to the logging file
        raise NotImplementedError

    def publish_qc(self,projects=None,location=None):
        # Copy the QC reports to the webserver
        #
        # projects: specify a pattern to match one or more projects to
        #           publish the reports for (default is to publish all reports)
        # location: override the target location specified in the settings
        #           can be of the form '[[user@]server:]directory' 
        #
        # Process pattern matching
        if projects is None:
            project_pattern = '*'
        else:
            project_pattern = projects
        # Get location to publish qc reports to
        if location is None:
            user = auto_process_settings.qc_web_server.user
            server = auto_process_settings.qc_web_server.server
            dirn = auto_process_settings.qc_web_server.webdir
        else:
            user,server,dirn = auto_process_utils.split_user_host_dir(location)
        # Check the settings
        if dirn is None:
            raise Exception, "No target directory specified"
        dirn = os.path.join(dirn,os.path.basename(self.analysis_dir))
        # Get project data
        projects = self.get_analysis_projects(project_pattern)
        # Make an index page
        title = "QC reports for %s" % os.path.basename(self.params.analysis_dir)
        index_page = qcreporter.HTMLPageWriter(title)
        # Add CSS rules
        index_page.addCSSRule("h1 { background-color: #42AEC2;\n"
                              "     color: white;\n"
                              "     padding: 5px 10px; }")
        index_page.addCSSRule("table { margin: 10 10;\n"
                              "        border: solid 1px grey;\n"
                              "        background-color: white; }")
        index_page.addCSSRule("th    { background-color: grey;\n"
                              "        color: white;\n"
                              "        padding: 2px 5px; }")
        index_page.addCSSRule("td    { text-align: right; \n"
                              "        padding: 2px 5px;\n"
                              "        border-bottom: solid 1px lightgray; }")
        index_page.addCSSRule("p.footer { font-style: italic;\n"
                              "           font-size: 70%; }")
        # Build the page
        index_page.add("<h1>%s</h1>" % title)
        index_page.add("<table>")
        index_page.add("<tr><th>Project</th><th>User</th><th>Library</th><th>Organism</th><th>PI</th><th>Samples</th><th colspan='2'>Reports</th></tr>")
        # Make a directory for the QC reports
        if server is None:
            # Local directory
            bcf_utils.mkdir(dirn)
        else:
            # Remote directory
            try:
                mkdir_cmd = applications.general.ssh_command(user,server,('mkdir',dirn))
                mkdir_cmd.run_subprocess()
            except Exception, ex:
                raise Exception, "Exception making remote directory for QC reports: %s" % ex
        # Deal with QC for each project
        for project in projects:
            # Initial information
            index_page.add("<tr>")
            index_page.add("<td>%s</td>" % project.name)
            index_page.add("<td>%s</td>" % project.metadata.user)
            index_page.add("<td>%s</td>" % project.metadata.library_type)
            index_page.add("<td>%s</td>" % project.metadata.organism)
            index_page.add("<td>%s</td>" % project.metadata.PI)
            index_page.add("<td>%s</td>" % project.prettyPrintSamples())
            # Copy QC report, if QC was verified
            if project.verify_qc():
                # Copy the qc zip file
                qc_zip = os.path.join(project.dirn,project.qc_report)
                report_copied = True
                if server is None:
                    # Local directory
                    shutil.copy(qc_zip,dirn)
                    # Unpack
                    unzip_cmd = applications.Command('unzip','-q','-o','-d',dirn,
                                                     os.path.join(dirn,project.qc_report))
                    print "Running %s" % unzip_cmd
                    unzip_cmd.run_subprocess()
                else:
                    try:
                        # Remote directory
                        scp = applications.general.scp(user,server,qc_zip,dirn)
                        print "Running %s" % scp
                        scp.run_subprocess()
                        # Unpack at the other end
                        unzip_cmd = applications.general.ssh_command(
                            user,server,
                            ('unzip','-q','-o','-d',dirn,
                             os.path.join(dirn,project.qc_report)))
                        print "Running %s" % unzip_cmd
                        unzip_cmd.run_subprocess()
                    except Exception, ex:
                        print "Failed to copy QC report: %s" % ex
                        index_page.add("<tr><td>%s</td><td>Missing</td><td>Missing</td></tr>"
                                       % project.name)
                        report_copied = False
                # Append info to the index page
                if report_copied:
                    index_page.add("<td><a href='%s/qc_report.html'>[Report]</a></td>"
                                   % os.path.splitext(project.qc_report)[0])
                    index_page.add("<td><a href='%s'>[Zip]</a></td>"
                                   % project.qc_report)
            else:
                # QC did not run
                index_page.add("<td colspan='2'>QC reports not available</td>")
            # Finish table row for this project
            index_page.add("</tr>")
        # Finish index page
        index_page.add("</table>")
        index_page.add("<p class='footer'>Generated by auto_process.py %s on %s</p>" % \
                       (__version__,time.asctime()))
        # Copy to server
        index_html = os.path.join(self.tmp_dir,'index.html')
        index_page.write(index_html)
        if server is None:
            # Local directory
            shutil.copy(index_html,dirn)
        else:
            # Remote directory
            scp = applications.general.scp(user,server,index_html,dirn)
            print "Running %s" % scp
            scp.run_subprocess()

    def report(self,logging=False,summary=False,full=False):
        # Report the contents of the run in various formats
        # Short form "logging"-style report
        if logging:
            report = self.report_logging_format()
            print report
        if summary:
            report = self.report_summary_format()
            print report
        if full:
            report = self.report_full_format()
            print report
        if not (logging or summary or full):
            illumina_data = IlluminaData.IlluminaData(self.params.analysis_dir,
                                                      unaligned_dir=self.params.unaligned_dir)
            for project in illumina_data.projects:
                print "%s" % analyse_illumina_run.describe_project(project)

    def report_logging_format(self):
        # Generate short form "logging"-style report
        # e.g. Paired end: 'PJB': Peter Briggs, Mouse ChIP-seq (PI: P Briggs) (6 samples); ...
        #
        # Acquire data
        illumina_data = IlluminaData.IlluminaData(self.params.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        project_metadata = ProjectMetadataFile(os.path.join(self.analysis_dir,
                                                            self.params.project_metadata))
        # Generate report text
        report = []
        for p in project_metadata:
            project = illumina_data.get_project(p['Project'])
            report.append("'%s': %s, %s %s (PI: %s) (%d samples)" % \
                          (p['Project'],
                           p['User'],
                           p['Organism'],
                           p['Library'],
                           p['PI'],
                           len(project.samples)))
        report = '; '.join(report)
        # Paired end run?
        if illumina_data.paired_end:
            report = "Paired end: " + report
        return report

    def report_summary_format(self):
        # Generate summary form "email"-style report for record-keeping
        # Includes:
        # - Platform
        # - Run name
        # - Project subdirectory
        # - Researcher (aka user)
        # - PI
        # - Application (aka library type)
        # - Organism
        # - Number of samples
        #
        # Acquire data
        illumina_data = IlluminaData.IlluminaData(self.params.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        project_metadata = ProjectMetadataFile(os.path.join(self.analysis_dir,
                                                            self.params.project_metadata))
        # Gather information
        datestamp = None
        instrument = None
        run_number = None
        if self.params.data_dir is not None:
            run_name = os.path.basename(self.params.data_dir)
            try:
                datestamp,instrument,run_number = IlluminaData.split_run_name(run_name)
            except Exception, ex:
                logging.warning("Unable to extract information from run name '%s'" \
                                % run_name)
                logging.warning("Exception: %s" % ex)
        else:
            run_name = os.path.basename(self.analysis_dir)
            if run_name.endswith('_analysis'):
                # Strip trailing _analysis
                run_name = run_name[:-len('_analysis')]
        if self.params.platform is not None:
            platform = self.params.platform.upper()
        else:
            platform = 'unknown'
        # Generate report text
        report = []
        if datestamp and instrument and run_number:
            report.append("%s run #%s datestamped %s\n" % (self.params.platform.upper(),
                                                           int(run_number),
                                                           datestamp))
        else:
            report.append("%s\n" % os.path.basename(self.analysis_dir))
        report.append("Run name : %s" % run_name)
        report.append("Platform : %s" % platform)
        report.append("Directory: %s" % self.params.analysis_dir)
        report.append("")
        n_projects = len(project_metadata)
        report.append("%d project%s:" % (n_projects,
                                         '' if n_projects == 1 else 's'))
        for p in project_metadata:
            project = illumina_data.get_project(p['Project'])
            report.append("- '%s':\t%s\t(PI %s)\t%s\t(%s)\t%d samples" % \
                          (p['Project'],
                           p['User'],
                           p['PI'] if p['PI'] != '?' else 'unknown',
                           p['Library'],
                           p['Organism'] if p['Organism'] != '?' else 'unknown organism',
                           len(project.samples)))
        report = '\n'.join(report)
        return report

    def report_full_format(self):
        # Generate long form "full"-style report suitable for sending
        # to bioinformaticians
        # e.g. Peter Briggs data is now available at
        #
        # /path/to/data/140204_SQ12345_0001_AB12CDXYZ/PJB/
        #
        # The samples are:
        #
        # PJB1, PJBA1-4 (5 paired end samples, multiple fastqs per sample)
        #
        # Additionally:
        # Some extra information.
        #
        # Acquire data
        illumina_data = IlluminaData.IlluminaData(self.params.analysis_dir,
                                                  unaligned_dir=self.params.unaligned_dir)
        project_metadata = ProjectMetadataFile(os.path.join(self.analysis_dir,
                                                            self.params.project_metadata))
        # Generate report text
        report = []
        for p in project_metadata:
            project = auto_process_utils.AnalysisProject(p['Project'],
                                                         os.path.join(self.params.analysis_dir,
                                                                      p['Project']))
            title = "%s %s %s data from %s run %s" % \
                          (project.metadata.user,
                           project.metadata.library_type,
                           project.metadata.organism,
                           project.metadata.platform.upper(),
                           os.path.basename(self.params.analysis_dir).split('_')[0])
            report.append("%s\n%s\n" % (title,'-'*len(title)))
            report.append("The data for %(user)s's %(org)s %(lib)s is now "
                          "available at\n\n%(dirn)s\n" % \
                          dict(user=project.metadata.user,
                               dirn=project.dirn,
                               org=project.metadata.organism,
                               lib=project.metadata.library_type))
            report.append("The samples are:\n\n%s (%d%s samples%s)" % \
                          (project.prettyPrintSamples(),
                           len(project.samples),
                           " paired end" if project.metadata.paired_end else '',
                           ", multiple fastqs per sample" if project.multiple_fastqs else ''))
            report.append("\nAdditional comments:\n\t%s" % project.metadata.comments)
        report = '\n'.join(report)
        return report

class ProjectMetadataFile(TabFile.TabFile):
    def __init__(self,filen=None):
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
        # Add project info to the metadata file
        self.append(data=(project_name,
                          bcf_utils.pretty_print_names(sample_names),
                          '.' if user is None else user,
                          '.' if library_type is None else library_type,
                          '.' if organism is None else organism,
                          '.' if PI is None else PI,
                          '.' if comments is None else comments))

    def save(self,filen=None):
        # Save the data back to file
        if filen is not None:
            self.__filen = filen
        self.write(filen=self.__filen,include_header=True)

#######################################################################
# Functions
#######################################################################

def make_custom_sample_sheet(input_sample_sheet,output_sample_sheet=None):
    # Read sample sheet info from input_sample_sheet
    # Do clean up
    # Write to output_sample_sheet (if specified)
    # Return CasavaSampleSheet object
    sample_sheet = IlluminaData.get_casava_sample_sheet(input_sample_sheet)
    for line in sample_sheet:
        if not line['SampleProject']:
            line['SampleProject'] = line['SampleID']
    sample_sheet.fix_illegal_names()
    sample_sheet.fix_duplicated_names()
    if output_sample_sheet is not None:
        sample_sheet.write(output_sample_sheet)
    return sample_sheet

def get_bases_mask(run_info_xml,sample_sheet_file):
    # Return bases mask string generated from data in RunInfo.xml and
    # sample sheet files
    # Get initial bases mask
    bases_mask = IlluminaData.IlluminaRunInfo(run_info_xml).bases_mask
    print "Bases mask: %s (from RunInfo.xml)" % bases_mask
    # Update bases mask from sample sheet
    example_barcode = IlluminaData.get_casava_sample_sheet(sample_sheet_file)[0]['Index']
    bases_mask = IlluminaData.fix_bases_mask(bases_mask,example_barcode)
    print "Bases mask: %s (updated for barcode sequence '%s')" % (bases_mask,
                                                                  example_barcode)
    return bases_mask

def list_available_commands(cmds):
    # Pretty-print available commands
        print ""
        print "Available commands are:"
        for cmd in cmds:
            print "\t%s" % cmd
        print ""

def set_debug(flag):
    # Turn on debug output
    if flag:
        logging.getLogger().setLevel(logging.DEBUG)

# Command line parsers

def setup_parser():
    p = optparse.OptionParser(usage="%prog setup [OPTIONS] DIR",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequencing "
                              "data from DIR.")
    p.add_option('--analysis-dir',action='store',dest='analysis_dir',default=None,
                 help="Make new directory called ANALYSIS_DIR (otherwise default is "
                 "'DIR_analysis')")
    p.add_option('--fastq-dir',action='store',dest='fastq_dir',default=None,
                 help="Import fastq.gz files from FASTQ_DIR (which should be a "
                 "subdirectory of DIR with the same structure as that produced "
                 "by CASAVA/bcl2fastq i.e. 'Project_<name>/Sample_<name>/<fastq>')")
    return p

def make_fastqs_parser():
    p = optparse.OptionParser(usage="%prog make_fastqs [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequence from "
                              "ANALYSIS_DIR.")
    p.add_option('--skip-rsync',action='store_true',
                 dest='skip_rsync',default=False,
                 help="don't rsync the primary data at the beginning of processing")
    p.add_option('--ignore-missing-bcl',action='store_true',
                 dest='ignore_missing_bcl',default=False,
                 help="Use the --ignore-missing-bcl option for bcl2fastq (treat "
                 "missing bcl files as no call)")
    p.add_option('--ignore-missing-stats',action='store_true',
                 dest='ignore_missing_stats',default=False,
                 help="Use the --ignore-missing-stats option for bcl2fastq (fill in "
                 "with zeroes when *.stats files are missing)")
    p.add_option('--keep-primary-data',action='store_true',
                 dest='keep_primary_data',default=False,
                 help="Don't delete the primary data at the end of processing")
    p.add_option('--generate-stats',action='store_true',
                 dest='generate_stats',default=False,
                 help="(Re)generate statistics for fastq files")
    return p

def run_qc_parser():
    p = optparse.OptionParser(usage="%prog run_qc [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Automatically process Illumina sequence from "
                              "ANALYSIS_DIR.")
    p.add_option('--projects',action='store',
                 dest='project_pattern',default=None,
                 help="simple wildcard-based pattern specifying a subset of projects "
                 "and samples to run the QC on. PROJECT_PATTERN should be of the form "
                 "'pname[/sname]', where 'pname' specifies a project (or set of "
                 "projects) and 'sname' optionally specifies a sample (or set of "
                 "samples).")
    p.add_option('--max-jobs',action='store',
                 dest='max_jobs',default=4,type='int',
                 help="explicitly specify maximum number of concurrent QC jobs to run "
                 "(default: 4)")
    return p

def publish_qc_parser():
    p = optparse.OptionParser(usage="%prog publish_qc [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Copy QC reports from ANALYSIS_DIR to local "
                              "or remote directory (e.g. web server).")
    p.add_option('--projects',action='store',
                 dest='project_pattern',default=None,
                 help="simple wildcard-based pattern specifying a subset of projects "
                 "and samples to publish the QC for. PROJECT_PATTERN can specify a "
                 "single project, or a set of projects.")
    p.add_option('--target',action='store',
                 dest='target_dir',default=None,
                 help="specify target directory to copy QC reports to. TARGET_DIR can "
                 "be a local directory, or a remote location in the form "
                 "'[[user@]host:]directory'. Overrides the default settings.")
    return p

def archive_parser():
    p = optparse.OptionParser(usage="%prog publish [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Copy sequencing analysis data directory "
                              "ANALYSIS_DIR to 'archive' destination.")
    p.add_option('--platform',action='store',
                 dest='platform',default=None,
                 help="specify the platform e.g. 'hiseq', 'miseq' etc (overrides "
                 "automatically determined platform, if any). Use 'other' for cases "
                 "where the platform is unknown.")
    p.add_option('--year',action='store',
                 dest='year',default=None,
                 help="specify the year e.g. '2014' (default is the current year)")
    return p

def report_parser():
    p  = optparse.OptionParser(usage="%prog report [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description="Report information on processed Illumina "
                               "sequence data in ANALYSIS_DIR.")
    p.add_option('--logging',action='store_true',dest='logging',default=False,
                 help="print short report suitable for logging file")
    p.add_option('--summary',action='store_true',dest='summary',default=False,
                 help="print full report suitable for bioinformaticians")
    p.add_option('--full',action='store_true',dest='full',default=False,
                 help="print summary report suitable for record-keeping")
    return p

def generic_parser(description=None):
    if description is None:
        description = "Automatically process Illumina sequence from ANALYSIS_DIR."
    p  = optparse.OptionParser(usage="%prog CMD [OPTIONS] [ANALYSIS_DIR]",
                              version="%prog "+__version__,
                              description=description)
    return p

def add_debug_option(p):
    # Add debug option to a parser
    p.add_option('--debug',action='store_true',dest='debug',default=False,
                 help="Turn on debugging output from Python libraries")
    return p

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Available commands and corresponding
    cmd_parsers = bcf_utils.OrderedDictionary()
    cmd_parsers['setup'] = setup_parser()
    cmd_parsers['make_fastqs'] = make_fastqs_parser()
    cmd_parsers['setup_analysis_dirs'] = generic_parser()
    cmd_parsers['run_qc'] = run_qc_parser()
    cmd_parsers['archive'] = archive_parser()
    cmd_parsers['publish_qc'] = publish_qc_parser()
    cmd_parsers['report'] = report_parser()

    # Process major command
    try:
        cmd = sys.argv[1]
    except IndexError:
        cmd = "help"
    if cmd == "help" or cmd == "--help" or cmd == "-h":
        list_available_commands(cmd_parsers)
        sys.exit(0)
    else:
        err = None
        if cmd not in cmd_parsers:
            err = "Unrecognised command '%s'" % cmd
        elif len(sys.argv) < 2:
            err = "Need to supply a command"
        if err is not None:
            print err
            list_available_commands(cmd_parsers)
            sys.stderr.write("%s\n" % err)
            sys.exit(1)
        p = cmd_parsers[cmd]

    # Add debug option (available for all commands)
    add_debug_option(p)
    
    # Process remaining command line arguments
    options,args = p.parse_args(sys.argv[2:])

    # Report name and version
    print "%s version %s" % (os.path.basename(sys.argv[0]),__version__)

    # Turn on debugging?
    set_debug(options.debug)

    # Setup the processing object and run the requested command
    if cmd == 'setup':
        if len(args) != 1:
            sys.stderr.write("Need to supply a data source location\n")
            sys.exit(1)
        d = AutoProcess()
        if options.fastq_dir is None:
            d.setup(args[0],analysis_dir=options.analysis_dir)
        else:
            d.setup_from_fastq_dir(args[0],options.fastq_dir)
    else:
        # For other options check if an analysis
        # directory was specified
        if len(args) > 0:
            analysis_dir = args[0]
        else:
            analysis_dir = os.getcwd()
        d = AutoProcess(analysis_dir)
        # Run the specified stage
        if cmd == 'make_fastqs':
            d.bcl_to_fastq(skip_rsync=options.skip_rsync,
                           keep_primary_data=options.keep_primary_data,
                           ignore_missing_bcl=options.ignore_missing_bcl,
                           ignore_missing_stats=options.ignore_missing_stats,
                           generate_stats=options.generate_stats)
        elif cmd == 'setup_analysis_dirs':
            d.setup_analysis_dirs()
        elif cmd == 'run_qc':
            d.run_qc(projects=options.project_pattern,
                     max_jobs=options.max_jobs)
        elif cmd == 'archive':
            d.copy_to_archive(platform=options.platform,
                              year=options.year)
        elif cmd == 'publish_qc':
            d.publish_qc(projects=options.project_pattern,
                         location=options.target_dir)
        elif cmd == 'report':
            d.report(logging=options.logging,
                     summary=options.summary,
                     full=options.full)
