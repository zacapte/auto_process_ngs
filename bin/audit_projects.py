#!/usr/bin/env python
#
#     audit_projects.py: summarise data for multiple sequencing projects
#     Copyright (C) University of Manchester 2015-2016 Peter Briggs
#
"""
Summarise information for multiple runs processed using 'auto_process'

"""

########################################################################
# Imports
#######################################################################

import optparse
import fnmatch
import os
import sys
import bcftbx.utils as utils
from auto_process_ngs.utils import AnalysisDir
from auto_process_ngs.utils import AnalysisProject

#######################################################################
# Classes
#######################################################################

class ProjectAuditor(object):
    """
    Class for auditing analysis projects from auto_process_ngs

    """
    def __init__(self):
        self._projects = []

    def add(self,project):
        """
        Add an AnalysisProject to the auditor

        Arguments:
          project (AnalysisProject): project to add to the
            auditor

        """
        self._projects.append(AuditProject(project))

    def PIs(self,name=None):
        """
        Return list of PI names

        Arguments:
          name (str): optional PI name to look for; can
            include wildcard characters

        """
        PIs = list(set(filter(lambda p: p is not None,
                              [p.PI for p in self._projects])))
        if name is not None:
            return filter(lambda x: fnmatch.fnmatch(x,name),PIs)
        return PIs

    def projects(self,name=None,PI=None,
                 no_PI=False,
                 exclude_undetermined=False):
        """
        Return list of projects

        By default all projects are returned; the list
        can be filtered by specifying:

        name: project name (can include wildcards)
        PI: PI name (can include wildcards)

        Optional modifiers:

        no_PI: if True then return projects with no PI
        assigned (ignored if 'PI' is specified)
        exclude_undetermined: if True then projects
        named 'undetermined' are excluded

        """
        projects = [p for p in self._projects]
        # Exclude undetermined
        if exclude_undetermined:
            projects = filter(lambda p: p.name != 'undetermined',
                              projects)
        # Filter on project name
        if name is not None:
            projects = filter(lambda p: fnmatch.fnmatch(p.name,name),
                              projects)
        # Filter on PI name
        if PI is not None:
            projects = filter(lambda p: (p.PI is not None) and (fnmatch.fnmatch(p.PI,PI)),
                              projects)
        elif no_PI:
            projects = filter(lambda p: p.PI is None,projects)
        return projects

    def disk_usage(self,PI=None):
        """
        Return disk usage (bytes) for specified criteria
        """
        if PI is not None:
            return sum([p.disk_usage
                        for p in self.projects(PI=PI)])
        # Return total usage for all projects
        return sum([p.disk_usage for p in self.projects()])

class AuditProject(object):
    """
    Class representing an AnalysisProject for auditing

    """
    def __init__(self,project):
        """
        Create a new AuditProject instance

        Arguments:
          project (AnalysisProject): analysis project to audit

        """
        self._size = None
        self._project = project
        self._parent = AnalysisDir(os.path.dirname(self._project.dirn))
        #
        self.name = self._project.name
        self.run = self._project.info.run
        self.PI = self._project.info.PI
        self.user = self._project.info.user
        self.platform = self._project.info.platform
        self.samples = len(self._project.samples)
        #
        self.run_number = self._parent.metadata.run_number

    @property
    def disk_usage(self):
        """
        Return disk_usage (in bytes) for project
        """
        if self._size is None:
            self._size = get_size(self._project.dirn)
            if self._project.fastqs_are_symlinks:
                # Actual fastq files are outside the project
                # and need to be explicitly added
                for fq in self._project.fastqs:
                    try:
                        self._size += get_size(
                            utils.Symlink(fq).resolve_target())
                    except Exception,ex:
                        raise OSError("Failed to get size for fastq file "
                                      "'%s': %s" % (fq,ex))
        return self._size

#######################################################################
# Functions
#######################################################################

def get_size(f):
    """
    Return size (in bytes) for file or directory
    
    This wraps the 'get_blocks' function and returns the
    number of blocks * 512.

    """
    return get_blocks(f)*512

def get_blocks(f):
    """
    Return number of 512-byte blocks for file or directory
    
    For a file, returns the 'st_blocks' value from the
    os.lstat() function.

    For a directory, returns the sum of all 'st_blocks' values
    for the directory contents (recursing into subdirectories
    as required).

    """
    blocks = os.lstat(f).st_blocks
    if not os.path.isfile(f):
        for dirpath,dirnames,filenames in os.walk(f):
            for d in dirnames:
                blocks += os.lstat(os.path.join(dirpath,d)).st_blocks
            for f in filenames:
                blocks += os.lstat(os.path.join(dirpath,f)).st_blocks
    return blocks

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    p = optparse.OptionParser(usage="%prog DIR [DIR...]",
                              description="Summarise the disk usage for runs that have "
                              "been processed using auto_process. The supplied DIRs are "
                              "directories holding the top-level analysis directories "
                              "corresponding to different runs. The program reports "
                              "total disk usage for projects assigned to each PI across "
                              "all DIRs.")
    p.add_option("--PI",action='store',dest="PI_name",default=None,
                 help="List data for specific PI(s) matching PI_NAME (can use glob-style "
                 "patterns)")
    p.add_option("--unassigned",action='store_true',dest="unassigned",default=False,
                 help="List data for projects where PI is not assigned")
    ##p.add_option("--report",type='choice',action='store',dest="report_type",
    ##             choices=['disk_usage',],default='disk_usage',
    ##             help="Type of data to report; must be one of 'disk_usage' (default), "
    ##             "'projects'")
    ##p.add_option("--xlsx",action='store',dest="xlxs_file",default=None,
    ##             help="Output report to Excel spreadsheet XLSX_FILE")
    opts,args = p.parse_args()
    # Collect data
    auditor = ProjectAuditor()
    for d in args:
        for dirn in utils.list_dirs(d):
            dirn = os.path.join(d,dirn)
            #print "Examining %s" % dirn
            try:
                run = AnalysisDir(dirn)
                for project in run.get_projects():
                    auditor.add(project)
            except Exception,ex:
                print "Failed to load directory '%s' as run: %s" % \
                    (dirn,ex)
    # Report unassigned projects
    if opts.unassigned:
        unassigned = auditor.projects(no_PI=True,
                                      exclude_undetermined=True)
        print "%d unassigned projects (no PI):" % len(unassigned)
        for project in unassigned:
            print "%s: %s" % (project.name,project.info.run)
            sys.exit(0)
    # Fetch PIs and sort into disk usage order
    PI_list = sorted(auditor.PIs(name=opts.PI_name),
                     key=lambda x: auditor.disk_usage(PI=x),
                     reverse=True)
    # Report if no PIs were found
    if len(PI_list) == 0:
        print "No projects assigned to PIs found"
        sys.exit(0)
    # Report PIs, projects etc
    print "Summary (PI, # of projects, total usage):"
    print "========================================="
    total_projects = 0
    total_disk_usage = 0
    for PI in PI_list:
        n_projects = len(auditor.projects(PI=PI))
        disk_usage = auditor.disk_usage(PI=PI)
        print "%s\t%d\t%s" % (PI,
                              n_projects,
                              utils.format_file_size(disk_usage))
        total_projects += n_projects
        total_disk_usage += disk_usage
    print "Total usage\t%d\t%s" % (total_projects,
                                   utils.format_file_size(total_disk_usage))
    print "\nBreakdown by PI/project:"
    print "========================"
    for PI in PI_list:
        print "%s:" % PI
        for project in auditor.projects(PI=PI):
            print "\t%s:\t%s\t%s" % (project.run,
                                     project.name,
                                     utils.format_file_size(project.disk_usage))
    undetermined = auditor.projects(name='undetermined')
    if undetermined:
        print "\nUsage for 'undetermined' reads:"
        print "==============================="
        total_size = 0
        for project in undetermined:
            print "%s\t%s" % (project.run,
                              utils.format_file_size(project.disk_usage))
            total_size += project.disk_usage
