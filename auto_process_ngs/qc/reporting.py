#!/usr/bin/env python
#
#     reporting: report QC from analysis projects
#     Copyright (C) University of Manchester 2018 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import sys
import os
import logging
import time
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import cmp_sample_names
from bcftbx.TabFile import TabFile
from bcftbx.qc.report import strip_ngs_extensions
from bcftbx.utils import AttributeDictionary
from bcftbx.utils import extract_prefix
from bcftbx.utils import extract_index
from ..applications import Command
from ..docwriter import Document
from ..docwriter import Section
from ..docwriter import Table
from ..docwriter import Img
from ..docwriter import Link
from ..docwriter import Target
from .fastqc import Fastqc
from .fastq_screen import Fastqscreen
from .illumina_qc import IlluminaQC
from .illumina_qc import fastqc_output
from .illumina_qc import fastq_screen_output
from .plots import uscreenplot
from .plots import ufastqcplot
from .plots import uboxplot
from .plots import encode_png
from .. import get_version

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Data
#######################################################################

from .illumina_qc import FASTQ_SCREENS

QC_REPORT_CSS_STYLES = """/* Headers */
h1 { background-color: #42AEC2;
     color: white;\n
     padding: 5px 10px; }
h2 { background-color: #8CC63F;
     color: white;
     display: inline-block;
     padding: 5px 15px;
     margin: 0;
     border-top-left-radius: 20px;
     border-bottom-right-radius: 20px; }
 h3, h4 { background-color: grey;
          color: white;
          display: block;
          padding: 5px 15px;
          margin: 0;
          border-top-left-radius: 20px;
          border-bottom-right-radius: 20px; }
/* Samples and Fastqs */
.sample { margin: 10 10;
          border: solid 2px #8CC63F;
          padding: 0;
          background-color: #ffe;
          border-top-left-radius: 25px;
          border-bottom-right-radius: 25px; }
.fastqs { border: 1px solid grey;
          padding: 5px;
          margin: 5px 20px; }
.fastq { border: 2px solid lightgray;
         padding: 5px;
         margin: 5px;
         float: left; }
.clear { clear: both; }
/* Metadata table */
table.metadata {
          margin: 10 10;
          border: solid 1px grey;
          background-color: white;
         font-size: 90%; }
table.metadata tr td:first-child {
          background-color: grey;
          color: white;
          padding: 2px 5px;
          font-weight: bold; }
/* Summary table */
table.summary { border: solid 1px grey;
          background-color: white;
          font-size: 80%; }
table.summary th { background-color: grey;
                   color: white;
                   padding: 2px 5px; }
table.summary td { text-align: right;
                   padding: 2px 5px;
                   border-bottom: solid 1px lightgray; }
table.fastq_summary tr td:first-child {
          background-color: grey;
          color: white;
          font-weight: bold; }
table.fastq_summary tr td:first-child a {
          color: white;
          font-weight: bold; }
/* FastQC summary table */
table.fastqc_summary span.PASS { font-weight: bold;
                                 color: green; }
table.fastqc_summary span.WARN { font-weight: bold;
                                 color: orange; }
table.fastqc_summary span.FAIL { font-weight: bold;
                                 color: red; }
/* Program versions */
table.programs th { text-align: left;
                    background-color: grey;
                    color: white;
                    padding: 2px 5px; }
table.programs td { padding: 2px 5px;
                    border-bottom: solid 1px lightgray; }
/* General styles */
p { font-size: 85%;
    color: #808080; }
/* Rules for printing */
@media print
{
a { color: black; text-decoration: none; }
.sample { page-break-before: always; }
table th { border-bottom: solid 1px lightgray; }
.no_print { display: none; }
}
"""

#######################################################################
# Classes
#######################################################################

class QCReporter(object):
    """
    Class describing QC results for an AnalysisProject

    Provides the follow properties:

    name: project name
    paired_end: True if project is paired-end
    samples: list of QCSample instances

    Provides the follow methods:

    verify: checks the QC outputs for the project
    report: generate a HTML report for the project
    """
    def __init__(self,project):
        """
        Initialise a new QCReporter instance

        Arguments:
           project (AnalysisProject): project to report QC for
        """
        self._project = project
        self._samples = []
        self._parent_dir = os.path.dirname(self._project.dirn)
        for sample in self._project.samples:
            self._samples.append(QCSample(sample))
        self._samples = sorted(self._samples,
                               cmp=lambda x,y: cmp_sample_names(x.name,y.name))
        logger.debug("Found %d samples" % len(self._samples))

    @property
    def name(self):
        return self._project.name

    @property
    def paired_end(self):
        return self._project.info.paired_end

    @property
    def samples(self):
        return self._samples

    def verify(self,qc_dir=None,illumina_qc=None):
        """
        Check the QC outputs are correct for the project

        Arguments:
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project being checked.
          illumina_qc (IlluminaQC): configured IlluminaQC
            instance to use for verification (optional)

        Returns:
          Boolean: Returns True if all expected QC products
            are present, False if not.
        """
        logger.debug("QCReporter.verify: qc_dir (initial): %s" % qc_dir)
        if qc_dir is None:
            qc_dir = self._project.qc_dir
        else:
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(self._project.dirn,
                                      qc_dir)
        logger.debug("QCReporter.verify: qc_dir (final)  : %s" % qc_dir)
        verified = True
        for sample in self._samples:
            if not sample.verify(qc_dir,illumina_qc=illumina_qc):
                verified = False
        return verified

    def report(self,title=None,filename=None,qc_dir=None,
               relative_links=False):
        """
        Report the QC for the project

        Arguments:
          title (str): optional, specify title for the report
            (defaults to '<PROJECT_NAME>: QC report')
          filename (str): optional, specify path and name for
            the output report file (defaults to
            '<PROJECT_NAME>.qc_report.html')
          qc_dir (str): path to the QC output dir
          relative_links (boolean): optional, if set to True
            then use relative paths for links in the report
            (default is to use absolute paths)

        Returns:
          String: filename of the output HTML report.
        """
        # Set title and output destination
        if title is None:
            title = "%s: QC report" % self.name
        if filename is None:
            filename = "%s.qc_report.html" % self.name
        # Use relative paths for links
        if relative_links:
            relpath = os.path.dirname(filename)
        else:
            relpath = None
        # Initialise report
        report = QCReport(self._project,title=title,qc_dir=qc_dir)
        # Styles
        report.add_css_rule(QC_REPORT_CSS_STYLES)
        # Write the report
        report.write(filename)
        # Return the output filename
        return filename

class QCSample(object):
    """
    Class describing QC results for an AnalysisSample

    Provides the follow properties:

    name: sample name
    fastq_pairs: list of FastqSet instances

    Provides the follow methods:

    verify: checks that QC outputs are present
    """
    def __init__(self,sample):
        """
        Initialise a new QCSample instance

        Arguments:
           sample (AnalysisSample): sample instance
        """
        self._sample = sample
        self._fastq_pairs = get_fastq_pairs(sample)

    @property
    def name(self):
        return self._sample.name

    @property
    def fastq_pairs(self):
        return self._fastq_pairs

    def verify(self,qc_dir,illumina_qc=None):
        """
        Check QC products for this sample

        Checks that expected QC outputs for the sample
        are present in the specified QC directory.

        Arguments:
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project being checked.
          illumina_qc (IlluminaQC): configured IlluminaQC
            instance to use for verification (optional)

        Returns:
          Boolean: returns True if the QC products are
            present, False otherwise.
        """
        logger.debug("QCSample.verify: qc_dir: %s" % qc_dir)
        for fq_pair in self.fastq_pairs:
            if not fq_pair.verify(qc_dir,illumina_qc=illumina_qc):
                return False
        return True

class FastqSet(object):
    """
    Class describing a set of Fastq files

    A set can be a single or a pair of fastq files.

    Provides the following properties:

    r1: R1 Fastq in the pair
    r2: R2 Fastq (will be None if no R2)
    fastqs: list of Fastq files

    Provides the following methods:

    verify: checks the QC outputs for the set
    """
    def __init__(self,fqr1,fqr2=None):
        """
        Initialise a new QCFastqSet instance

        Arguments:
           fqr1 (str): path to R1 Fastq file
           fqr2 (str): path to R2 Fastq file, or
             None if the 'set' is a single Fastq
        """
        self._fastqs = list((fqr1,fqr2))

    def __getitem__(self,key):
        return self.fastqs[key]

    @property
    def r1(self):
        """
        Return R1 Fastq file from pair
        """
        return self._fastqs[0]

    @property
    def r2(self):
        """
        Return R2 Fastq file from pair
        """
        return self._fastqs[1]

    @property
    def fastqs(self):
        """
        Return list of Fastq files in the set
        """
        return filter(lambda fq: fq is not None,
                      self._fastqs)

    def verify(self,qc_dir,illumina_qc=None):
        """
        Check QC products for this Fastq pair

        Checks that fastq_screens and FastQC files were found.
        Returns True if the QC products are present, False
        otherwise.

        Arguments:
          qc_dir (str): path to the location of the QC
            output directory
          illumina_qc (IlluminaQC): configured IlluminaQC
            instance to use for verification (optional)

        Returns:
          Boolean: returns True if the QC products are
            present, False otherwise.
        """
        logger.debug("FastqSet.verify: fastqs: %s" % (self._fastqs,))
        logger.debug("FastqSet.verify: qc_dir: %s" % qc_dir)
        if illumina_qc is None:
            illumina_qc = IlluminaQC()
        for fq in self._fastqs:
            if fq is None:
                continue
            present,missing = illumina_qc.check_outputs(fq,qc_dir)
            if missing:
                return False
        return True

class QCReport(Document):
    """
    Create a QC report document for a project

    Example usage:

    >>> report = QCReport(project)
    >>> report.write("qc_report.html")
    """
    def __init__(self,project,title=None,qc_dir=None):
        """
        Create a new QCReport instance

        Arguments:
          project (AnalysisProject): project to report QC for
          title (str): title for the report (defaults to
            "QC report: <PROJECT_NAME")
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project
        """
        logger.debug("QCReport: qc_dir (initial): %s" % qc_dir)
        # Store project
        self.project = project
        # Sort out target QC dir
        if qc_dir is None:
            qc_dir = self.project.qc_dir
        else:
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(self.project.dirn,
                                      qc_dir)
        self.qc_dir = qc_dir
        logger.debug("QCReport: qc_dir (final): %s" % self.qc_dir)
        # Set up title
        if title is None:
            title = "QC report: %s" % self.project.name
        # Initialise superclass
        Document.__init__(self,title)
        # Initialise tables
        self.metadata_table = self._init_metadata_table()
        self.summary_table = self._init_summary_table()
        # Initialise report sections
        self.preamble = self._init_preamble_section()
        self.summary = self._init_summary_section()
        # Add data
        self.report_metadata()
        for sample in self.project.samples:
            self.report_sample(sample)

    def _init_metadata_table(self):
        """
        Internal: set up a table for project metadata

        Associated CSS class is 'metadata'
        """
        metadata_tbl = Table(('item','value',))
        metadata_tbl.no_header()
        metadata_tbl.add_css_classes('metadata')
        return metadata_tbl

    def _init_summary_table(self):
        """
        Internal: set up a table for summarising samples

        Associated CSS classes are 'summary' and 'fastq_summary'
        """
        if self.project.info.paired_end:
            fields = ('sample',
                      'fastqs',
                      'reads',
                      'fastqc_r1',
                      'boxplot_r1',
                      'screens_r1',
                      'fastqc_r2',
                      'boxplot_r2',
                      'screens_r2',)
        else:
            fields = ('sample',
                      'fastq',
                      'reads',
                      'fastqc_r1',
                      'boxplot_r1',
                      'screens_r1')
        field_descriptions = { 'sample': 'Sample',
                               'fastq' : 'Fastq',
                               'fastqs': 'Fastqs (R1/R2)',
                               'reads': '#reads',
                               'fastqc_r1': 'FastQC',
                               'boxplot_r1': 'Boxplot',
                               'screens_r1': 'Screens',
                               'fastqc_r2': 'FastQC',
                               'boxplot_r2': 'Boxplot',
                               'screens_r2': 'Screens' }
        summary_tbl = Table(fields,**field_descriptions)
        summary_tbl.add_css_classes('summary','fastq_summary')
        return summary_tbl

    def _init_preamble_section(self):
        """
        Internal: set up a "preamble" section
        """
        preamble = self.add_section()
        preamble.add("Report generated by auto_process %s on %s" %
                     (get_version(),time.asctime()))
        return preamble

    def _init_summary_section(self):
        """
        Internal: set up a summary section for the report

        Associated name is 'summary'
        """
        summary = self.add_section("Summary",name='summary')
        summary.add(self.metadata_table)
        summary.add("%d samples | %d fastqs" % (len(self.project.samples),
                                                len(self.project.fastqs)))
        summary.add(self.summary_table)
        return summary

    def report_metadata(self):
        """
        Report the project metadata

        Adds entries for the project metadata to the "metadata"
        table in the report
        """
        metadata_items = ['user','PI','library_type','organism',]
        if self.project.info.single_cell_platform is not None:
            metadata_items.insert(3,'single_cell_platform')
        metadata_titles = {
            'user': 'User',
            'PI': 'PI',
            'library_type': 'Library type',
            'single_cell_platform': 'Single cell preparation platform',
            'organism': 'Organism',
        }
        for item in metadata_items:
            if self.project.info[item]:
                self.metadata_table.add_row(
                    item=metadata_titles[item],
                    value=self.project.info[item])

    def report_sample(self,sample):
        """
        Report the QC for a sample

        Reports the QC for the sample and Fastqs to
        the summary table and appends a section with
        detailed reports to the document.

        Arguments:
          sample (AnalysisSample): sample to report
        """
        sample = QCSample(sample)
        # Create a new section for the sample
        sample_report = self.add_section(
            "Sample: %s" % sample.name,
            name="sample_%s" % sample.name)
        sample_report.add_css_classes('sample')
        # Number of fastqs
        if self.project.info.paired_end:
            sample_report.add("%d fastq R1/R2 pairs" %
                              len(sample.fastq_pairs))
        else:
            sample_report.add("%d fastqs" %
                              len(sample.fastq_pairs))
        # Report each Fastq/Fastq pair
        sample_name = sample.name
        for fq_pair in sample.fastq_pairs:
            # Report Fastq pair
            fastq_pair = QCReportFastqPair(fq_pair.r1,
                                           fq_pair.r2,
                                           self.qc_dir)
            fastq_pair.report(sample_report)
            # Add line in summary table
            if sample_name is not None:
                idx = self.summary_table.add_row(sample=Link(sample_name,
                                                             sample_report))
            else:
                idx = self.summary_table.add_row(sample="&nbsp;")
            fastq_pair.update_summary_table(self.summary_table,idx=idx)

class QCReportFastqPair(object):
    """
    Utility class for reporting the QC for a Fastq pair

    Provides the following properties:

    r1: QCReportFastq instance for R1 Fastq
    r2: QCReportFastq instance for R2 Fastq

    Provides the following methods:

    report: 
    """
    def __init__(self,fastqr1,fastqr2,qc_dir):
        """
        Create a new QCReportFastqPair

        Arguments:
          fastqr1 (str): R1 Fastq file
          fastqr2 (str): R2 Fastq file (None if 'pair' is
            single-ended)
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project
        """
        self.fastqr1 = fastqr1
        self.fastqr2 = fastqr2
        self.qc_dir = qc_dir

    @property
    def paired_end(self):
        """
        True if pair consists of R1/R2 files
        """
        return (self.fastqr2 is not None)

    @property
    def r1(self):
        """
        QCReportFastq instance for R1 Fastq
        """
        return QCReportFastq(self.fastqr1,self.qc_dir)

    @property
    def r2(self):
        """
        QCReportFastq instance for R2 Fastq (None if not paired end)
        """
        if self.fastqr2 is not None:
            return QCReportFastq(self.fastqr2,self.qc_dir)
        return None

    def report(self,sample_report,attrs=None):
        """
        Add report for Fastq pair to a document section

        Creates a new subsection in 'sample_report' for the
        Fastq pair, within which are additional subsections
        for each Fastq file.

        The following 'attributes' can be reported for each
        Fastq:

        - fastqc
        - fastq_screen
        - program versions

        By default all attributes are reported.

        Arguments:
          sample_report (Section): section to add the report
            to
          attrs (list): optional list of custom 'attributes'
            to report
        """
        # Attributes to report
        if attrs is None:
            attrs = ('fastqc','fastq_screen','program_versions')
        # Add container section for Fastq pair
        fastqs_report = sample_report.add_subsection()
        fastqs_report.add_css_classes('fastqs')
        # Create sections for individual Fastqs
        for fq in (self.r1,self.r2):
            if fq is None:
                continue
            fq_report = fastqs_report.add_subsection(fq.name,
                                                     name=fq.safe_name)
            fq_report.add_css_classes('fastq')
            # Add reports for each requested 'attribute'
            for attr in attrs:
                if attr == "fastqc":
                    # FastQC outputs
                    fq.report_fastqc(fq_report)
                elif attr == "fastq_screen":
                    # FastQScreens
                    fq.report_fastq_screens(fq_report)
                elif attr == "program_versions":
                    # Versions of programs used
                    new_section = fq.report_program_versions(fq_report)
                else:
                    raise KeyError("'%s': unrecognised reporting element "
                                   % attr)
        # Add an empty section to clear HTML floats
        clear = fastqs_report.add_subsection()
        clear.add_css_classes("clear")

    def update_summary_table(self,summary_table,idx=None,fields=None):
        """
        Add a line to a summary table reporting a Fastq pair

        Creates a new line in 'summary_table' (or updates an
        existing line) for the Fastq pair, adding content for
        each specied field.

        The following fields can be reported for each Fastq
        pair:

        - fastqs (if paired-end)
        - fastq (if single-end)
        - reads
        - boxplot_r1
        - boxplot_r2
        - fastqc_r1
        - fastqc_r2
        - screens_r1
        - screens_r2

        Arguments:
          summary_table (Table): table to add the summary to
          idx (int): if supplied then indicates which existing
            table row to update (if None then a new row is
            appended)
          fields (list): list of custom fields to report
        """
        # Fields to report
        if fields is None:
            if self.paired_end:
                fields = ('fastqs',
                          'reads',
                          'boxplot_r1','boxplot_r2',
                          'fastqc_r1','fastqc_r2',
                          'screens_r1','screens_r2')
            else:
                fields = ('fastq',
                          'reads',
                          'boxplot_r1',
                          'fastqc_r1',
                          'screens_r1')
        # Add row to summary table
        if idx is None:
            idx = summary_table.add_row()
        # Populate with data
        for field in fields:
            if field == "fastq":
                summary_table.set_value(idx,'fastq',
                                        Link(self.r1.name,
                                             "#%s" % self.r1.safe_name))
            elif field == "fastqs":
                summary_table.set_value(idx,'fastqs',
                                        "%s<br />%s" %
                                        (Link(self.r1.name,
                                              "#%s" % self.r1.safe_name),
                                         Link(self.r2.name,
                                              "#%s" % self.r2.safe_name)))
            elif field == "reads":
                summary_table.set_value(idx,'reads',
                                        pretty_print_reads(
                                            self.r1.fastqc.data.basic_statistics(
                                                'Total Sequences')))
            elif field == "boxplot_r1":
                summary_table.set_value(idx,'boxplot_r1',
                                        Img(self.r1.uboxplot(),
                                            href="#boxplot_%s" %
                                            self.r1.safe_name))
            elif field == "boxplot_r2":
                summary_table.set_value(idx,'boxplot_r2',
                                        Img(self.r2.uboxplot(),
                                            href="#boxplot_%s" %
                                            self.r2.safe_name))
            elif field == "fastqc_r1":
                summary_table.set_value(idx,'fastqc_r1',
                                        Img(self.r1.ufastqcplot(),
                                            href="#fastqc_%s" %
                                            self.r1.safe_name))
            elif field == "fastqc_r2":
                summary_table.set_value(idx,'fastqc_r2',
                                        Img(self.r2.ufastqcplot(),
                                            href="#fastqc_%s" %
                                            self.r2.safe_name))
            elif field == "screens_r1":
                summary_table.set_value(idx,'screens_r1',
                                        Img(self.r1.uscreenplot(),
                                            href="#fastq_screens_%s" %
                                            self.r1.safe_name))
            elif field == "screens_r2":
                summary_table.set_value(idx,'screens_r2',
                                        Img(self.r2.uscreenplot(),
                                            href="#fastq_screens_%s" %
                                            self.r2.safe_name))
            else:
                raise KeyError("'%s': unrecognised field for summary "
                               "table" % field)

class QCReportFastq(object):
    """
    Provides interface to QC outputs for Fastq file

    Provides the following attributes:

    name: basename of the Fastq
    path: path to the Fastq
    safe_name: name suitable for use in HTML links etc
    fastqc: Fastqc instance
    fastq_screen.names: list of FastQScreen names
    fastq_screen.SCREEN.description: description of SCREEN
    fastq_screen.SCREEN.png: associated PNG file for SCREEN
    fastq_screen.SCREEN.txt: associated TXT file for SCREEN
    fastq_screen.SCREEN.version: associated version for SCREEN
    program_versions.NAME: version of package NAME

    Provides the following methods:

    report_fastqc
    report_fastq_screens
    report_program_versions
    uboxplot
    ufastqcplot
    uscreenplot
    """
    def __init__(self,fastq,qc_dir):
        """
        Create a new QCReportFastq instance

        Arguments:
          fastq (str): path to Fastq file
          qc_dir (str): path to QC directory
        """
        # Source data
        self.name = os.path.basename(fastq)
        self.path = os.path.abspath(fastq)
        self.safe_name = strip_ngs_extensions(self.name)
        # FastQC
        try:
            self.fastqc = Fastqc(os.path.join(
                qc_dir,fastqc_output(fastq)[0]))
        except Exception as ex:
            self.fastqc = None
        # Fastqscreen
        self.fastq_screen = AttributeDictionary()
        self.fastq_screen['names'] = list()
        for name in FASTQ_SCREENS:
            png,txt = fastq_screen_output(fastq,name)
            png = os.path.join(qc_dir,png)
            txt = os.path.join(qc_dir,txt)
            self.fastq_screen['names'].append(name)
            self.fastq_screen[name] = AttributeDictionary()
            self.fastq_screen[name]["description"] = name.replace('_',' ').title()
            self.fastq_screen[name]["png"] = png
            self.fastq_screen[name]["txt"] = txt
            self.fastq_screen[name]["version"] = Fastqscreen(txt).version
        # Program versions
        self.program_versions = AttributeDictionary()
        if self.fastqc is not None:
            fastqc_version = self.fastqc.version
        else:
            fastqc_version = '?'
        self.program_versions['fastqc'] = fastqc_version
        fastq_screen_versions = list(
            set([self.fastq_screen[s].version
                 for s in self.fastq_screen.names]))
        if fastq_screen_versions:
            fastq_screen_versions = ','.join(sorted(fastq_screen_versions))
        else:
            fastq_screen_versions = '?'
        self.program_versions['fastq_screen'] = fastq_screen_versions

    def report_fastqc(self,document,relpath=None):
        """
        Report the FastQC outputs to a document

        Creates a new subsection called "FastQC" with
        a copy of the FastQC sequence quality boxplot and
        a summary table of the results from each FastQC
        module.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        fastqc_report = document.add_subsection("FastQC")
        if self.fastqc:
            # FastQC quality boxplot
            fastqc_report.add("Per base sequence quality boxplot:")
            boxplot = Img(self.fastqc.quality_boxplot(inline=True),
                          height=250,
                          width=480,
                          href=self.fastqc.summary.link_to_module(
                              'Per base sequence quality',
                              relpath=relpath),
                          name="boxplot_%s" % self.safe_name)
            fastqc_report.add(boxplot)
            # FastQC summary table
            fastqc_report.add("FastQC summary:")
            fastqc_tbl = Target("fastqc_%s" % self.safe_name)
            fastqc_report.add(fastqc_tbl,
                              self.fastqc.summary.html_table(relpath=relpath))
            if relpath:
                fastqc_html_report = os.path.relpath(self.fastqc.html_report,
                                                     relpath)
            else:
                fastqc_html_report = self.fastqc.html_report
            fastqc_report.add("%s for %s" % (Link("Full FastQC report",
                                                  fastqc_html_report),
                                             self.name))
        else:
            fastqc_report.add("!!!No FastQC data available!!!")
        return fastqc_report

    def report_fastq_screens(self,document,relpath=None):
        """
        Report the FastQScreen outputs to a document

        Creates a new subsection called "Screens" with
        copies of the screen plots for each screen and
        links to the "raw" text files.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        screens_report = document.add_subsection("Screens",
                                                 name="fastq_screens_%s" %
                                                 self.safe_name)
        raw_data = list()
        for name in self.fastq_screen.names:
            # Title for the screen
            screens_report.add(self.fastq_screen[name].description)
            # Gather data for screen
            png = self.fastq_screen[name].png
            txt = self.fastq_screen[name].txt
            if relpath:
                png_href = os.path.relpath(png,relpath)
                txt_href = os.path.relpath(txt,relpath)
            else:
                png_href = png
                txt_href = txt
            # Screen plot (PNG)
            if os.path.exists(png):
                screens_report.add(Img(encode_png(png),
                                       height=250,
                                       href=png_href))
            else:
                screens_report.add("!!!No FastqScreen plot available!!!")
            # Link to raw data (TXT)
            if os.path.exists(txt):
                raw_data.append(
                    Link(self.fastq_screen[name].description,
                         txt_href).html())
            else:
                screens_report.add("!!!No FastqScreen data available!!!")
        # Add links to raw data (TXT files)
        screens_report.add("Raw screen data: " +
                           "| ".join(raw_data))
        return screens_report

    def report_program_versions(self,document):
        """
        Report the program versions to a document

        Creates a new subsection called "Program versions"
        with a table listing the versions of the QC
        programs.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        versions = document.add_subsection("Program versions")
        programs = Table(("Program","Version"))
        programs.add_css_classes("programs","summary")
        programs.add_row(Program='fastqc',
                         Version=self.program_versions.fastqc)
        programs.add_row(Program='fastq_screen',
                         Version=self.program_versions.fastq_screen)
        versions.add(programs)
        return versions

    def uboxplot(self,inline=True):
        """
        Return a mini-sequence quality boxplot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return uboxplot(self.fastqc.data.path,inline=inline)

    def ufastqcplot(self,inline=True):
        """
        Return a mini-FastQC summary plot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return ufastqcplot(self.fastqc.summary.path,inline=inline)

    def uscreenplot(self,inline=True):
        """
        Return a mini-FastQScreen summary plot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        screen_files = list()
        for name in self.fastq_screen.names:
            screen_files.append(self.fastq_screen[name].txt)
        return uscreenplot(screen_files,inline=inline)

#######################################################################
# Functions
#######################################################################

def get_fastq_pairs(sample):
    """
    Return pairs of Fastqs for an AnalysisSample instance

    Arguments:
       sample (AnalysisSample): sample to get Fastq pairs for

    Returns:
       list: list of FastqSet instances, sorted by R1 names
    """
    pairs = []
    fastqs_r1 = sample.fastq_subset(read_number=1)
    fastqs_r2 = sample.fastq_subset(read_number=2)
    for fqr1 in fastqs_r1:
        # Split up R1 name
        logger.debug("fqr1 %s" % os.path.basename(fqr1))
        dir_path = os.path.dirname(fqr1)
        # Generate equivalent R2 file
        fqr2 = sample.fastq_attrs(fqr1)
        fqr2.read_number = 2
        fqr2 = os.path.join(dir_path,"%s%s" % (fqr2,fqr2.extension))
        logger.debug("fqr2 %s" % os.path.basename(fqr2))
        if fqr2 in fastqs_r2:
            pairs.append(FastqSet(fqr1,fqr2))
        else:
            pairs.append(FastqSet(fqr1))
    pairs = sorted(pairs,cmp=lambda x,y: cmp(x.r1,y.r1))
    return pairs

def pretty_print_reads(n):
    """
    Print the number of reads with commas at each thousand

    For example:

    >>> pretty_print_reads(10409789)
    10,409,789

    Arguments:
      n (int): number of reads

    Returns:
      String: representation with commas for every thousand.
    """
    n = str(int(n))[::-1]
    n0 = []
    while len(n) >= 3:
        n0.append(n[0:3])
        n = n[3:]
    if n: n0.append(n)
    return (','.join(n0))[::-1]
